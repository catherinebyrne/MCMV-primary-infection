library(pomp)
library(ggplot2)
library(reshape2)
library(dplyr)
library(fitdistrplus)
library(tidyr)
library(plyr)
library(tidyverse)
library(broom)
library(reshape2)


data = read.csv("mcmv_sg_inoc_luminescence_results.csv")
#determine what bioimaging background is
background = subset(data,Treatment=="control")
#since we will fit to TOTAL ROI and not MEAN, we need to see how these relate
factor = round(data$ROI_TOTAL.Photons.s/(data$MEAN.Photons.s.cm.2.sr*data$AREA.Pixels..cm.2),2)[1]
area_sg = subset(data,Site=="Salivary Gland")$AREA.Pixels..cm.2[1]
area_b = subset(data,Site=="Body")$AREA.Pixels..cm.2[1]
#thus, we will draw background from gamma distribution following shape and rate, multiply by factor, and multiply by area of region
background_init = mean(background$MEAN.Photons.s.cm.2.sr)

init <- Csnippet("
                 Vs_r = V_0;
                 Vb_r = 0;
                 Vs = V_0+background_init*factor*area_sg;
                 Vb = background_init*factor*area_b;
                 Ib = 0;
                 Is = 0;
                 C = 0;
                 T=0;
                 ")

dmeas <- Csnippet("
lik=0;
if (!ISNA(N_T)){
lik+=dpois(N_T,T*rho2+1e-6,1);
}
if(!ISNA(N_Vs)){
lik+=dnorm(log(N_Vs),log(Vs),rho1,1);
}
if(!ISNA(N_Vb)){
lik+=dnorm(log(N_Vb),log(Vb),rho1,1);
}
lik = (give_log) ? lik: exp(lik);
")

model <- vectorfield(
  Csnippet("DC = beta*Vs_r-z*C;
            DIs =eta*Vs_r-delta*Is;
            DIb = (eta2+pow(10,-3))*Vb_r-delta*Ib-m*Ib*T;
            DT = alpha*(Ib+Is)/(Ib+Is+w+1)-(d+0.01)*T;
            DVs_r = p*(exp(-y*C))*Is-c*Vs_r-mu*Vs_r+mu*Vb_r;
            DVb_r = p*Ib-c*Vb_r+mu*Vs_r-mu*Vb_r;
            DVs = DVs_r;
            DVb = DVb_r;"))

param_names = c("p","c","delta","rho1","rho2","alpha","d","m","mu","eta","eta2","z","y","beta","w","background_init","factor","area_sg","area_b","V_0")
state_names = c("Vs_r","Vs","Is","T","Vb_r","Ib","Vb","C")
par_trans = parameter_trans(log = c("p","c","delta","rho1","alpha","w","eta","eta2","z","beta","mu"),logit = c("rho2","y","d","m"))

rproc <- Csnippet("
                  #define MAX(x, y) (((x) > (y)) ? (x) : (y))
                  double rate[15];
                  double dN[15];
                  rate[0] = beta*Vs_r;
                  dN[0] = rpois(rate[0]*dt);
                  rate[1] = z*C;
                  dN[1] = rpois(rate[1]*dt);
                  rate[2] = eta*Vs_r;
                  dN[2] = rpois(rate[2]*dt);
                  rate[3] = delta*Is;
                  dN[3] = rpois(rate[3]*dt);
                  rate[4] = (eta2+pow(10,-3))*Vb_r;
                  dN[4] = rpois(rate[4]*dt);
                  rate[5] = delta;
                  rate[6] = (m+0.0001)*T;
                  reulermultinom(2,Ib,&rate[5],dt,&dN[5]);
                  rate[7] = alpha*(Ib+Is)/(Ib+Is+w+1);
                  dN[7] = rpois(rate[7]*dt);
                  rate[8] = (d+0.01)*T;
                  dN[8] = rpois(rate[8]*dt);
                  rate[9] = p*(exp(-y*C))*Is;
                  dN[9] = rpois(rate[9]*dt);
                  rate[10] = c;
                  rate[11] = mu;
                  reulermultinom(2,Vs_r,&rate[10],dt,&dN[10]);
                  rate[12] = p*Ib;
                  dN[12] = rpois(rate[12]*dt);
                  rate[13] = c;
                  rate[14] = mu;
                  reulermultinom(2,Vb_r,&rate[13],dt,&dN[13]);
                  C = MAX((C+dN[0]-dN[1]),0);
                  Is = MAX((Is+dN[2]-dN[3]),0);
                  Ib = MAX((Ib+dN[4]-dN[5]-dN[6]),0);
                  T  = MAX((T+dN[7]-dN[8]),0);
                  Vs_r = MAX((Vs_r+dN[9]-dN[10]-dN[11]+dN[14]),0);
                  Vb_r = MAX((Vb_r+dN[12]-dN[13]+dN[11]-dN[14]),0);
                  Vs = Vs_r+background_init*factor*area_sg;
                  Vb = Vb_r+background_init*factor*area_b;
                  ")

rmeas = Csnippet("
            N_Vs = exp(rnorm(log(Vs),rho1));
            N_Vb = exp(rnorm(log(Vb),rho1));
            N_T = rbinom(T,rho2);")

mcmv <- pomp(
  data=data.frame(time = seq(0,8,0.1),N_Vs = NA, N_Vb = NA,N_T = NA),
  times="time",t0=0,
  skeleton = model,
  dmeasure=dmeas,
  rmeasure = rmeas,
  rprocess=euler(rproc,delta.t=0.001),
  rinit = init,
  partrans=par_trans,
  paramnames=param_names,
  statenames=state_names)

sim_fun = function(params_chosen){
  simul <- as.data.frame(simulate(mcmv,params = params_chosen,
                                  nsim=1000))
  colnames(simul)[13]="number"
  
  simul = melt(simul[,c(1:4,5,7,9,10,12,13)],id = c("time","number"))
  
  colnames(simul) = c("time","number","species","counts")
  
  
  simul$species = revalue(simul$species,c("N_Vb"="Body","N_Vs"="Salivary Gland","N_T"="IE1"))
  simul = data.frame(simul)
  return (simul)
  
}
params= c(p =10^2,c = 8.8,delta = 1,rho1 =0.50,rho2 = 0.95,V_0 = 1000,
          m =0.633-0.0001,alpha =193,d =0.0838-0.01,mu = 0.533*1000, eta = 0.261,eta2 = 0.0574-10^(-3),
          z = 10^(-2),y = 5.32*10^(-5),beta = 2.05*10^(-3),w = 1.21*10^7-1,
          background_init = background_init,factor = factor,area_sg = area_sg,area_b = area_b)


#Simulate model for all potential viral loads
all_sim = NULL
for (i in seq(10,300,10)){
  params_chosen = params
  params_chosen[6]=i
  hi = data.frame(sim_fun(params_chosen),V_init = i)
  all_sim = rbind(all_sim,hi)
}

all_sim$number_Vinit = paste0(all_sim$number,"_",all_sim$V_init)

#Define different types of infections - based on whether viral loads increase and the site of where replication is occurring (not mutually exclusive)
#active viral replication defined as an increase in viral load
change = subset(all_sim,species%in%c("Vs_r","Vb_r"))
change = dcast(change, time+number+V_init+number_Vinit ~ species, value.var="counts")
change$total_v = change$Vs_r+change$Vb_r
change$change = c(NA,(change$total_v[2:nrow(change)]-change$total_v[1:(nrow(change)-1)]))
change$change[which(change$time==0)]=0

pos_change = unique(subset(change,change>0)$number_Vinit)

#successful infection in the body is where you simply have infected cells appear in the body at some time
frac_successful_body = unique(subset(all_sim,species=="Ib"&counts>0)$number_Vinit)
frac_successful_body = subset(all_sim,number_Vinit%in%frac_successful_body)
frac_successful_body = subset(frac_successful_body,number_Vinit%in%pos_change)

#successful infection in the salivary glands is where you simply have infected cells appear in the salivary gland at some time
frac_successful_sg = unique(subset(all_sim,species=="Is"&counts>0)$number_Vinit)
frac_successful_sg = subset(all_sim,number_Vinit%in%frac_successful_sg)
frac_successful_sg = subset(frac_successful_sg,number_Vinit%in%pos_change)

#successful infection is where you have infected cells appear in both the salivary glands and the body 
frac_successful = unique(frac_successful_sg$number_Vinit)[unique(frac_successful_sg$number_Vinit)%in%unique(frac_successful_body$number_Vinit)]
frac_successful= subset(all_sim,number_Vinit%in%frac_successful)

#sucessful long infection is where you have sucessful infection in body and salivary glands and there are infected cells in the salivary glands past 7 days post infection
frac_successful_long = subset(frac_successful,number_Vinit%in%unique(subset(frac_successful,species=="Ib"&time>7&counts>0)$number_Vinit))
frac_successful_long = unique(subset(frac_successful_long,species=="Is"&time>7&counts>0)$number_Vinit)
frac_successful_long = subset(frac_successful,number_Vinit%in%frac_successful_long)
frac_successful_long$type = "successful"

#unsuccessful infection is where you have neither infected cells in the salivary gland or the body at any time
frac_unsuccessful = subset(all_sim,!number_Vinit%in%unique(c(frac_successful_sg$number_Vinit,frac_successful_body$number_Vinit)))
frac_unsuccessful$type = "unsuccessful"

#transient infection is where you have infection somewhere, but it dies out before 7 days
frac_transient = subset(all_sim,!number_Vinit%in%unique(c(frac_unsuccessful$number_Vinit,frac_successful_long$number_Vinit)))

#transient both sites is where you have infection within both sites but it dies out before 7 days
frac_transient_both_site = subset(frac_transient,number_Vinit%in%unique(frac_successful$number_Vinit))
frac_transient_both_site$type = "transient_both_site"

#transient one site is where you have infection at only one site but it dies out before 7 days
frac_transient_one_site = subset(frac_transient,!number_Vinit%in%unique(frac_successful$number_Vinit))

#transient infection in sg is where you have infection only in the salivary glands but it dies out before 7 days
frac_transient_sg = subset(frac_transient_one_site,number_Vinit%in%unique(frac_successful_sg$number_Vinit))
frac_transient_sg$type = "transient_sg"

#transient infection in the body is where you have infection only in the body but it dies out before 7 days
frac_transient_body = subset(frac_transient_one_site,number_Vinit%in%unique(frac_successful_body$number_Vinit))
frac_transient_body$type = "transient_body"

#check that nothing has been classified into two of final categories 
check = c(unique(frac_successful_long$number_Vinit),
          unique(frac_unsuccessful$number_Vinit),
          unique(frac_transient_sg$number_Vinit),
          unique(frac_transient_body$number_Vinit),
          unique(frac_transient_both_site$number_Vinit))
length(check)
sum(duplicated(check))


#combine all data with assigned infection type
all_data = rbind(frac_successful_long,frac_unsuccessful,frac_transient_sg,frac_transient_body,frac_transient_both_site)


all_frac = NULL
for (i in unique(all_sim$V_init)){
  new = data.frame(type = c("successful","unsuccessful","transient","transient_both","transient_sg","transient_body"),
                   percent = c(length(unique(subset(frac_successful_long,V_init==i)$number_Vinit)),
                               length(unique(subset(frac_unsuccessful,V_init==i)$number_Vinit)),
                               length(unique(subset(frac_transient,V_init==i)$number_Vinit)),
                               length(unique(subset(frac_transient_both_site,V_init==i)$number_Vinit)),
                               length(unique(subset(frac_transient_sg,V_init==i)$number_Vinit)),
                               length(unique(subset(frac_transient_body,V_init==i)$number_Vinit)))/1000,
                   V_init = i)
  all_frac = rbind(all_frac,new)
}

trans_transformed = NULL

for (i in unique(all_sim$V_init)){
  successful = subset(all_frac,V_init ==i&type == "successful")$percent
  transient = subset(all_frac,V_init==i&type=="transient_sg")$percent
  new = transient/(successful+transient)
  trans_transformed = rbind(trans_transformed,data.frame(type = "trans_transformed",percent = new,V_init = i))
}

all_frac = rbind(all_frac,trans_transformed)

sum_fit = as.data.frame(all_frac%>%dplyr::filter(type%in%c("successful","trans_transformed"))%>%dplyr::group_by(type)%>%do(fit = nls(percent ~ SSasymp(V_init, yf, y0, log_alpha), data = .)) %>% 
                          tidy(fit)%>%dplyr::select(type,term,estimate)%>%
                          spread(term, estimate) %>% 
                          mutate(alpha = exp(log_alpha)))


alpha_successful= sum_fit$alpha[1]
start_successful = sum_fit$y0[1]
end_successful = sum_fit$yf[1]

alpha_frac_trans = sum_fit$alpha[2]
start_frac_trans = sum_fit$y0[2]
end_frac_trans = sum_fit$yf[2]


successful_fit = data.frame(V_init= seq(0,500,0.01),fit = end_successful+(start_successful-end_successful)*exp(-alpha_successful*seq(0,500,0.01)),variable = "successful")

id_50_val = log((0.5- end_successful)/(start_successful-end_successful))/(-alpha_successful)

(successful_frac = ggplot()+geom_point(data = subset(all_frac,type=="successful"),aes(x = V_init,y = percent),colour = "blue",alpha = 0.25)+
    geom_line(data = successful_fit,aes(x = V_init,y = fit),colour = "blue",size = 1)+
    theme_bw()+labs(x = "Initial Dose (pfu)",y="Fraction",title ="A")+
    theme(text = element_text(size=10),
          legend.position = "none",
          legend.title=element_text(size = 10),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),
          axis.title = element_text(size=10),
          strip.text = element_text(size = 10))+
    scale_x_continuous(limits = c(0,300))+
    #geom_hline(yintercept = 0.5,colour = "red",size = 1)+
    geom_label(aes(120,0.5,label=paste0("52,0.50")))+
    geom_point(aes(x = 51.56918,y = 0.5),size = 2,colour = "red")+
    scale_colour_discrete("Type of Infection",labels = c("Successful Salivary Gland\nInfection","Successful Systemic\nInfection")))

frac_trans_fit = data.frame(V_init = seq(0,500,0.01),fit = end_frac_trans+(start_frac_trans-end_frac_trans)*exp(-alpha_frac_trans*seq(0,500,0.01)),variable= "frac_trans")

frac_trans_og_fit = data.frame(V_init = seq(0,500,0.01),fit = frac_trans_fit$fit*successful_fit$fit/(1-frac_trans_fit$fit))

max_trans = frac_trans_og_fit[which(frac_trans_og_fit$fit==max(frac_trans_og_fit$fit)),]

(trans_frac = ggplot()+geom_point(data = subset(all_frac,type=="transient_sg"),aes(x = V_init,y = percent),colour = "purple",alpha = 0.25)+
    geom_line(data = frac_trans_og_fit,aes(x = V_init,y = fit),colour = "purple",size = 1)+
    theme_bw()+labs(x = "Initial Dose (pfu)",y="Fraction",title ="B")+
    theme(text = element_text(size=10),
          legend.position = "none",
          legend.title=element_text(size = 10),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),
          axis.title = element_text(size=10),
          strip.text = element_text(size = 10))+
    geom_label(aes(120,0.08881065,label=paste0("38,0.089")))+
    geom_point(aes(x = 37.57,y =0.08881065),size = 2,colour = "red")+
    scale_x_continuous(limits = c(0,300))+
    scale_colour_discrete("Type of Infection",labels = c("Successful Salivary Gland\nInfection","Successful Systemic\nInfection")))


#Repeat simulation for just the viral load at which transient infection in the salivary glands is most likely to occur (38 pfu)

all_sim = NULL
for (i in 38){
  params_chosen = params
  params_chosen[6]=i
  hi = data.frame(sim_fun(params_chosen),V_init = i)
  all_sim = rbind(all_sim,hi)
}

all_sim$number_Vinit = paste0(all_sim$number,"_",all_sim$V_init)

#active viral replication defined as an increase in viral load
change = subset(all_sim,species%in%c("Vs_r","Vb_r"))
change = dcast(change, time+number+V_init+number_Vinit ~ species, value.var="counts")
change$total_v = change$Vs_r+change$Vb_r
change = arrange(change,number,time)
change$change = c(NA,(change$total_v[2:nrow(change)]-change$total_v[1:(nrow(change)-1)]))
change$change[which(change$time==0)]=0

pos_change = unique(subset(change,change>0)$number_Vinit)

frac_successful_body = unique(subset(all_sim,species=="Ib"&counts>0)$number_Vinit)
frac_successful_body = subset(all_sim,number_Vinit%in%frac_successful_body)
frac_successful_body = subset(frac_successful_body,number_Vinit%in%pos_change)

frac_successful_sg = unique(subset(all_sim,species=="Is"&counts>0)$number_Vinit)
frac_successful_sg = subset(all_sim,number_Vinit%in%frac_successful_sg)
frac_successful_sg = subset(frac_successful_sg,number_Vinit%in%pos_change)

frac_successful = unique(frac_successful_sg$number_Vinit)[unique(frac_successful_sg$number_Vinit)%in%unique(frac_successful_body$number_Vinit)]
frac_successful= subset(all_sim,number_Vinit%in%frac_successful)


frac_successful_long = subset(frac_successful,number_Vinit%in%unique(subset(frac_successful,species=="Ib"&time>7&counts>0)$number_Vinit))
frac_successful_long = unique(subset(frac_successful_long,species=="Is"&time>7&counts>0)$number_Vinit)
frac_successful_long = subset(frac_successful,number_Vinit%in%frac_successful_long)
frac_successful_long$type = "successful"

frac_unsuccessful = subset(all_sim,!number_Vinit%in%unique(c(frac_successful_sg$number_Vinit,frac_successful_body$number_Vinit)))
frac_unsuccessful$type = "unsuccessful"

frac_transient = subset(all_sim,!number_Vinit%in%unique(c(frac_unsuccessful$number_Vinit,frac_successful_long$number_Vinit)))

frac_transient_both_site = subset(frac_transient,number_Vinit%in%unique(frac_successful$number_Vinit))
frac_transient_both_site$type = "transient_both_site"

frac_transient_one_site = subset(frac_transient,!number_Vinit%in%unique(frac_successful$number_Vinit))
frac_transient_sg = subset(frac_transient_one_site,number_Vinit%in%unique(frac_successful_sg$number_Vinit))
frac_transient_sg$type = "transient_sg"

frac_transient_body = subset(frac_transient_one_site,number_Vinit%in%unique(frac_successful_body$number_Vinit))
frac_transient_body$type = "transient_body"

check = c(unique(frac_successful_long$number_Vinit),
          unique(frac_unsuccessful$number_Vinit),
          unique(frac_transient_sg$number_Vinit),
          unique(frac_transient_body$number_Vinit),
          unique(frac_transient_both_site$number_Vinit))
length(check)
sum(duplicated(check))

#combine all data with assigned infection type
all_data = rbind(frac_successful_long,frac_unsuccessful,frac_transient_sg,frac_transient_body,frac_transient_both_site)

#plot infection dynamics of transient sg infections

max_Is= as.data.frame(frac_transient_sg%>%dplyr::filter(species=="Is")%>%dplyr::group_by(number,species)%>%dplyr::slice(which.max(counts)))

ggplot(max_Is,aes(x = species,y = counts,colour = type))+geom_boxplot()+scale_y_continuous(limits = c(0,6))

ggplot(max_Is,aes(x = species,y = time,colour = type))+geom_boxplot()+scale_y_continuous(limits = c(0,4))

#what time does Is reach its maximum?
quantile(max_Is$time,c(0.05,0.25,0.5,0.75,0.95))
#what is the maximum number of Is?
quantile(max_Is$counts,c(0.05,0.25,0.5,0.75,0.95))

Is_dynamics = as.data.frame(frac_transient_sg%>%dplyr::filter(species=="Is")%>%dplyr::group_by(time)%>%dplyr::summarize(lo = quantile(counts,0.05),
                                                                                                                        med = quantile(counts,0.5),
                                                                                                                        high = quantile(counts,0.95)))

(Is_dynamics_graph = ggplot(Is_dynamics)+geom_line(aes(x = time,y = med),colour = "pink",size = 1)+
  geom_ribbon(aes(x = time,ymin = lo,ymax = high),alpha = 0.5,fill = "pink")+
  theme_bw()+
  scale_x_continuous(limits = c(0,10),breaks = c(0,2,4,6,8,10))+
  labs(y="Infected Cells",x = "Days Post Infection",title = "C")+
  theme(text = element_text(size=10),
        legend.position = "none",
        legend.title=element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),
        axis.title = element_text(size=10),
        strip.text = element_text(size = 10)))








