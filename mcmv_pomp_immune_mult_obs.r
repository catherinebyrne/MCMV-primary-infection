library(pomp)
library(reshape2)
library(dplyr)
library(fitdistrplus)
library(tidyr)
library(plyr)

#read in bioimaging data
data = read.csv("mcmv_sg_inoc_luminescence_results.csv")

#determine what bioimaging background is
background = subset(data,Treatment=="control")
#since we will fit to TOTAL ROI and not MEAN, we need to see how these relate
factor = round(data$ROI_TOTAL.Photons.s/(data$MEAN.Photons.s.cm.2.sr*data$AREA.Pixels..cm.2),2)[1]
area_sg = subset(data,Site=="Salivary Gland")$AREA.Pixels..cm.2[1]
area_b = subset(data,Site=="Body")$AREA.Pixels..cm.2[1]
#thus, we will draw background from gamma distribution following shape and rate, multiply by factor, and multiply by area of region
background_init = mean(background$MEAN.Photons.s.cm.2.sr)

sg_rho = 0.5
body_rho = 0.5

#read in flow data
data_flow = read.csv("mcmv_sg_inoc_flow_results.csv")
#will fit to IE1 CD8 T cell data
data_flow = subset(data_flow,Cell_pop=="IE1")
data_flow = subset(data_flow,Site=="Blood")
data_flow = data_flow[,c(1,3,5,6,7)]
#convert percentages to be value out of 1000
data_flow$Statistic = data_flow$Statistic*10
data_flow$Statistic = round(data_flow$Statistic)
#want to find initial amount of IE1 in absence of infection - fit all control data
background_flow = subset(data_flow,Treatment=="control")

#combine bioimaging and immune data
data2 =  data[,c("ROI_TOTAL.Photons.s","Treatment","day","Site","Mouse")]
colnames(data2) = c("Statistic","Treatment","day","Site","Mouse")
data_total = merge(data2,data_flow,all = TRUE)
mouse = data_total$Mouse
data_total = data_total%>%separate(Mouse,c("Group","Mouse_number"),sep="-")
data_total = cbind(data_total,mouse)

dat2 = subset(data_total,Treatment=="infected"&Group%in%c("4B","4D"))
dat3 = dat2
dat3$Site = revalue(dat3$Site,c("Salivary Gland"="N_Vs","Body"="N_Vb","Blood"="N_T"))
dat3$Mouse_number[dat3$Group=="4D"]=as.character(as.numeric(dat3$Mouse_number[dat3$Group=="4D"])+5)
dat3$Site = paste0(dat3$Site,"_",dat3$Mouse_number)
dat3 = dat3[,c(1,3,4)]

dat3 = spread(dat3,Site,Statistic)
colnames(dat3)[1]="time"

extra = as.data.frame(matrix(NA,nrow = nrow(dat3),ncol = 20))
colnames(extra)=c(paste0("N_Vsr_",1:10),paste0("N_Vbr_",1:10))
dat3 = cbind(dat3,extra)
dat3 = rbind(dat3,c(60,rep(NA,30),rep(0,20)))

dat4 = dat3
dat4 = melt(dat4,id = "time")
dat4 = dat4%>%separate(variable,c("N","Species"),sep="_")
dat4 = dat4[,c(1,3,4)]
colnames(dat4)=c("time","species","dat")

dat4$species = revalue(dat4$species,c("Vb"="Body","Vs"="Salivary Glands","T"="IE1"))
dat4$species <- factor(dat4$species, levels = c("Body","Salivary Glands","IE1"))
dat4 = na.omit(dat4)

model <- vectorfield(
  Csnippet("DIs =eta*Vsr-delta*Is-m1*Is*T;
            DIb = eta2*Vbr-delta*Ib-m2*10*Ib*T;
            DT = alpha*(Ib+Is)/(Ib+Is+w+1)-(d+0.01)*T;
            DVsr = p1*Is-c*Vsr-mu*4.54/V_0*Vsr+mu/V_0*Vbr;
            DVbr = p2*Ib-c*Vbr+mu*4.54/V_0*Vsr-mu/V_0*Vbr;
            DVs = DVsr;
            DVb = DVbr;"))

rproc <- Csnippet("
                  #define MAX(x, y) (((x) > (y)) ? (x) : (y))
                  double rate[14];
                  double dN[14];

                  
                  rate[0] = eta*Vsr;
                  dN[0] = rpois(rate[0]*dt);
                  rate[1] = delta;
                  rate[2] = m1*T;
                  reulermultinom(2,Is,&rate[1],dt,&dN[1]);
                  rate[3] = eta2*Vbr;
                  dN[3] = rpois(rate[3]*dt);
                  rate[4] = delta;
                  rate[5] = m2*T;
                  reulermultinom(2,Ib,&rate[4],dt,&dN[4]);
                  rate[6] = alpha*(Ib+Is)/(Ib+Is+w+1);
                  dN[6] = rpois(rate[6]*dt);
                  rate[7] = (d+0.01)*T;
                  dN[7] = rpois(rate[7]*dt);
                  rate[8] = p1*Is;
                  dN[8] = rpois(rate[8]*dt);
                  rate[9] = c;
                  rate[10] = mu*4.54/V_0;
                  reulermultinom(2,Vsr,&rate[9],dt,&dN[9]);
                  rate[11] = p2*Ib;
                  dN[11] = rpois(rate[11]*dt);
                  rate[12] = c;
                  rate[13] = mu/V_0;
                  reulermultinom(2,Vbr,&rate[12],dt,&dN[12]);
                  Is = MAX((Is+dN[0]-dN[1]-dN[2]),0);
                  Ib = MAX((Ib+dN[3]-dN[4]-dN[5]),0);
                  T  = MAX((T+dN[6]-dN[7]),0);
                  Vsr = MAX((Vsr+dN[8]-dN[9]-dN[10]+dN[13]),0);
                  Vbr = MAX((Vbr+dN[11]-dN[12]+dN[10]-dN[13]),0);
                  Vs = Vsr+background_init*factor*area_sg;
                  Vb = Vbr+background_init*factor*area_b;
                  ")

init <- Csnippet("
                 Vsr = V_0;
                 Vbr = 0;
                 Vs = V_0+background_init*factor*area_sg;
                 Vb = background_init*factor*area_b;
                 Ib = 0;
                 Is = 0;
                 T=0;
                 ")

dmeas <- Csnippet("
lik=0;
if(!ISNA(N_T_1)){
lik+=dpois(N_T_1,T*rho2+1e-6,1);
}if(!ISNA(N_T_2)){
lik+=dpois(N_T_2,T*rho2+1e-6,1);
}if(!ISNA(N_T_3)){
lik+=dpois(N_T_3,T*rho2+1e-6,1);
}if(!ISNA(N_T_4)){
lik+=dpois(N_T_4,T*rho2+1e-6,1);
}if(!ISNA(N_T_5)){
lik+=dpois(N_T_5,T*rho2+1e-6,1);
}if(!ISNA(N_T_6)){
lik+=dpois(N_T_6,T*rho2+1e-6,1);
}if(!ISNA(N_T_7)){
lik+=dpois(N_T_7,T*rho2+1e-6,1);
}if(!ISNA(N_T_8)){
lik+=dpois(N_T_8,T*rho2+1e-6,1);
}if(!ISNA(N_T_9)){
lik+=dpois(N_T_9,T*rho2+1e-6,1);
}if(!ISNA(N_T_10)){
lik+=dpois(N_T_10,T*rho2+1e-6,1);
}

if(!ISNA(N_Vs_1)){
lik+=dnorm(log(N_Vs_1),log(Vs),sg_rho,1);
}if(!ISNA(N_Vs_2)){
lik+=dnorm(log(N_Vs_2),log(Vs),sg_rho,1);
}if(!ISNA(N_Vs_3)){
lik+=dnorm(log(N_Vs_3),log(Vs),sg_rho,1);
}if(!ISNA(N_Vs_4)){
lik+=dnorm(log(N_Vs_4),log(Vs),sg_rho,1);
}if(!ISNA(N_Vs_5)){
lik+=dnorm(log(N_Vs_5),log(Vs),sg_rho,1);
}if(!ISNA(N_Vs_6)){
lik+=dnorm(log(N_Vs_6),log(Vs),sg_rho,1);
}if(!ISNA(N_Vs_7)){
lik+=dnorm(log(N_Vs_7),log(Vs),sg_rho,1);
}if(!ISNA(N_Vs_8)){
lik+=dnorm(log(N_Vs_8),log(Vs),sg_rho,1);
}if(!ISNA(N_Vs_9)){
lik+=dnorm(log(N_Vs_9),log(Vs),sg_rho,1);
}if(!ISNA(N_Vs_10)){
lik+=dnorm(log(N_Vs_10),log(Vs),sg_rho,1);
}

if(!ISNA(N_Vb_1)){
lik+=dnorm(log(N_Vb_1),log(Vb),body_rho,1);
}if(!ISNA(N_Vb_2)){
lik+=dnorm(log(N_Vb_2),log(Vb),body_rho,1);
}if(!ISNA(N_Vb_3)){
lik+=dnorm(log(N_Vb_3),log(Vb),body_rho,1);
}if(!ISNA(N_Vb_4)){
lik+=dnorm(log(N_Vb_4),log(Vb),body_rho,1);
}if(!ISNA(N_Vb_5)){
lik+=dnorm(log(N_Vb_5),log(Vb),body_rho,1);
}if(!ISNA(N_Vb_6)){
lik+=dnorm(log(N_Vb_6),log(Vb),body_rho,1);
}if(!ISNA(N_Vb_7)){
lik+=dnorm(log(N_Vb_7),log(Vb),body_rho,1);
}if(!ISNA(N_Vb_8)){
lik+=dnorm(log(N_Vb_8),log(Vb),body_rho,1);
}if(!ISNA(N_Vb_9)){
lik+=dnorm(log(N_Vb_9),log(Vb),body_rho,1);
}if(!ISNA(N_Vb_10)){
lik+=dnorm(log(N_Vb_10),log(Vb),body_rho,1);
}

if(!ISNA(N_Vbr_1)){
lik+=dnorm(log(N_Vbr_1+1),log(Vbr+1),body_rho,1);
}if(!ISNA(N_Vbr_2)){
lik+=dnorm(log(N_Vbr_2+1),log(Vbr+1),body_rho,1);
}if(!ISNA(N_Vbr_3)){
lik+=dnorm(log(N_Vbr_3+1),log(Vbr+1),body_rho,1);
}if(!ISNA(N_Vbr_4)){
lik+=dnorm(log(N_Vbr_4+1),log(Vbr+1),body_rho,1);
}if(!ISNA(N_Vbr_5)){
lik+=dnorm(log(N_Vbr_5+1),log(Vbr+1),body_rho,1);
}if(!ISNA(N_Vbr_6)){
lik+=dnorm(log(N_Vbr_6+1),log(Vbr+1),body_rho,1);
}if(!ISNA(N_Vbr_7)){
lik+=dnorm(log(N_Vbr_7+1),log(Vbr+1),body_rho,1);
}if(!ISNA(N_Vbr_8)){
lik+=dnorm(log(N_Vbr_8+1),log(Vbr+1),body_rho,1);
}if(!ISNA(N_Vbr_9)){
lik+=dnorm(log(N_Vbr_9+1),log(Vbr+1),body_rho,1);
}if(!ISNA(N_Vbr_10)){
lik+=dnorm(log(N_Vbr_10+1),log(Vbr+1),body_rho,1);
}

if(!ISNA(N_Vsr_1)){
lik+=dnorm(log(N_Vsr_1+1),log(Vsr+1),sg_rho,1);
}if(!ISNA(N_Vsr_2)){
lik+=dnorm(log(N_Vsr_2+1),log(Vsr+1),sg_rho,1);
}if(!ISNA(N_Vsr_3)){
lik+=dnorm(log(N_Vsr_3+1),log(Vsr+1),sg_rho,1);
}if(!ISNA(N_Vsr_4)){
lik+=dnorm(log(N_Vsr_4+1),log(Vsr+1),sg_rho,1);
}if(!ISNA(N_Vsr_5)){
lik+=dnorm(log(N_Vsr_5+1),log(Vsr+1),sg_rho,1);
}if(!ISNA(N_Vsr_6)){
lik+=dnorm(log(N_Vsr_6+1),log(Vsr+1),sg_rho,1);
}if(!ISNA(N_Vsr_7)){
lik+=dnorm(log(N_Vsr_7+1),log(Vsr+1),sg_rho,1);
}if(!ISNA(N_Vsr_8)){
lik+=dnorm(log(N_Vsr_8+1),log(Vsr+1),sg_rho,1);
}if(!ISNA(N_Vsr_9)){
lik+=dnorm(log(N_Vsr_9+1),log(Vsr+1),sg_rho,1);
}if(!ISNA(N_Vsr_10)){
lik+=dnorm(log(N_Vsr_10+1),log(Vsr+1),sg_rho,1);
}

lik = (give_log) ? lik : exp(lik);
")

param_names = c("p1","p2","c","delta","body_rho","sg_rho","rho2","alpha","d","m1","m2","mu","eta","eta2","w","background_init","factor","area_sg","area_b","V_0")
state_names = c("Vsr","Vs","Is","T","Vbr","Ib","Vb")
par_trans = parameter_trans(log = c("p1","p2","c","delta","alpha","w","d","eta","eta2","mu"),logit = c("m1","m2"))

mcmv <- pomp(
  data=merge(dat3,data.frame(time = 0:32),all = TRUE),
  times="time",t0=0,
  skeleton = model,
  rprocess=euler(rproc,delta.t=0.001),
  dmeasure=dmeas, 
  rinit = init,
  partrans=par_trans,
  paramnames=param_names,
  statenames=state_names
)

enames = c("d","alpha","m1","m2","mu","eta","eta2","w","p1","p2")

params_init = c(p1 =10^2,p2 = 10^2,c = 8.8,delta = 1,body_rho = body_rho,sg_rho = sg_rho,rho2 = 0.95,V_0 = 1000,
                m1 =0.25,m2 = 0.025,alpha =37,d =0.1,mu = 0.85, eta = 0.174,eta2 = 0.608,w = 1*10^7,
                background_init = background_init,factor = factor,area_sg = area_sg,area_b = area_b)

for (i in 1:10){
  ofun <- mcmv %>%
    traj_objfun(
      est=enames,
      dmeasure=dmeas,
      paramnames=param_names,
      statenames=state_names,
      params =params_init
    )
  
  #fit - optimizer searches parameter space to find parameters under which the likelihood of the data, given a trajectory of the deterministic skeleton, is maximized.
  fit <- optim(
    fn=ofun,
    par=coef(ofun,enames,transform=TRUE),
    method="Nelder-Mead",
    control=list(trace=0)
  )
  #minimized value of the negative log likelihood
  fit$value
  
  #parameters of fit
  pars_fit = coef(ofun)
  params_init = pars_fit
}

#extract fit
tdat <- mcmv %>%
  trajectory(params = pars_fit,format="data.frame")

tdat = melt(tdat[,c(2,4,7,8)],id="time")
colnames(tdat)=c("time","species","fit")

tdat$species = revalue(tdat$species,c("Vb"="Body","Vs"="Salivary Glands","T"="IE1"))
tdat$species <- factor(tdat$species, levels = c("Body","Salivary Glands","IE1"))
median_data = as.data.frame(dat4%>%dplyr::group_by(time,species)%>%dplyr::summarize(med_dat = median(dat)))
fit_ode = merge(tdat,dat4,all = TRUE)
fit_ode = merge(fit_ode,median_data,all = TRUE)
fit_ode = arrange(fit_ode,time)

#Now stochastically simulate
rmeas = Csnippet("
            N_Vs_1 = exp(rnorm(log(Vs),sg_rho));
            N_Vs_2 = exp(rnorm(log(Vs),sg_rho));
            N_Vs_3 = exp(rnorm(log(Vs),sg_rho));
            N_Vs_4 = exp(rnorm(log(Vs),sg_rho));
            N_Vs_5 = exp(rnorm(log(Vs),sg_rho));
            N_Vs_6 = exp(rnorm(log(Vs),sg_rho));
            N_Vs_7 = exp(rnorm(log(Vs),sg_rho));
            N_Vs_8 = exp(rnorm(log(Vs),sg_rho));
            N_Vs_9 = exp(rnorm(log(Vs),sg_rho));
            N_Vs_10 = exp(rnorm(log(Vs),body_rho));
            N_Vb_1 = exp(rnorm(log(Vb),body_rho));
            N_Vb_2 = exp(rnorm(log(Vb),body_rho));
            N_Vb_3 = exp(rnorm(log(Vb),body_rho));
            N_Vb_4 = exp(rnorm(log(Vb),body_rho));
            N_Vb_5 = exp(rnorm(log(Vb),body_rho));
            N_Vb_6 = exp(rnorm(log(Vb),body_rho));
            N_Vb_7 = exp(rnorm(log(Vb),body_rho));
            N_Vb_8 = exp(rnorm(log(Vb),body_rho));
            N_Vb_9 = exp(rnorm(log(Vb),body_rho));
            N_Vb_10 = exp(rnorm(log(Vb),body_rho));
            N_T_1 = rbinom(T,rho2);
            N_T_2 = rbinom(T,rho2);
            N_T_3 = rbinom(T,rho2);
            N_T_4 = rbinom(T,rho2);
            N_T_5 = rbinom(T,rho2);
            N_T_6 = rbinom(T,rho2);
            N_T_7 = rbinom(T,rho2);
            N_T_8 = rbinom(T,rho2);
            N_T_9 = rbinom(T,rho2);
            N_T_10 = rbinom(T,rho2);")

mcmv <- pomp(
  data=merge(dat3,data.frame(time = 0:32),all = TRUE),
  times="time",t0=0,
  skeleton = model,
  rprocess=euler(rproc,delta.t=0.001),
  dmeasure=dmeas, 
  rmeasure = rmeas,
  rinit = init,
  partrans=par_trans,
  paramnames=param_names,
  statenames=state_names
)

simul <- as.data.frame(simulate(mcmv,params = pars_fit,
                                nsim=100))

simul = melt(simul[,c(1:31,59)],id = c("time",".id"))
simul = simul%>%separate(variable,c("N","Species","Number"),sep="_")
simul = simul[,c(1,4,6)]
colnames(simul) = c("time","species","counts")

simul_quantile = as.data.frame(simul%>%dplyr::group_by(time,species)%>%dplyr::summarize(
  lo=quantile(counts,prob=0.05),
  med = quantile(counts,prob = 0.5),
  high = quantile(counts,prob = 0.95)))

simul_quantile$species = revalue(simul_quantile$species,c("Vb"="Body","Vs"="Salivary Glands","T"="IE1"))
simul_quantile$species <- factor(simul_quantile$species, levels = c("Body","Salivary Glands","IE1"))
simul_quantile = merge(simul_quantile,fit_ode,all = TRUE)

simul_quantile$dat[simul_quantile$species!="IE1"]=log10(simul_quantile$dat[simul_quantile$species!="IE1"])
simul_quantile$lo[simul_quantile$species!="IE1"]=log10(simul_quantile$lo[simul_quantile$species!="IE1"])
simul_quantile$high[simul_quantile$species!="IE1"]=log10(simul_quantile$high[simul_quantile$species!="IE1"])
simul_quantile$med[simul_quantile$species!="IE1"]=log10(simul_quantile$med[simul_quantile$species!="IE1"])
simul_quantile$fit[simul_quantile$species!="IE1"]=log10(simul_quantile$fit[simul_quantile$species!="IE1"])
simul_quantile$med_dat[simul_quantile$species!="IE1"]=log10(simul_quantile$med_dat[simul_quantile$species!="IE1"])

get_AIC <- function(optim_fit) {
  2 * length(optim_fit$par) + 2 * optim_fit$value
}

AIC = get_AIC(fit)

simul_quantile$AIC = AIC

simul_quantile = cbind(simul_quantile,t(pars_fit))
