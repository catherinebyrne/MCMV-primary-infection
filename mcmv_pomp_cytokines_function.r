library(pomp)
library(fitdistrplus)
library(reshape2)
library(dplyr)
library(scales)
library(tidyr)
library(plyr)

#read in bioimaging data
data = read.csv("mcmv_sg_inoc_luminescence_results.csv")
#determine what bioimaging background is
background = subset(data,Treatment=="control")
#since we will fit to TOTAL ROI and not MEAN, we need to see how these relate
factor = round(data$ROI_TOTAL.Photons.s/(data$MEAN.Photons.s.cm.2.sr*data$AREA.Pixels..cm.2),2)[1]#12.57 steradians in a sphere
area_sg = subset(data,Site=="Salivary Gland")$AREA.Pixels..cm.2[1]
area_b = subset(data,Site=="Body")$AREA.Pixels..cm.2[1]
#thus, we will draw background from gamma distribution following shape and rate, multiply by factor, and multiply by area of region
background_init = mean(background$MEAN.Photons.s.cm.2.sr)

#read in flow data
data_flow_all = read.csv("mcmv_sg_inoc_flow_results.csv")
#will fit to IE1 CD8 T cell data

data_flow_chosen = function(cell_pop,site,factor = 1,mean=TRUE,chosen_mouse=NA){
  data_flow = subset(data_flow_all,Site==site&Cell_pop%in%cell_pop&Treatment=="infected")
  if(mean==TRUE){
    result = as.data.frame(data_flow%>%dplyr::group_by(day,Cell_pop)%>%dplyr::summarize(Statistic = mean(Statistic)))
  }else{
    result = subset(data_flow,Mouse==chosen_mouse)
    result = result[,c(1,8,7)]
  }
  result = arrange(result,day)
  colnames(result)=c("time","species","Statistic")
  result$Statistic = result$Statistic*factor
  result$species = revalue(result$species,c("IE1"="N_T","NK_KLRG1"="N_N"))
  result = spread(result,species,Statistic)
  return(result)
}

data_roi_chosen = function(mean = TRUE,chosen_mouse = NA,group=NULL,site=NULL){
  if(mean==TRUE){
    data_roi = subset(data,Treatment=="infected")
    if(!is.null(group)){
      data_roi = subset(data_roi,Group%in%group)
    }
    data_roi = as.data.frame(data_roi%>%dplyr::group_by(day,Site)%>%dplyr::summarize(ROI_TOTAL.Photons.s = mean(ROI_TOTAL.Photons.s)))
  }else{
    data_roi = subset(data,Mouse==chosen_mouse)
    data_roi = data_roi[,c("day","Site","ROI_TOTAL.Photons.s")]
  }
  if(!is.null(site)){
    data_roi = subset(data_roi,Site%in%site)
  }
  data_roi$Site = revalue(data_roi$Site,c("Salivary Gland"="N_Vs","Body"="N_Vb"))
  data_roi = spread(data_roi,Site,ROI_TOTAL.Photons.s)
  colnames(data_roi)[1]="time"
  return(data_roi)
}

data_all = function(mouse_chosen,n_times){
  print(mouse_chosen)
  flow_chosen = data_flow_chosen(cell_pop=c("IE1"),site = "Blood",mean = FALSE,chosen_mouse = mouse_chosen,factor  = 10)
  roi_chosen = data_roi_chosen(mean = FALSE,chosen_mouse = mouse_chosen,group = NA)
  
  dat3 = merge(flow_chosen,roi_chosen,all = TRUE)
  dat3 = arrange(dat3,time)
  dat3 = round(dat3)
  #dat3[1,2]=1
  dat3$N_Vsr = NA
  dat3$N_Vbr = NA
  dat3 = rbind(dat3,data.frame(time = 60,N_T= NA,N_Vb = NA,N_Vs = NA,N_Vsr = 0,N_Vbr = 0 ))
  
  dat4 = dat3
  dat4 = melt(dat4,id = "time")
  dat4 = dat4%>%separate(variable,c("N","Species"),sep="_")
  dat4 = dat4[,c(1,3,4)]
  colnames(dat4)=c("time","species","dat")
  
  dat4$species = revalue(dat4$species,c("Vb"="Body","Vs"="Salivary Gland","T"="IE1"))
  dat4$species <- factor(dat4$species, levels = c("Body","Salivary Gland","IE1"))
  dat4 = na.omit(dat4)
  
  sg_rho = 0.5
  body_rho = 0.5
  
  
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
if (!ISNA(N_T)){
lik+=dpois(N_T,T*rho2+1e-6,1);
}
if(!ISNA(N_Vs)){
lik+=dnorm(log(N_Vs),log(Vs),sg_rho,1);
}
if(!ISNA(N_Vb)){
lik+=dnorm(log(N_Vb),log(Vb),body_rho,1);
}
if(!ISNA(N_Vbr)){
lik+=dnorm(log(N_Vbr+1),log(Vbr+1),body_rho,1);
}
if(!ISNA(N_Vsr)){
lik+=dnorm(log(N_Vsr+1),log(Vsr+1),sg_rho,1);
}
lik = (give_log) ? lik: exp(lik);
")
  
  model <- vectorfield(
    Csnippet("
            double fun;
              fun = (1-b1*t/(b2+t));
            DIs =eta*Vsr-delta*Is;
            DIb = (eta2+0.02)*Vbr-delta*Ib-m*4*Ib*T;
            DT = alpha*(Ib+Is)/(Ib+Is+w+1)-(d+0.01)*T;
            DVsr = p1*fun*Is-c*Vsr-mu*4.54/V_0*Vsr+mu/V_0*Vbr;
            DVbr = (p2*pow(10,4)+5)*Ib-c*Vbr+mu*4.54/V_0*Vsr-mu/V_0*Vbr;
            DVs = DVsr;
            DVb = DVbr;"))
  
  param_names = c("p1","p2","c","delta","sg_rho","body_rho","rho2","alpha","d","m","mu","eta","eta2","b1","b2","w","background_init","factor","area_sg","area_b","V_0")
  state_names = c("Vsr","Vs","Is","T","Vbr","Ib","Vb")
  par_trans = parameter_trans(log = c("p1","c","delta","alpha","w","d","eta","eta2","mu","b2"),logit = c("b1","m","p2"))
  
  rproc <- Csnippet("
                  #define MAX(x, y) (((x) > (y)) ? (x) : (y))
                  double rate[13];
                  double dN[13];
                  double fun;
                    fun = (1-b1*t/(b2+t));
                  rate[0] = eta*Vsr;
                  dN[0] = rpois(rate[0]*dt);
                  rate[1] = delta*Is;
                  dN[1] = rpois(rate[1]*dt);
                  rate[2] = (eta2+0.02)*Vbr;
                  dN[2] = rpois(rate[2]*dt);
                  rate[3] = delta;
                  rate[4] = m*2*T;
                  reulermultinom(2,Ib,&rate[3],dt,&dN[3]);
                  rate[5] = alpha*(Ib+Is)/(Ib+Is+w+1);
                  dN[5] = rpois(rate[5]*dt);
                  rate[6] = (d+0.01)*T;
                  dN[6] = rpois(rate[6]*dt);
                  rate[7] = p1*fun*Is;
                  dN[7] = rpois(rate[7]*dt);
                  rate[8] = c;
                  rate[9] = mu/V_0*4.54;
                  reulermultinom(2,Vsr,&rate[8],dt,&dN[8]);
                  rate[10] = (p2*pow(10,4)+5)*Ib;
                  dN[10] = rpois(rate[10]*dt);
                  rate[11] = c;
                  rate[12] = mu/V_0;
                  reulermultinom(2,Vbr,&rate[11],dt,&dN[11]);
                  Is = MAX((Is+dN[0]-dN[1]),0);
                  Ib = MAX((Ib+dN[2]-dN[3]-dN[4]),0);
                  T  = MAX((T+dN[5]-dN[6]),0);
                  Vsr = MAX((Vsr+dN[7]-dN[8]-dN[9]+dN[12]),0);
                  Vbr = MAX((Vbr+dN[10]-dN[11]+dN[9]-dN[12]),0);
                  Vs = Vsr+background_init*factor*area_sg;
                  Vb = Vbr+background_init*factor*area_b;
                  ")
  
  
  rmeas = Csnippet("
            N_Vs = exp(rnorm(log(Vs),sg_rho));
            N_Vb = exp(rnorm(log(Vb),body_rho));
            N_T = rbinom(T,rho2);")
  
  mcmv <- pomp(
    data=merge(data.frame(time = seq(0,max(dat3$time),1)),dat3,all = TRUE),
    times="time",t0=0,
    skeleton = model,
    dmeasure=dmeas,
    rmeasure = rmeas,
    rprocess=euler(rproc,delta.t=0.001),
    rinit = init,
    partrans=par_trans,
    paramnames=param_names,
    statenames=state_names)
  
  #Now fit model to data
  enames = c("d","alpha","m","mu","eta","eta2","b1","b2","w","p1","p2")
  
  params_init = c(p1 =5.613523e+01,p2=2.237391e+01/10^4,c = 8.8,delta = 1,sg_rho = sg_rho, body_rho = body_rho,rho2 = 0.95,V_0 = 1000,
                  m =0.25/4,alpha =8.631208e+02,d =7.984945e-02,mu = 5.906199e+02, eta = 3.056654e+00,eta2 = 4.274397e-01 ,b1 = 9.667376e-01,b2 = 3.753485e-01,w = 5.487836e+08,
                  background_init = background_init,factor = factor,area_sg = area_sg,area_b = area_b)
  
  #create objective function -  quantifies the mismatch between model predictions and data
  #ofun saves information each time it is evaluated
  
  for (i in 1:n_times){
    print(i)
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
    likelihood  = fit$value
    
    #parameters of fit
    pars_fit = coef(ofun)
    pars_fit2 = t(as.data.frame(pars_fit))
    params_init = pars_fit
  }
  
  #extract fit and plot
  tdat <- mcmv %>%
    trajectory(params = pars_fit,format="data.frame")
  
  #tdat = melt(tdat[,c(2,4,7,9)],id="time")
  tdat = melt(tdat[1:8],id="time")
  colnames(tdat)=c("time","species","fit")
  
  tdat$species = revalue(tdat$species,c("Vb"="Body","Vs"="Salivary Gland","T"="IE1"))
  fit_ode = merge(tdat,dat4,all = TRUE)
  fit_ode = arrange(fit_ode,time)
  
  #Now stochastically simulate
  
  simul <- as.data.frame(simulate(mcmv,params = pars_fit,
                                  nsim=100))
  
  simul = melt(simul[,c(1:4,7,9,11,12)],id = c("time"))
  
  colnames(simul) = c("time","species","counts")
  
  simul_quantile = as.data.frame(simul%>%dplyr::group_by(time,species)%>%dplyr::summarize(
    lo=quantile(counts,prob=0.025),
    med = quantile(counts,prob = 0.5),
    high = quantile(counts,prob = 0.975)))
  
  simul_quantile$species = revalue(simul_quantile$species,c("N_Vb"="Body","N_Vs"="Salivary Gland","N_T"="IE1"))
  simul_quantile = merge(simul_quantile,fit_ode,all = TRUE)
  simul_quantile = data.frame(simul_quantile,data.frame(pars_fit2,likelihood,like_norm = as.numeric(likelihood)/nrow(dat4),mouse = mouse_chosen))
  return(simul_quantile)
}

times = data.frame(mouse = unique(subset(data,Treatment=="infected")$Mouse),times = 2)

sim_all_init = NULL
for (i in 1:nrow(times)){
  new = data_all(as.character(times$mouse[i]),times$times[i])
  sim_all_init = rbind(sim_all_init,new)
}
#Save output to file
write.csv(sim_all_init,"sim_all_init_v_0.csv","row.names"=FALSE)
