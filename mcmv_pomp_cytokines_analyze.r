library(ggplot2)
library(reshape2)
library(dplyr)
library(scales)
library(tidyr)
library(plyr)

scientific_10_1 <- function(x) {
  parse(text=gsub("\\+","",gsub("1e", "10^", scientific_format()(x))))
}

#this line only needs to be run once - it fits model to each individual mouse's data (saves results as sim_all_init_v_0.csv)
source("mcmv_pomp_cytokines_function.r")

sim_all_init = read.csv("sim_all_init_v_0.csv")

sim_all = sim_all_init

sim_all$dat[sim_all$species%in%c("Body","Salivary Gland")]=log10(sim_all$dat[sim_all$species%in%c("Body","Salivary Gland")])
sim_all$lo[sim_all$species%in%c("Body","Salivary Gland")]=log10(sim_all$lo[sim_all$species%in%c("Body","Salivary Gland")])
sim_all$high[sim_all$species%in%c("Body","Salivary Gland")]=log10(sim_all$high[sim_all$species%in%c("Body","Salivary Gland")])
sim_all$med[sim_all$species%in%c("Body","Salivary Gland")]=log10(sim_all$med[sim_all$species%in%c("Body","Salivary Gland")])
sim_all$fit[sim_all$species%in%c("Body","Salivary Gland")]=log10(sim_all$fit[sim_all$species%in%c("Body","Salivary Gland")])

sim_all$dat[sim_all$species%in%c("Vsr","Vbr","Ib","C","Is")]=log10(sim_all$dat[sim_all$species%in%c("Vsr","Vbr","Ib","C","Is")]+1)
sim_all$lo[sim_all$species%in%c("Vsr","Vbr","Ib","C","Is")]=log10(sim_all$lo[sim_all$species%in%c("Vsr","Vbr","Ib","C","Is")]+1)
sim_all$high[sim_all$species%in%c("Vsr","Vbr","Ib","C","Is")]=log10(sim_all$high[sim_all$species%in%c("Vsr","Vbr","Ib","C","Is")]+1)
sim_all$med[sim_all$species%in%c("Vsr","Vbr","Ib","C","Is")]=log10(sim_all$med[sim_all$species%in%c("Vsr","Vbr","Ib","C","Is")]+1)
sim_all$fit[sim_all$species%in%c("Vsr","Vbr","Ib","C","Is")]=log10(sim_all$fit[sim_all$species%in%c("Vsr","Vbr","Ib","C","Is")]+1)

sim_all$species = revalue(sim_all$species,c("Salivary Gland"="Salivary\nGlands"))

ribbon = subset(sim_all[,c(1:3,5,31)],species%in%c("Body","Salivary\nGlands","IE1"))
ribbon$species = factor(ribbon$species,levels = c("Salivary\nGlands","Body","IE1"))
fit = subset(sim_all[,c(1,2,4,6,7,31)],species%in%c("Body","Salivary\nGlands","IE1"))
fit$species = factor(fit$species,levels = c("Salivary\nGlands","Body","IE1"))
fit = na.omit(melt(fit,id = c("time","species","mouse")))
fit$variable = factor(fit$variable,levels = c("dat","fit","med"))
fit$variable = revalue(fit$variable,c("dat"="Data","fit"="ODE Fit","med"="Stochastic Simulation\nMedian"))

fit_chosen = subset(fit,mouse%in%c(paste0("4D-",c(1,2,3,4,5)))&variable%in%c("Data","ODE Fit"))
fit_chosen$mouse = revalue(fit_chosen$mouse,c("4D-1"="i","4D-2"="ii","4D-3"="iii","4D-4"="iv","4D-5"="v"))

#Plot example fits for subset of mouse data
(plot = ggplot()+
    geom_line(data =subset(fit_chosen,time<=32),aes(x = time,y = value,linetype = variable,colour = species),size=1)+
    facet_grid(species~mouse,scales = "free")+theme_bw()+
    scale_x_continuous(breaks = c(0,8,16,24,32))+
    labs(x = "Days Post Infection",y = "Units",colour = "Units",fill = "Units",linetype = "Line Type",title ="A")+theme_bw()+
    scale_fill_manual(labels=c(bquote(log[10](photons/s)), bquote(log[10](photons/s)), "Cells/1000 CD8+ T Cells"),values=c("#009999","#6600CC","#33CC33"))+
    scale_colour_manual(labels=c(bquote(log[10](photons/s)), bquote(log[10](photons/s)), "Cells/1000 CD8+ T Cells"),values=c("#009999","#6600CC","#33CC33"))+
    theme(text = element_text(size=10),
          legend.title=element_text(size = 10),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),
          axis.title = element_text(size=10),
          strip.text = element_text(size = 10)))


#Look at median dynamics of fits across all compartments
sim_all_summarized = as.data.frame(
  sim_all%>%dplyr::group_by(species,time)%>%dplyr::summarize(hi = quantile(fit,0.95),
                                                             med_hi = quantile(fit,0.75),
                                                             med = quantile(fit,0.5),
                                                             med_low = quantile(fit,0.25),
                                                             low = quantile(fit,0.05)))

dat = arrange(na.omit(unique(sim_all[,c("time","dat","mouse","species")])),time)
dat_summarized = as.data.frame(dat%>%dplyr::group_by(time,species)%>%dplyr::summarize(dat = median(dat)))
dat_summarized = merge(dat_summarized,sim_all_summarized[,c("time","species","med")],all = TRUE)
dat_summarized = melt(dat_summarized,id = c("time","species"))

dat_summarized$species = revalue(dat_summarized$species,c("Vsr"="Vs","Vbr"="Vb"))
sim_all_summarized$species = revalue(sim_all_summarized$species,c("Vsr"="Vs","Vbr"="Vb"))


dat_summarized$species = revalue(dat_summarized$species,c("Salivary\nGlands"="Salivary Glands"))
sim_all_summarized$species = revalue(sim_all_summarized$species,c("Salivary\nGlands"="Salivary Glands"))

dat_summarized$species = factor(dat_summarized$species,levels = c("Salivary Glands","Body","Vs","Vb","Is","Ib","C","IE1"))
sim_all_summarized$species = factor(sim_all_summarized$species,levels = c("Salivary Glands","Body","Vs","Vb","Is","Ib","C","IE1"))

dat_summarized$variable = revalue(dat_summarized$variable,c("dat"="Data Median","med"="ODE Fit Median"))

(all_ribbons = ggplot(subset(sim_all_summarized,time<=32))+
    geom_ribbon(aes(x = time,ymin = low,ymax = hi,fill= species),alpha = 0.25)+
    facet_wrap(~species)+
    scale_x_continuous(breaks = c(0,8,16,24,32))+
    geom_ribbon(aes(x = time,ymin = med_low,ymax = med_hi,fill= species),alpha = 0.25)+facet_wrap(~species,scales = "free",ncol = 2)+
    geom_line(data = na.omit(subset(dat_summarized,time<=32)),aes(x = time,y = value,colour = species,linetype = variable),size = 1)+theme_bw()+
    theme(text = element_text(size=10),
          legend.title=element_text(size = 10),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),
          axis.title = element_text(size=10),
          strip.text = element_text(size = 10))+
    labs(x = "Days Post Infection",y = "Units",colour = "Units",fill = "Units",linetype = "Line Type",title ="B")+
    scale_fill_manual(labels=c(bquote(log[10](photons/s)), bquote(log[10](photons/s)),bquote(log[10](Virions+1)),bquote(log[10](Virions+1)),bquote(log[10](Cells+1)),bquote(log[10](Cells+1)),bquote(log[10]("?")),"Cells/1000 CD8 T Cells"),
                      values=c("#009999","#6600CC","#0066CC","#FF9933","#CC3399","#FF3333","#FF6633","#33CC33"))+
    scale_colour_manual(labels=c(bquote(log[10](photons/s)), bquote(log[10](photons/s)),bquote(log[10](Virions+1)),bquote(log[10](Virions+1)),bquote(log[10](Cells+1)),bquote(log[10](Cells+1)),bquote(log[10]("?")),"Cells/1000 CD8 T Cells"),
                        values=c("#009999","#6600CC","#0066CC","#FF9933","#CC3399","#FF3333","#FF6633","#33CC33")))

#Examine how well model fit each mouse's data
likelihood_sg = unique(sim_all_init[,c(28:30)])#normalized by dividing likelihood by number of data points for a mouse 
likelihood_sg = arrange(likelihood_sg,like_norm)

#Examine parameter values of each mouse's fits
params = unique(sim_all_init[,c(8:31)])
params2 = melt(params)


#adjust to correct for fitting transformations
params2$value[params2$variable=="d"]=params2$value[params2$variable=="d"]+0.01
params2$value[params2$variable=="w"]=params2$value[params2$variable=="w"]+1
params2$value[params2$variable=="mu"]=params2$value[params2$variable=="mu"]/1000
params2$value[params2$variable=="p2"]=params2$value[params2$variable=="p2"]*10^4+5
params2$value[params2$variable=="eta2"]=params2$value[params2$variable=="eta2"]+0.02
params2$value[params2$variable=="m"]=params2$value[params2$variable=="m"]*4

params2 = params2 %>% separate(mouse, c("group", "Mouse"), sep = "-")

(par_plot = ggplot(params2,aes(y = value,x = variable,colour = variable))+scale_y_log10(label = scientific_10_1)+geom_boxplot(outlier.shape = NA)+theme_bw()+
    theme(text = element_text(size=12),
          legend.position = "none",
          legend.title=element_text(size = 12),
          legend.text = element_text(size = 12),
          axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),
          axis.title = element_text(size=12),
          strip.text = element_text(size = 12)))

#Calculate R_0 for each mouse and look at median
params_R0 = dcast(params2, group+Mouse~ variable, value.var="value")
R_0 = params_R0%>%dplyr::group_by(group,Mouse)%>%dplyr::summarize(R0 = mu/(mu+4.54*mu)*eta*p1/(delta*c)+4.54*mu/(mu+4.54*mu)*eta2*p2/(delta*c))
median(R_0$R0)
quantile(R_0$R0,0.05)
quantile(R_0$R0,0.95)
