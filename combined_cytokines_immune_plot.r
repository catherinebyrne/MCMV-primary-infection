library(ggplot2)
library(extrafont)
font_import()
loadfonts(device = "win")

cytokines_fun = function(){
  source("mcmv_pomp_cytokines_mult_obs.r")
return(simul_quantile)
}

no_cytokines_fun = function(){
  source("mcmv_pomp_immune_mult_obs.r")
  return(simul_quantile)
}

no_cytokines = no_cytokines_fun()
cytokines = cytokines_fun()

all_dat = rbind(data.frame(no_cytokines[,1:9],model = "Model 1"),
                data.frame(cytokines[,1:9],model = "Model 2"))

all_dat$species = revalue(all_dat$species,c("Salivary Glands"="Salivary\nGlands"))
all_dat$species = factor(all_dat$species,levels = c("Salivary\nGlands","Body","IE1"))

dat_point = na.omit(all_dat[,c(1,2,7,10)])
fit = all_dat[,c(1,2,4,6,8,10)]
fit = na.omit(melt(fit,id = c("time","species","model")))
fit$variable = factor(fit$variable,levels = c("med_dat","fit","med"))
fit$variable = revalue(fit$variable,c("med_dat"="Data Median","fit"="ODE Fit","med"="Stochastic Simulation\nMedian"))
 
AIC  =unique(all_dat[,c("AIC","model","species")])
AIC = subset(AIC, species=="Salivary\nGlands")

(plot_fit = ggplot()+
  geom_point(data = dat_point,aes(x = time,y =dat,colour = species),size = 2,alpha = 0.1)+
  geom_line(data = subset(fit,variable!="Stochastic Simulation\nMedian"),aes(x = time,y = value,linetype = variable,colour = species),linewidth=1)+
  facet_grid(species~model,scales = "free")+theme_bw()+
    scale_x_continuous(breaks = c(0,8,16,24,32),limits = c(0,32))+
  geom_text(data    = AIC,
            mapping = aes(x = 25, y = 6.25, label = paste0("AIC=",round(AIC,0))),size = 4)+
  labs(x = "Days Post Infection",y = "Units",colour = "Units",fill = "Units",linetype = "Line Type")+theme_bw()+
  scale_fill_discrete(labels=c(bquote(log[10](photons/s)), bquote(log[10](photons/s)), "Cells/100 CD8+ T Cells"))+
  scale_colour_manual(labels=c(bquote(log[10](photons/s)), bquote(log[10](photons/s)), "Cells/1000 CD8+ T Cells"),values=c("#009999","#6600CC","#33CC33"))+
  theme(text = element_text(size=10),
        legend.title=element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),
        axis.title = element_text(size=10),
        strip.text = element_text(size = 10)))





