###############################################################################
###############################################################################
###       Modelo Gamma para predecir colapsos y ventanas de oportunidad
###
###       Created on Wednesday March 16 2023
###       Bootstrap Confidence Intervals and MSE
###                                                     
###       @author: mariabugalloporto

rm(list=ls())
library(glmmTMB)
library(ggplot2)
options(warn=-1)

set.seed(1605)

########################################
######## PART 1 
########################################
# Load data
datos <- read.csv("features_prov_week.csv")
datos$anhon <- datos$anho
datos$semanan <- datos$semana
datos$anho <- factor(datos$anho)
datos$semana <- factor(datos$semana)
datos$id_provincia <- factor(datos$id_provincia)
datos$provincia <- factor(datos$provincia)
datos$domain <- factor(datos$domain)
datos$dist_1_nearest_buildup <- datos$dist_1_nearest_buildup/1000
datos$dist_10_nearest_buildup <- datos$dist_10_nearest_buildup/1000

stand <- function(x){ med = x - mean(x)
return(med/sd(med))
}

datos[,5:16] <- apply(datos[,5:16], 2, stand)
attach(datos)

I <- 9; J <- 18; K <- 41

#####
# Semicontinuous data (many exact zeros and continuous positive outcomes)
# 1. To model the total burned area:
#####

modeloA <- glmmTMB(total ~ medios_autobomba + medios_bulldozer + medios_tractores + 
                     aereos_aviones_anfibios + aereos_aviones_carga + 
                     aereos_helicopteros_transporte + tmax + tmed + 
                     velmedia + (1|domain), 
                   
                   ziformula = ~ sol + tmax + tmed + velmedia + (1|domain),
                   data = datos,
                   family = ziGamma(link = "log"))

mu_PIAA <- predict(modeloA, datos, type="response")

#####
# 2. To model the average burned area:
#####

datos$promedio <- total/num_incendios
datos$promedio[is.na(datos$promedio)] <- 0
datos$est_incendios <- (num_incendios- mean (num_incendios)) / sd(num_incendios- mean (num_incendios))

modeloB <- glmmTMB(promedio ~ est_incendios + dist_1_nearest_buildup + tmax + 
                     tmed + velmedia + (1|domain), 
                   
                   ziformula = ~ sol + tmax + tmed + velmedia + (1|domain),
                   data = datos,
                   family = ziGamma(link = "log"))

mu_PIBB <- predict(modeloB, datos, type="response")

########################################
######## PART 2: Bootstrap resampling
########################################	

R = 500 # number of simulations
nobs <- dim(datos)[1]

nparamA <- 10+5+1+2
nparamB <- 6+5+1+2
est_boot_parametrosA <- matrix(ncol=nparamA, nrow=R)
est_boot_parametrosB <- matrix(ncol=nparamB, nrow=R)

difA <- matrix(ncol=6642, nrow=R)
difB <- matrix(ncol=6642, nrow=R)

t <- proc.time() # time

for(r in 1:R){
  print(r)
  
  y_bootA <- simulate(modeloA)[[1]]
  datos$y_bootA <- y_bootA
  try(modeloAboot <- glmmTMB(y_bootA ~ medios_autobomba + medios_bulldozer + medios_tractores + 
                               aereos_aviones_anfibios + aereos_aviones_carga + 
                               aereos_helicopteros_transporte + tmax + tmed + 
                               velmedia + (1|domain), 
                             
                             ziformula = ~ sol + tmax + tmed + velmedia + (1|domain),
                             data = datos,
                             family = ziGamma(link = "log")))
      
  est_boot_parametrosA[r,]=modeloAboot$sdr$par.fixed    
  mu_PIA <- predict(modeloAboot, datos, type="response")
  
  difA[r, ] <- (y_bootA - mu_PIA)^2
  
  y_bootB <- simulate(modeloB)[[1]]
  datos$y_bootB <- y_bootB
  try(modeloBboot <- glmmTMB(y_bootB ~ est_incendios + dist_1_nearest_buildup + tmax + 
                               tmed + velmedia + (1|domain), 
                             
                             ziformula = ~ sol + tmax + tmed + velmedia + (1|domain),
                             data = datos,
                             family = ziGamma(link = "log")))
  
  est_boot_parametrosB[r,]=modeloBboot$sdr$par.fixed
  mu_PIB <- predict(modeloBboot, datos, type="response")	
  
  difB[r, ] <- (y_bootB - mu_PIB)^2
}		

proc.time()-t # time	

# Confidence intervals
apply(est_boot_parametrosA, 2, quantile, probs=c(0.025, 0.975))
apply(est_boot_parametrosB, 2, quantile, probs=c(0.025, 0.975))

rmseA <- sqrt(colMeans(difA, na.rm=F))
rmseB <- sqrt(colMeans(difB, na.rm=F))

rrmseA <- rmseA / mu_PIAA
rrmseB <- rmseB / mu_PIBB

save_results <- data.frame(rmseA, rmseB, rrmseA, rrmseB)
write.csv(save_results,"rmse500.csv", row.names = FALSE)

# Usual plot
plot_results <- save_results[,3:4]

zero <- ifelse(datos$total > 0, 0, 1)
plot_results$Outcome <- '>0'
plot_results$Outcome[zero==1] <- '=0'

plot_results <- plot_results[datos$anhon==2015, ]

ggplot(data = plot_results, aes(x = 1:dim(plot_results)[1], y = rrmseA)) + 
  geom_point(aes(colour = Outcome), size=3) +
  scale_colour_manual(values=c("#2F4F4F", "#56B4E9")) + 
  theme_bw(base_size = 25) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="top") +
  xlab("Domain") + ylab("RRMSE (%)")

ggplot(data = plot_results, aes(x = 1:dim(plot_results)[1], y = rrmseB)) + 
  geom_point(aes(colour = Outcome), size=3) +
  scale_colour_manual(values=c("#2F4F4F", "#56B4E9")) + 
  theme_bw(base_size = 25) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="top") +
  xlab("Domain") + ylab("RRMSE (%)")

ggplot(data = plot_results, aes(x = rrmseB, y = rrmseA)) + 
  geom_point(aes(colour = Outcome), size=2.7) +
  scale_colour_manual(values=c("#2F4F4F", "#56B4E9")) + 
  theme_bw(base_size = 25) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="top") + xlim(0, 272) +
  xlab("aZIGA (%)") + ylab("aZIGT (%)")


round(quantile(plot_results[,1], probs=seq(0,1,by=0.25)), 3) ## RRMSEA 2015
round(quantile(plot_results[,2], probs=seq(0,1,by=0.25)), 3) ## RRMSEB 2015

