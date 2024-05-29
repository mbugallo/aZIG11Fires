###############################################################################
### 
###       Gamma model to predict collapse points and windows of opportunity
###
###       Created on Tuesday February 14 2023                                                  
###                                                     
###       @author: mariabugalloporto

rm(list=ls())
library(confintr)
library(ggplot2)
library(tidyverse)
library(glmmTMB)
library(car)

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


datos[,5:16] <- apply(datos[,5:16], 2, function(x){ med = x - mean(x);  return(med/sd(med)) })
attach(datos)


index_prov <- rep('Peninsular Center', dim(datos)[1])
index_prov[datos$id_provincia %in% c(1, 15, 20, 24, 27, 32, 33, 36, 39, 48, 49)] <- 'Northwest Spain' #Norte
index_prov[datos$id_provincia %in% c(3, 4, 8, 11, 12, 14, 17, 18, 21, 23, 25, 29, 30, 41, 43, 46)]  <- 'Mediterranean Coast' # Sur



########################################
######## PART 2 
########################################

id <- 1:dim(datos)[1]

#####
# Semicontinuous data (many exact zeros and continuous positive outcomes)
# 1. To model the total burned area:
#####

modeloA <- glmmTMB(total ~ -1+ medios_autobomba + medios_bulldozer + medios_tractores + 
                     aereos_aviones_anfibios + aereos_aviones_carga + 
                     aereos_helicopteros_transporte + prec + tmax+  tmed  +velmedia+ hrmedia + (1|domain), 
                
                   ziformula = ~ -1 + prec + tmax + tmed + velmedia + hrmedia + (1|domain),
                   data = subset(datos, anhon < 2015),
                   family = ziGamma(link = "log"))
summary(modeloA)
round(confint(modeloA),3)

#####
# 2. To model the average burned area:
#####
  
datos$promedio <- total/num_incendios
datos$promedio[is.na(datos$promedio)] <- 0

# Pearson's correlation coefficients between the logarithm of the count and the auxiliary variables
j <- (total > 0)
ci_cor(log(datos$promedio[j]), hrmedia[j])

datos$est_incendios <- (num_incendios- mean (num_incendios)) / sd(num_incendios- mean (num_incendios))

modeloB <- glmmTMB(promedio ~ -1 + est_incendios + dist_1_nearest_buildup + prec + tmax + 
                     tmed + velmedia + hrmedia + (1|domain), 
                   
                   ziformula = ~ -1 +prec + tmax + tmed + velmedia + hrmedia + (1|domain),
                   data = subset(datos, anhon < 2015),
                   family = ziGamma(link = "log"))
summary(modeloB) 
round(confint(modeloB),3)

########################################
######## PART 3 : Residuals 
########################################
#####
# 2. To model the average burned area:
#####

datosfit <- subset(datos, anhon < 2015)
datosfit$zero <- (subset(datos, anhon < 2015)$promedio==0)*1

# Standardized residuals
res.btB = datosfit$total/datosfit$num_incendios - fitted(modeloB)
cd <- ( datosfit$num_incendios==0 )
res.btB[cd] <- (- fitted(modeloB) )[cd]
res.stB = (res.btB - mean(res.btB))/sd(res.btB)

datosfit$res.stB <- res.stB
datosfit$res.btB <- res.btB

par(mfrow=c(1,3))

# Usual plot
datosfit$am <- 1
datosfit$am[datosfit$promedio<500] <- 0
datosfit$am <- as.factor(datosfit$am)
datosfit$Average <- datosfit$promedio
datosfit$zero2 <- 'Non zero'
datosfit$zero2[datosfit$zero==1] <- 'Zero'
ggplot(data = datosfit, aes(x = 1:(dim(datosfit)[1]), y = res.stB)) + 
  geom_point(aes(colour = am, size=Average)) +
  facet_grid(~zero2) + scale_colour_manual(values=c("black", "darkblue"), guide = "none") + 
  theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
  axis.line = element_line(colour = "black"), text = element_text(size=16)) +
  xlab("Domain") + ylab("Standarized residuals")


# Random Effects
u1 <- ranef(modeloB)$zi[[1]]$`(Intercept)` / sqrt(summary(modeloB)$varcor$zi$domain[1])
u2 <- ranef(modeloB)$cond[[1]]$`(Intercept)` / sqrt(summary(modeloB)$varcor$cond$domain[1])


########################################
######## PART 4 : Fitting
########################################

predA <- log(predict(modeloA, newdata = subset(datos, anhon == 2015), type='response'))
predB <- log(predict(modeloB, newdata = subset(datos, anhon == 2015),  type='response')*num_incendios[datos$anhon == 2015])

pred2015 <- data.frame('id' = id[datos$anhon == 2015], predA, predB)
pred2015$quemado <- total[datos$anhon == 2015]
predB[pred2015$quemado == 0] <- 0
pred2015$predB <- predB
pred2015$logquemado <- log(pred2015$quemado)
pred2015$logquemado[pred2015$quemado == 0] <- 0
pred2015$semana <- semana[datos$anhon == 2015]
pred2015$index_prov <- factor(index_prov[datos$anho == '2015'], levels=c('Northwest Spain', 'Peninsular Center', 'Mediterranean Coast'))


datos_index <- aggregate(data.frame(predA, predB, pred2015$logquemado), by=list(index_prov[datos$anho == '2015'], pred2015$semana ), mean)  

par(mar=c(4, 4.1, 4, 4), mfrow=c(1,3))
plot(datos_index$predA[datos_index[,1]=='Northwest Spain'], col='#999999', type='l', lwd=2, ylim=c(-1.5,6),
		 ylab='log Ha', xlab='Week', cex.axis=1.6, cex.lab=1.6, cex.main=1.6, main='Northwest Spain')
lines(datos_index$predB[datos_index[,1]=='Northwest Spain'], col='#E69F00', type='l', lwd=2)
lines(datos_index$pred2015.logquemado[datos_index[,1]=='Northwest Spain'], col='#56B4E9', type='l', lwd=2)

plot(datos_index$predA[datos_index[,1]=='Peninsular Center'], col='#999999', type='l', lwd=2, ylim=c(-1.5,6),
		 ylab='log Ha', xlab='Week', cex.axis=1.6, cex.lab=1.6, cex.main=1.6, main='Peninsular Center')
lines(datos_index$predB[datos_index[,1]=='Peninsular Center'], col='#E69F00', type='l', lwd=2)
lines(datos_index$pred2015.logquemado[datos_index[,1]=='Peninsular Center'], col='#56B4E9', type='l', lwd=2)


plot(datos_index$predA[datos_index[,1]=='Mediterranean Coast'], col='#999999', type='l', lwd=2, ylim=c(-1.5,6),
		 ylab='log Ha', xlab='Week', cex.axis=1.6, cex.lab=1.6, cex.main=1.6, main='Mediterranean Coast')
lines(datos_index$predB[datos_index[,1]=='Mediterranean Coast'], col='#E69F00', type='l', lwd=2)
lines(datos_index$pred2015.logquemado[datos_index[,1]=='Mediterranean Coast'], col='#56B4E9', type='l', lwd=2)

legend("topright", legend=c("aZIGA","aZIGT ", "Observed"), lty=c(1,1), col=c('#999999', '#E69F00', '#56B4E9'),
       cex=1.6, lwd=2, bty = "n")	

 

########################################
######## PART 5 : Clustering. K-means
########################################

set.seed(65)
kmeans_out <- kmeans(datos[, c(5:16, 22:31, 35:36)], 3, iter.max = 1000, nstart = 10)

index_voraces1 <- (1:6642)[kmeans_out$cluster==1]
voraces1 <- datos[index_voraces1, ]

index_voraces2 <- (1:6642)[kmeans_out$cluster==2]
voraces2 <- datos[index_voraces2, ]

index_voraces3 <- (1:6642)[kmeans_out$cluster==3]
voraces3 <- datos[index_voraces3, ]
   