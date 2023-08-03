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

stand <- function(x){ med = x - mean(x)
         return(med/sd(med))		}

datos[,5:16] <- apply(datos[,5:16], 2, stand)
attach(datos)

# 1
zeros_year <- data.frame(aggregate((total == 0)*1, by=list(anhon), sum)) 
nonzeros_year <- data.frame(aggregate((total > 0)*1, by=list(anhon), sum)) 
zeros_year2 <- data.frame(rep(zeros_year$Group.1, times=2))
zeros_year2$values <- c(nonzeros_year$x, zeros_year$x)
zeros_year2$Type <- rep(c('>0', '=0'), each = 9)
names(zeros_year2) <- c('Year', 'Values', 'Outcome')

ggplot(zeros_year2, aes(x = Year, y = Values, fill = Outcome)) + geom_bar(alpha = 1, stat='identity')  +
  theme_bw(base_size = 16) +  scale_fill_manual( values = c("#2F4F4F", "#56B4E9")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="top") + xlab("Year") + ylab("")

# 2
zeros_week <- data.frame(aggregate((total == 0)*1, by=list(semanan), sum)) 
nonzeros_week <- data.frame(aggregate((total > 0)*1, by=list(semanan), sum)) 
zeros_week2 <- data.frame(rep(zeros_week$Group.1, times=2))
zeros_week2$values <- c(nonzeros_week$x, zeros_week$x)
zeros_week2$Type <- rep(c('>0', '=0'), each = 18)
names(zeros_week2) <- c('Week', 'Values', 'Outcome')

ggplot(zeros_week2, aes(x = Week, y = Values, fill = Outcome)) + geom_bar(alpha = 1, stat='identity')  +
  theme_bw(base_size = 16) + scale_fill_manual( values = c("#2F4F4F", "#56B4E9")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="top") +
  xlab("Week") + ylab("")

# 3
zeros_prov <- data.frame(aggregate((total == 0)*1, by=list(id_provincia), sum)) 
nonzeros_prov <- data.frame(aggregate((total > 0)*1, by=list(id_provincia), sum)) 
zeros_prov2 <- data.frame(rep(zeros_prov$Group.1, times=2))
zeros_prov2$values <- c(nonzeros_prov$x, zeros_prov$x)
zeros_prov2$Type <- rep(c('>0', '=0'), each = 41)
names(zeros_prov2) <- c('Province', 'Values', 'Outcome')

ggplot(zeros_prov2, aes(x = Province, y = Values, fill = Outcome)) + geom_bar(alpha = 1, stat='identity')  +
  theme_bw(base_size = 16) +  scale_fill_manual( values = c("#2F4F4F", "#56B4E9")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="top") + scale_x_discrete(breaks = function(x){x[c(F, F, F, T, F, F, F, F)]}) +
  xlab("Province") + ylab("")

########################################
######## PART 2 
########################################
#####
# 0. To model the zeros separately from the non-zeros with a Bernoulli mixed model
#####

id <- 1:dim(datos)[1]
non_zero <- ifelse(datos$total > 0, 1, 0)
zero <- ifelse(datos$total > 0, 0, 1)
datos$non_zero <- non_zero
datos$zero <- zero

# Only AEMET variables (we eliminate the non significant ones).
# prec, tmin, dir
modelo0 <- glmmTMB(zero ~  sol + tmax + tmed + velmedia + (1|domain),
		 data = subset(datos, anhon < 2015), family = binomial(link = "logit"))
summary(modelo0)
round(confint(modelo0),3)


#####
# Semicontinuous data (many exact zeros and continuous positive outcomes)
# 1. To model the total burned area:
#####

modeloA <- glmmTMB(total ~ medios_autobomba + medios_bulldozer + medios_tractores + 
                   aereos_aviones_anfibios + aereos_aviones_carga + 
                   aereos_helicopteros_transporte + tmax + tmed + 
                   velmedia + (1|domain), 
                    
                   ziformula = ~ sol + tmax + tmed + velmedia + (1|domain),
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
ci_cor(log(datos$promedio[j]), dir[j])

datos$est_incendios <- (num_incendios- mean (num_incendios)) / sd(num_incendios- mean (num_incendios))

modeloB <- glmmTMB(promedio ~ est_incendios + dist_1_nearest_buildup + tmax + 
                     tmed + velmedia + (1|domain), 
                   
                   ziformula = ~ sol + tmax + tmed + velmedia + (1|domain),
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

# Standardized residuals
res.btB = datosfit$total/datosfit$num_incendios - fitted(modeloB)
cd <- ( datosfit$num_incendios==0 )
res.btB[cd] <- (- fitted(modeloB) )[cd]
res.stB = (res.btB - mean(res.btB))/sd(res.btB)

datosfit$res.stB <- res.stB
datosfit$res.btB <- res.btB

# Usual plot
datosfit$am <- 1
datosfit$am[datosfit$promedio<500] <- 0
datosfit$am <- as.factor(datosfit$am)
datosfit$Average <- datosfit$promedio
datosfit$zero2 <- 'Non zero'
datosfit$zero2[datosfit$zero==1] <- 'Zero'
ggplot(data = datosfit, aes(x = 1:(dim(datosfit)[1]), y = res.stB)) + 
  geom_point(aes(colour = am, size=Average)) +
  facet_grid(~zero2) + scale_colour_manual(values=c("black", "#E69F00"), guide = "none") + 
  theme_bw(base_size = 13.5) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("Domain") + ylab("Residuals")


# Boxplots (deleting GFFs)
subset(datosfit, promedio<500) %>% ggplot(aes( x= anho, y = res.stB)) + geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 17.5) + xlab("Year") + ylab(" ") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(breaks = function(x){x[c(F, T)]})
  
subset(datosfit, promedio<500) %>% ggplot(aes( x= semana, y = res.stB)) + geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 17.5) + xlab("Week") + ylab(" ") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(breaks = function(x){x[c(F, F, F, T, F)]})

subset(datosfit, promedio<500) %>% ggplot(aes( x= id_provincia, y = res.stB)) + geom_boxplot(alpha = 0.8) +
  theme_bw(base_size = 17.5) + xlab("Province") + ylab(" ") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(breaks = function(x){x[c(F, F, F, T, F, F, F, F)]})

# GFFs
summary(subset(datosfit, promedio>=500)$res.stB)
summary(subset(datosfit, promedio>=500)$res.btB)

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
pred2015$count <- num_incendios[datos$anhon == 2015]
pred2015$non_zero <- non_zero[datos$anhon == 2015]
pred2015$semana <- semana[datos$anhon == 2015]
pred2015$id_provincia <- id_provincia[datos$anhon == 2015]

pred2015 %>% rename('aZIGT    ' = predA, 'aZIGA    ' = predB, 'Observed    ' = logquemado) %>%
  pivot_longer(c('aZIGT    ', 'aZIGA    ', 'Observed    '), names_to = "Source") %>%
  ggplot(aes( x= semana, y = value, fill = Source)) +
  geom_boxplot(alpha = 0.8) + coord_cartesian(ylim=c(-3, 15)) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme_bw(base_size = 13.5) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = c(0.5, 0.9), legend.direction="horizontal",
        legend.text = element_text(size=16), legend.title = element_text(size=0)) +
  xlab("Week") + ylab("log Ha")

pred2015 %>% rename('aZIGT  ' = predA, 'aZIGA  ' = predB, 'Observed  ' = logquemado) %>%
  pivot_longer(c('aZIGT  ', 'aZIGA  ', 'Observed  '), names_to = "Source") %>%
  ggplot(aes( x= id_provincia, y = value, fill = Source)) +
  geom_boxplot(alpha = 0.8) + coord_cartesian(ylim=c(-3, 15))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme_bw(base_size = 13.5) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
  	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
  	legend.position = 'c(0.5, 1.85)', legend.text = element_text(size=15),
  	legend.title = element_text(size=0)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  xlab("Province") + ylab("log Ha")
 

########################################
######## PART 5 : Clustering. K-means
########################################

set.seed(65)
kmeans_out <- kmeans(datos[, c(5:16, 22:30, 36:37)], 3, iter.max = 1000, nstart = 10)

index_voraces1 <- (1:6642)[kmeans_out$cluster==1]
voraces1 <- datos[index_voraces1, ]

index_voraces2 <- (1:6642)[kmeans_out$cluster==2]
voraces2 <- datos[index_voraces2, ]

index_voraces3 <- (1:6642)[kmeans_out$cluster==3]
voraces3 <- datos[index_voraces3, ]

   