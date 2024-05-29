###############################################################################
###############################################################################
###       Modelo Gamma para predecir colapsos y ventanas de oportunidad
###
###       Created on Tuesday April 21 2023                                                  
###                                                     
###       @author: mariabugalloporto

           
rm(list=ls())           
library(maptools)
library(RColorBrewer)

########################################
######## PART 1 : Main functions
########################################

GroupClassification <- function(data,datacompare,intervals)
{
   n = length(data)
   group = matrix(0,nrow=n) 
   ninterv = length(intervals)
   for (i in 1:n)
   {
      for (j in 1:ninterv)
         if (datacompare[i]<intervals[j])
         {
            group[i]= intervals[j]
            break
         }
   }
   result = split(data,group)
   return (result)
}
 
PrintSpainMap <- function(pathmap,datos,colors,titlemap,textlegend)
{

   m <- matrix(c(1,1,1,2),2,2)
   layout(m, widths=c(1.5, 1), heights=c(1.5, 1), respect=F)

   xName <- readShapePoly(pathmap, IDvar="NAME", proj4string=CRS("+proj=longlat +ellps=clrk66"))     

   xName$datos <- NA
   for (i in 1:length(colors))
      xName$datos[datos[[i]]] <- colors[i]
 
   xSC <- xName[xName$ESP_PROV_I < 35 | xName$ESP_PROV_I >38 | xName$ESP_PROV_I==36 | xName$ESP_PROV_I ==37,]
   plot(xSC,  xlab="",  col=xSC$datos, axes=F)
  title(titlemap, line=-0.1, cex.main=1.7)
  legend( "topright", textlegend, col = 1, pt.bg = colors, pch=21, bty = "n", cex=1.2, pt.cex = 2.5)  #  cex=1.3

   xC <- xName[xName$ESP_PROV_I == 35 | xName$ESP_PROV_I == 38,]
   plot(xC,  xlab="",  col=xC$datos)
   box()
}

########################################
######## PART 2 : Paper maps
########################################
pathmap    <- "SPAIN_MAP/esp_prov.shp"

##############################################################

library(mapSpain)
library(tidyverse)

provincias <- esp_get_prov(year='2021')

cpro.index <- as.numeric(levels(as.factor(read.csv("original_prov_week.csv", header=T, sep='|')$id_provincia)))

colors <- rep('No data*', 52); 
colors[cpro.index] <- 'Peninsular Center'
colors[c(1, 15, 20, 24, 27, 32, 33, 36, 39, 48, 49)[c(1, 15, 20, 24, 27, 32, 33, 36, 39, 48, 49) %in% cpro.index]] <- 'Northwest Spain' #Norte

colors[c(3, 4, 8, 11, 12, 14, 17, 18, 21, 23, 25, 29, 30, 41, 43, 46)[c(3, 4, 8, 11, 12, 14, 17, 18, 21, 23, 25, 29, 30, 41, 43, 46) %in% cpro.index]] <- 'Mediterranean Coast' # Sur

colors.df <- data.frame('cpro' = 1:52, 'colors' = factor(colors, levels=c('No data*', 'Peninsular Center', 'Mediterranean Coast', 'Northwest Spain')))
colors.df$cpro <- ifelse(colors.df$cpro < 10, paste0('0', as.factor(colors.df$cpro)), as.factor(colors.df$cpro)) 

provincias <- merge(provincias, colors.df, by='cpro')

ggplot(provincias) +
  geom_sf(data = esp_get_can_box(), color = "grey70") +
  geom_sf(data = esp_get_can_provinces(), color = "grey70") +
  geom_sf(aes(fill = colors), color = "grey70") +
  labs(title = "") +
  scale_fill_discrete(' ', type = hcl.colors(4, "Blues")[4:1])+
  theme_bw() + theme(text = element_text(size=15), legend.position=c(.833,.166))
 


##############################################################
# Prediction Intervals
results <- read.csv("predBboot.csv",  header=T)

dom <- 1:50
id_prov <- cpro.index
                                                
results32 <- results[ results$semana==29, ]
	
	estML <- rep(NA, 50)
	estML[ id_prov ] <- 1
	estML[ id_prov ][ results32$Dicot4 ==1 ] <- 3
	estML[ id_prov ][ results32$Dicot3 ==1 ] <- 5
	estML[ id_prov ][ results32$Dicot2 ==1 ] <- 7
	estML[ id_prov ][ results32$Dicot1 ==1 ] <- 9
	estML[is.na(estML)] <- 11
	
	# Intervals  
	intervals_prop <- c(2, 4, 6, 8, 10, Inf)   

	# Colors
	colorsprop <- c(  brewer.pal(7,"GnBu")[c(2,3,4,5,6)] , hcl.colors(4, "Blues")[4])
	
	# Legend
	legend_prop <- expression("< 70 %", " 70 - 80 %", "80 - 85 %", "85 - 95 %", "> 95 %") 

	result = GroupClassification(dom,estML,intervals_prop)
	PrintSpainMap(pathmap,result,colorsprop,"13/07/2015 - 19/07/2015",legend_prop)


##############################################################
# RRMSE	Model B and Year = 2015
RRMSE <- data.frame('RRMSEB' = tail(read.csv("rmse500.csv")$rrmseB, 738)) 
RRMSE$anho <- results$anho
RRMSE$semana <- results$semana
RRMSE$idprovincia <- results$idprovincia

estML <- rep(100, 50)
estML[ id_prov ] <- RRMSE$RRMSEB[ RRMSE$semana==28 ]	

# Intervals  
intervals_prop <- c(3, 5, 10, 15, 90, Inf)

# Colors
colorsprop <- c(brewer.pal(8,'YlOrRd')[c(2,3,4,5,6)], hcl.colors(4, "Blues")[4])
	
# Legend
legend_prop <- expression("under 3 %", " 3 - 5 %", "5 - 10 %", "10 - 15 %", " 15- 25 %") 

result = GroupClassification(dom,estML,intervals_prop)
PrintSpainMap(pathmap,result,colorsprop,"06/07/2015 - 12/07/2015",legend_prop)


