###############################################################################
###############################################################################
### 
###       Gamma model to predict collapse points and windows of opportunity
###
###       Created on Tuesday February 14 2023                                              
###                                                     
###       @author: mariabugalloporto

rm(list=ls())

########################################
######## PART 1 
########################################

# Delete Cantabria, Madrid, Euskadi, Navarra, Baleares, Canarias, Ceuta y Melilla

feat <- read.csv("original_prov_week.csv", header=T, sep='|') 
id_year <- as.numeric(levels(as.factor(feat$anho)))
id_prov <- as.numeric(levels(as.factor(feat$id_provincia)))

AEMETfeat <- read.csv("original_AEMET_prov_week.csv", header=T, sep=',') 
AEMETfeat <- AEMETfeat[ AEMETfeat$anho %in% id_year, ]
AEMETfeat <- AEMETfeat[ AEMETfeat$id_provincia %in% id_prov, ]
AEMETfeat$domain <- as.factor(AEMETfeat$semana) : as.factor(AEMETfeat$id_provincia)

median_notNaN <- function(x){ median(x[x!=-9999]) }
median_domain <- data.frame(aggregate(AEMETfeat[, 5:11], by=list(AEMETfeat$domain), median_notNaN))
colnames(median_domain)[1] <- 'domain'

for (k in levels(AEMETfeat$domain)){
	AEMETfeat$dir[ AEMETfeat$dir==-9999 & AEMETfeat$domain==k] = median_domain$dir[ median_domain$domain == k ]
	AEMETfeat$prec[ AEMETfeat$prec==-9999 & AEMETfeat$domain==k] = median_domain$prec[ median_domain$domain == k ]
	AEMETfeat$sol[ AEMETfeat$sol==-9999 & AEMETfeat$domain==k] = median_domain$sol[ median_domain$domain == k ]
	AEMETfeat$tmax[ AEMETfeat$tmax==-9999 & AEMETfeat$domain==k] = median_domain$tmax[ median_domain$domain == k ]
	AEMETfeat$tmed[ AEMETfeat$tmed==-9999 & AEMETfeat$domain==k] = median_domain$tmed[ median_domain$domain == k ]
	AEMETfeat$tmin[ AEMETfeat$tmin==-9999 & AEMETfeat$domain==k] = median_domain$tmin[ median_domain$domain == k ]
	AEMETfeat$velmedia[ AEMETfeat$velmedia==-9999 & AEMETfeat$domain==k] = median_domain$velmedia[ median_domain$domain == k ]
	AEMETfeat$hrmedia [ AEMETfeat$vhrmedia ==-9999 & AEMETfeat$domain==k] = 
median_domain$hrmedia [ median_domain$domain == k ]
	
}

preproc_feat <- cbind(feat, AEMETfeat[, 5:13] )

# Weeks 27 (1st July) to 44 (4th October).
preproc_feat_Jul_Oct <- preproc_feat[ preproc_feat$semana %in% 27:44, ]

preproc_feat_Jul_Oct$dist_1_nearest_buildup[is.na(preproc_feat_Jul_Oct$dist_1_nearest_buildup)] <- 0 
preproc_feat_Jul_Oct$dist_10_nearest_buildup[is.na(preproc_feat_Jul_Oct$dist_10_nearest_buildup)] <- 0 

write.csv(preproc_feat_Jul_Oct, "features_prov_week.csv", row.names = FALSE)




