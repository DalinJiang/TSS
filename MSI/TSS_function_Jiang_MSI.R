#******************************************************************************************************************************************************************
#
# This R code is the function for estimating TSS from Sentinel 2 MSI bands using the method proposed by Jiang et al. (2023)
#
# *Note: this file only include the R function to estimate TSS, to process your own data (either satellite image or in situ Rrs data),
# 		 you will need another R code to call this function.
#
#
#-------------
# inputs are:
# Remote sensing reflectance (Rrs, sr-1) at MSI bands
# i.e., 443,490,560,665,705,740,783,865 nm
#
#-------------
# outputs are:
# 	(1) estimated absorption coefficient (a, m-1) at reference band
#	(2) particulate backscatering coefficient (bbp, m-1) at reference band
# 	(3) the reference band (band, nm) used for TSS estimation
#	(4) estimated TSS (g/m3)
# 
#-------------
# for more informaiton, please refer to the publication:
# Jiang et al. (2023). Estimating the concentration of total suspended solids in inland and coastal waters from Sentinel-2 MSI: a semi-analytical approach. 
#						ISPRS International Journal of Photogrammetry and Remote sensing.
#
#
#
# Version = V1.0
# Dalin JIANG, University of Stirling, UK
# 21 August, 2023
#
#
#******************************************************************************************************************************************************************
#

Estimate_TSS_Jiang_MSI <- function(site_Rrs){
	aw <- c(0.00515124,0.01919594,0.06299986,0.41395333,0.70385758,2.71167020,2.62000141,4.61714226)
	bbw <- c(0.00215037,0.00138116,0.00078491,0.00037474,0.00029185,0.00023499,0.00018516,0.00012066)
	MSI_wave <- c(443,490,560,665,705,740,783,865)
	names(aw) <- paste("aw",MSI_wave,sep="")
	names(bbw) <- paste("bbw",MSI_wave,sep="")
	names(site_Rrs) <- paste("Rrs",MSI_wave,sep="")
	#
	#----------560nm---------------------
	QAA_560 <- function(site_Rrs){
		rrs <- site_Rrs/(0.52+1.7*site_Rrs)
		u <- (-0.089+sqrt((0.089^2)+4*0.125*rrs))/(2*0.125)
		x <- log((rrs["Rrs443"]+rrs["Rrs490"])/(rrs["Rrs560"]+5*rrs["Rrs665"]*rrs["Rrs665"]/rrs["Rrs490"]),10)
		a560 <- aw["aw560"]+10^(-1.146-1.366*x-0.469*(x^2))            
		bbp560 <- ((u["Rrs560"]*a560)/(1-u["Rrs560"]))-bbw["bbw560"]
		one_tss <- 94.48785*bbp560
		bbp_wave <- 560
		return(c(a560,bbp560,bbp_wave,one_tss))
	}
	#----------665nm---------------------
	QAA_665 <- function(site_Rrs){
		rrs <- site_Rrs/(0.52+1.7*site_Rrs)
		u <- (-0.089+sqrt((0.089^2)+4*0.125*rrs))/(2*0.125)
		a665 <- aw["aw665"]+0.39*((site_Rrs["Rrs665"]/(site_Rrs["Rrs443"]+site_Rrs["Rrs490"]))^1.14)    
		bbp665 <- ((u["Rrs665"]*a665)/(1-u["Rrs665"]))-bbw["bbw665"]
		one_tss <- 113.87498*bbp665
		bbp_wave <- 665
		return(c(a665,bbp665,bbp_wave,one_tss))
	}  
	#----------740nm---------------------
	QAA_740 <- function(site_Rrs){
		rrs <- site_Rrs/(0.52+1.7*site_Rrs)
		u <- (-0.089+sqrt((0.089^2)+4*0.125*rrs))/(2*0.125)
		bbp740 <- ((u["Rrs740"]*aw["aw740"])/(1-u["Rrs740"]))-bbw["bbw740"]
		one_tss <- 134.91845*bbp740
		bbp_wave <- 740      
		return(c(aw["aw740"],bbp740,bbp_wave,one_tss))
	}
	#-----------865nm--------------------
	QAA_865 <- function(site_Rrs){
		rrs <- site_Rrs/(0.52+1.7*site_Rrs)
		u <- (-0.089+sqrt((0.089^2)+4*0.125*rrs))/(2*0.125)
		bbp865 <- ((u["Rrs865"]*aw["aw865"])/(1-u["Rrs865"]))-bbw["bbw865"]
		one_tss <- 166.07382*bbp865
		bbp_wave <- 865
		return(c(aw["aw865"],bbp865,bbp_wave,one_tss))
	}
	#-----------estimate Rrs620 from Rrs665----------
	estimate_Rrs620 <- function(in665){
		a <- 1.693846e+02                   
		b <- -1.557556e+01                  
		c <- 1.316727e+00                   
		d <- 1.484814e-04                  
	
		est620 <- a*in665^3+b*in665^2+c*in665+d
	
		return(est620)
	}
  
	Rrs620 <- estimate_Rrs620(site_Rrs["Rrs665"])   
	       
	#------------band selection----------------------------
	# Rrs490, Rrs560, Rrs665 and Rrs740 cannot be NA, otherwise there will be errors
	if (all(site_Rrs == 0) | all(is.na(site_Rrs) == TRUE) | any(is.na(c(site_Rrs["Rrs490"],site_Rrs["Rrs560"],site_Rrs["Rrs665"],site_Rrs["Rrs740"])))){        
		tmp_tss <- rep(NA,4)
	}else {
		if (site_Rrs["Rrs490"] > site_Rrs["Rrs560"]){
			tmp_tss <- QAA_560(site_Rrs)
		}else if (site_Rrs["Rrs490"] > Rrs620){
			tmp_tss <- QAA_665(site_Rrs)    
		}else if (site_Rrs["Rrs740"] > site_Rrs["Rrs490"] & site_Rrs["Rrs740"] > 0.010){
			tmp_tss <- QAA_865(site_Rrs)      
		}else{
			tmp_tss <- QAA_740(site_Rrs)     
		}  	
	}
	tmp_tss <- as.numeric(tmp_tss)   
	names(tmp_tss) <- c("a","bbp","band","TSS")
	#
	return(tmp_tss)
}