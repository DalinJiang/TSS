#******************************************************************************************************************************************************************
#
# This R code is the function for estimating TSS from MERIS or Sentinel-3 OLCI bands using the method proposed by Jiang et al. (2021)
#
# *Note: this file only include the R function to estimate TSS, to process your own data (either satellite image or in situ Rrs data),
# 		 you will need another R code to call this function.
#
#
#-------------
# inputs are:
# Remote sensing reflectance (Rrs, sr-1) at MERIS or OLCI bands
# i.e., 443,490,560,620,665,754,865 nm
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
# Jiang et al.(2021), Remotely estimating total suspended solids concentration in clear to extremely turbid waters using a novel semi-analytical method. 
#					  Remote Sensing of Environment.
#
#
#
# Version = V1.0
# Dalin JIANG, University of Tsukuba, Japan
# 09 November, 2021
#
#
#******************************************************************************************************************************************************************


Estimate_TSS_Jiang_MERIS_OLCI <- function(site_Rrs){
	aw <- c(0.005046443,0.013589323,0.062122106,0.276193682,0.42748488,2.868335728,4.639441062)
	bbw <- c(0.00214135,0.001381358,0.000778527,0.000502851,0.000372427,0.000217139,0.000120218)
	wave <- c(443,490,560,620,665,754,865)
	names(aw) <- paste("aw",wave,sep="")
	names(bbw) <- paste("bbw",wave,sep="")
	names(site_Rrs) <- paste("Rrs",wave,sep="")
	#
	#----------560nm---------------------
	QAA_560 <- function(site_Rrs){
		rrs <- site_Rrs/(0.52+1.7*site_Rrs)
		u <- (-0.089+sqrt((0.089^2)+4*0.125*rrs))/(2*0.125)
		x <- log((rrs["Rrs443"]+rrs["Rrs490"])/(rrs["Rrs560"]+5*rrs["Rrs665"]*rrs["Rrs665"]/rrs["Rrs490"]),10)
		a560 <- aw["aw560"]+10^(-1.146-1.366*x-0.469*(x^2))            
		bbp560 <- ((u["Rrs560"]*a560)/(1-u["Rrs560"]))-bbw["bbw560"]
		bbp_wave <- 560
		one_tss <- 94.6074*bbp560
		return(c(a560,bbp560,bbp_wave,one_tss))
	}
	#----------665nm---------------------
	QAA_665 <- function(site_Rrs){
		rrs <- site_Rrs/(0.52+1.7*site_Rrs)
		u <- (-0.089+sqrt((0.089^2)+4*0.125*rrs))/(2*0.125)
		a665 <- aw["aw665"]+0.39*((site_Rrs["Rrs665"]/(site_Rrs["Rrs443"]+site_Rrs["Rrs490"]))^1.14)   
		bbp665 <- ((u["Rrs665"]*a665)/(1-u["Rrs665"]))-bbw["bbw665"]
		bbp_wave <- 665
		one_tss <- 114.0121*bbp665
		return(c(a665,bbp665,bbp_wave,one_tss))
	}  
	#----------754nm---------------------
	QAA_754 <- function(site_Rrs){
		rrs <- site_Rrs/(0.52+1.7*site_Rrs)
		u <- (-0.089+sqrt((0.089^2)+4*0.125*rrs))/(2*0.125)
		bbp754 <- ((u["Rrs754"]*aw["aw754"])/(1-u["Rrs754"]))-bbw["bbw754"]
		bbp_wave <- 754
		one_tss <- 137.6652*bbp754
		return(c(aw["aw754"],bbp754,bbp_wave,one_tss))
	}
	#-----------865nm--------------------
	QAA_865 <- function(site_Rrs){
		rrs <- site_Rrs/(0.52+1.7*site_Rrs)
		u <- (-0.089+sqrt((0.089^2)+4*0.125*rrs))/(2*0.125)
		bbp865 <- ((u["Rrs865"]*aw["aw865"])/(1-u["Rrs865"]))-bbw["bbw865"]
		bbp_wave <- 865
		one_tss <- 166.1682*bbp865
		return(c(aw["aw865"],bbp865,bbp_wave,one_tss))
	}
	#
	#------------band selection----------------------------
	# in case any band used for water type classification is NA  
	if (all(site_Rrs == 0) | all(is.na(site_Rrs) == TRUE) | any(is.na(c(site_Rrs["Rrs490"],site_Rrs["Rrs560"],site_Rrs["Rrs620"],site_Rrs["Rrs665"],site_Rrs["Rrs754"])))){        
		# for those pixels with NA
		tmp_tss <- rep(NA,4)
	}else {
		if (site_Rrs["Rrs490"] > site_Rrs["Rrs560"]){
			tmp_tss <- QAA_560(site_Rrs)
		}else if (site_Rrs["Rrs490"] > site_Rrs["Rrs620"]){
			tmp_tss <- QAA_665(site_Rrs)
		}else if (site_Rrs["Rrs754"] > site_Rrs["Rrs490"] & site_Rrs["Rrs754"] > 0.010){
			tmp_tss <- QAA_865(site_Rrs)
		}else{
			tmp_tss <- QAA_754(site_Rrs)
		}
	}
	tmp_tss <- as.numeric(tmp_tss) 
	names(tmp_tss) <- c("a","bbp","band","TSS")
	#
	return(tmp_tss)
}
