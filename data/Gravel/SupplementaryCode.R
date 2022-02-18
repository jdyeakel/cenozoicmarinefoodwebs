##############################################
##############################################
#
# R Code supplementing the paper 
# Inferring food web structure form predator-prey body size relationship
# by Gravel, Poisot, Albouy, Velez, Mouillot
# Methods in Ecology and Evolution
# PUTS THE VOLUME/ISSUE/PAGES HERE
# February 2013
# 
##############################################
##############################################
#

#### READ ME #################################
#
# The code has the following structure
# 1. Loading of an artificial dataset of i) predator-prey interactions, ii) body size distribution, iii) hypothetical species presence/absence lists
# 2. Definition of functions run for the analysis
# 3. Calculation of the observed parameters for the niche model
# 4. Inferrence of the parameters for all species in the metaweb
#
# Name of the parameters of the niche model
# n: niche position along the niche axis
# c: centroid of the niche
# r: range of the niche
#
##############################################


##############################################
# 1.  The mediterranean Application ####### 

# data_loading: please load the R data file
# This file is a R list object with 4 objects
# Diet_matrix : predator prey claibration matrix that is a compilation of two resolved food webs, the Catalan Sea and Corsica.
# BodySize_cal : a vector of species log body size represented in the previous matrix for the model calibration
# Body_Size_med : log Body size of 557 Mediterranean fish species 
# Mat_Cooc : Co-ocurence martix 

### Data loading
load("data.rdata")


##############################################
# 2. Useful functions


##############################################
# Get regression parameters
# Input arguments:
# Bprey = log10 biomass of the prey
# Bpred = log10 biomass of the predator
# Quartil = a vector of the inferior and the superio quartile c(0.03,0.97)
# Returns a list of regression objectis
# Requires the quantreg package

reg_fn = function(Bprey,Bpred,quartil) {

	library(quantreg)
	mean_reg = lm(Bprey~Bpred)			# For the n parameter
	qrsup = rq(Bprey~Bpred,tau = quartil[2])	# For the higher limit of the range
	qrinf = rq(Bprey~Bpred,tau = quartil[1])	# For the lower limit of the range

	return(list(mean_reg$coef,qrsup$coef,qrinf$coef))

}

##############################################
# Estimate the niche parameters for all species of a list 
# Input arguments: 
# pars = resulting parameters of the function reg_Niche
# Ball = list of body size
# Returns a matrix with four parameters for each species

get_pars_Niche = function(pars,Ball) {
	
	# Unwrap the input parameters
	mean_reg = pars[[1]]
	qrsup = pars[[2]]
	qrinf = pars[[3]]
	
	# Estimate parameters for the allometric relationships
	delta = mean_reg[2]
	b1 = mean_reg[1]
	b2 = delta	
	
	# Estimate the parameters for the niche model
	n = Ball						# The niche n
	c = b1 + b2*Ball				# The centroid c
	low = qrinf[1] + qrinf[2]*Ball	# The lower limit of the range
	high = qrsup[1] + qrsup[2]*Ball	# The higher limit of the range
	
	return(cbind(n,c,low,high))	
}

##############################################
# Transform the parameters into an interaction matrix (the metaweb)
# Input:
# n = vector of size S with the parameter n for each of the S species
# c = vector of size S with the parameter c for each of the S species
# low = vector of size S with the parameter low for each of the S species
# high = vector of size S with the parameter high for each of the S species
# Returns a SxS matrix with 0 indicating absence of a link and 1 indicating the presence of a link
# Predators on columns, preys on rows

L_fn = function(n,c,low,high) {
	
	S = length(n)   	
   	L = matrix(0,nr=S,nc=S)

	for(i in 1:S)
    	for(j in 1:S)
    	  if(n[j]>low[i] && n[j]<high[i]) L[j,i] = 1
    	  
  	return(L)	
		
}

##############################################
#  3. Calculate parameters for the artificial dataset

# Calculate the parameters
pars_reg = reg_fn(Bprey,Bpred,quartil = c(0.05,0.95))
pars_niche = get_pars_Niche(pars_reg,Ball)
		
# Calculate the metaweb structure
L = L_fn(pars_niche[,1],pars_niche[,2],pars_niche[,3],pars_niche[,4])

# Illustrate the interaction matrix


##############################################
#  4. Application 
#
# The diet matrix
Mat_diet <- data$diet_matrix

# The log body size
BodySize <- sort(data$BodySize_cal)
S <- length(BodySize)

### Prey/predator interaction calculation
# Formating the data as a vector of prey body size and a vector for corresponding predator body size
pairs <- NULL
	
	for(i in 1:S){
		for(j in 1:S){
		
			if( Mat_diet[i,j] == 1){pairs <- rbind(pairs,c(BodySize[i],BodySize[j]))} # end of if
		} # end of j
	} # end of i

Bprey <- pairs[,1]
Bpred <- pairs[,2]

### Get regression parameters
Param_reg <- reg_fn (Bprey,Bpred,quartil = c(0.03,0.97))

### The log body size for the regional pool of 557 fish species
BS_med  <- data$Body_Size_med

### Estimation of the 3 niche parameters for all species
pars_Niche <- get_pars_Niche(Param_reg, BS_med)

### Put a threshold size to take into account that small fishes in the data do not have any prey (they are planktonivores)
Treshold  <-  which(pars_Niche[,"n"]<=1.367915)
pars_Niche[,"low"][Treshold] <- 0
pars_Niche[,"high"][Treshold] <- 0
    
### Computation of the final interaction matrix
mat_med <- L_fn(pars_Niche[,1],pars_Niche[,2],pars_Niche[,3],pars_Niche[,4])

### Correction by the range overlap
mat_Cooc <-data$mat_Cooc
Cor_mat = (mat_med * mat_Cooc)

### A representation of the metaweb

tiff(filename = "Figure_2.tiff",width =7.5, height = 10, units = "cm", res=800,compression = "lzw")
par(mfrow=c(2,1), mar=c(1,4,3.5,4))
	image(x=c(1:length(BS_med)),y=c(1:length(BS_med)),z=mat_med,xlab = "", ylab = "Prey body size rank",cex.axis = 0.4,cex.lab = 0.5,col=c("Black","white"),mgp=c(1,0.1,0),tcl=-0.2,xlim=c(25,557),ylim=c(25,557))
	mtext(text="(a)",side=2,line=1,las=2,at=558,cex=0.6)

par(mar=c(3.5,4,1,4))
	image(x=c(1:length(BS_med)),y=c(1:length(BS_med)),z=Cor_mat,xlab = "Predator body size rank", ylab = "Prey body size rank",cex.axis = 0.4,cex.lab = 0.5,col=c("black","white"),mgp=c(1,0.1,0),tcl=-0.2, xlim=c(25,557),ylim=c(25,557))
	mtext(text="(b)",side=2,line=1,las=2,at=558,cex=0.6)
dev.off()
