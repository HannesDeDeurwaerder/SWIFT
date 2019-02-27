SWIFT (0.1) [Release notes](https://github.com/HannesDeDeurwaerder/SWIFT/).
----------

The SWIFT model provides multiple functions for tracking stable water ([2H](https://en.wikipedia.org/wiki/Deuterium) and [18O](https://en.wikipedia.org/wiki/Isotopes_of_oxygen)) isotopic fluctuations and variance along the lenght of a plant. Functions ***SWIFT _ SB*** and ***SWIFT _ H*** can respectively be used to calculate the isotopic signature over time at the stem base of the plant, or at a by the user defined height and time. Additional functions are ***SoilRootCond***, ***PSI0calc*** and ***Bprep*** which in this order help the user to define *(i)* the soil to root conductivity for every defined soil layer, *(ii)* the water potential at stem base for every timestep, and *(iii)* the root length distribution for every defined soil layer.



### Quicklinks

-   [Quick start and tutorials](#quick-start-and-tutorials)
-   [Installation](#installation)
-   [Dependencies](#dependencies)
-   [Examples](#examples)
-   [Thanks](#thanks)
  

## Quick start and tutorials

An in depth description of the formula used within the SWIFT model can be found in **De Deurwaerder et al (In Review).** The example provided, in combination with the provided comments in the R-script should enable the user to succesfully run the model. 


## Installation

Use the **devtools** package to install the current development version from R.

### Install from Github repository

	# PACKAGE DEVTOOLS REQUIRED
	#--------------------------
	require(devtools)

	# INSTALL PACKAGE FROM GITHUB REPOSITORY
	#---------------------------------------
	install_github("HannesDeDeurwaerder/SWIFT")
	
	
### Dependencies

SWIFT has no other dependencies other than the base R installation.


## Examples

The code below defines a simple example of the use of the SWIFT model, based on the example provided in the paper of **De Deurwaerder et al, (In Review)** <Link will be added upon acceptance>


### Run the SWIFT model

	# PARAMETERS INITIALISATION
	#--------------------------
	n   <- 20     # Multiple number of days studied, needed for	spin up of the model
	tF  <- 60     # Time frequence of measurements per hour [in measurments per h] 
	t   <- seq(0,24*n,length.out = 24*tF*n)     # Discrete time vector [in h]
	dZ  <- 0.001    # Thickness of sampled layer [in m]	
	rho <- 1000     # Density of water [kg m^(-3)]
	L   <- 1        # maximum soil depth [in m]
	Z   <- seq(dZ,L,dZ)    # Discrete depth vector centered [in m]
	kr  <- 10*10^(-10) 	   # root membrane permeability [in s-1] 
						   # (source: Leuschner et al, 2004)
	As  <- 1      # Soil surface area covering the roots [in m^2]
	r   <- 0.0005 # Effective root radius [in m]		
	DBH <- 0.213  # Diameter at breast height [in m]
	LA  <- 0.136  # Lumen Fraction, i.e. lumen area/sapwood area [%]
				  # (derived from Zanne et al, 2010) [table 2, F-value]            
	Ax  <- LA * ((1.582*((DBH*100)^1.764)) /10^4)    
				# Total lumen area  of the studied tree [m²],i.e. the lumen area 
                # fraction multiplied by sapwood area estimated from the DBH 
 				# (Meinzer et al, 2001)  
	R0  <- -438688      # Derived from Huang et al, 2017
	beta<- 0.966        # Derived from Jackson et al., 1996 [table 1, beta-term]

		
	# SOURCE THE SWIFT MODEL
	#-----------------------
  	require(SWIFT)

  	# PREPARE ADITIONAL DATASETS
	#---------------------------
		# a. Diurnal sapflow pattern
		
			data(SFday) 	# SF for one day, time frequency: every min, expressed in
					# [kg h-1]. This data is derived from Huang et al, 2017
     			uch <- 60*60   	# unit conversion: h to sec
			SF  <- c( rep(SFday,n) )/(uch)		# repetition of SF day over n prefered days

		# b. Soil water potential curve with depth, from Meissner et al, 2012
			CTpsi <- 101.97		# Conversion factor between MPa and m H2O
			Apsi  <- 198.4455;  Bpsi <- 448.9092;  Cpsi <- 255.5937;
								# Apsi, Bpsi and Cpsi are values derived from the Meissner
								# et al 2012 data from optimizing the function defined in
								# the paper of De Deurwaerder et al. (In Review)
			PSIs  <- ((Apsi + Bpsi *log(Z) - Cpsi*Z^2) /10000) * CTpsi

		# c.Isotopic water isotopic signature with depth from Meissner et al, 2012
			lsoil   <- -73.98008;  msoil <- 0.148735;  soildev <- 0.001 	
					# lsoil and msoil are values derived from the Meissner et al 2012 
					# data from optimizing the function defined in the paper of 
					# De Deurwaerder et al. (In Review) 
			D2Hsoil <- lsoil * (Z+soildev)^msoil


	# REMARK
	#-------
	# NOTE that the data of Meissner et al 2012 is accessible and 
	# can be sourced as provided below:
	data(MeissnerData)


	# USING THE SWIFT MODEL FUNCTIONS
	#--------------------------------
		# a. Root density distribution 
        	B <- Bprep(beta, R0, Z) 		# output is in [m m^-3]
        
		# b. The soil-root conductance
        	k <- SoilRootCond(B, kr, PSIs, Z, 'Silt Loam')
        
		# c. Isotopic signature fluctuation at stem base
        	StemBase <- SWIFT_SB(As, B, D2Hsoil, dZ, k,  PSIs,  r, rho, SF, t, Z)
        
		# d. Isotopic signature at specific height and time
        	tstud <- seq(2*(24*tF) , 3*(24*tF) ,1)  
					# Provides the data of day 3 of the modeled tree 
					# (We selected day 3 to assure proper spin up of the model). 
        	hom   <- 1.3  	# measured at standard coring height, 1.3 m
        	D2Htree <- SWIFT_H( Ax, StemBase,  hom, SF, tstud, tF)
        
		# e. The water potential at stem base
        	PSI0vec <- PSI0calc(As, B, dZ, k,  PSIs, r, rho, t, Z)/CTpsi 	# CTpsi to convert from m to MPa     



	# MAKING A SIMPLE PLOT OF THE SWIFT OUTPUT
	#-----------------------------------------

	plot(t,StemBase, xlim= c( min(tstud)/tF, max(tstud)/tF ), type='l', lty=2,
	 	 col='grey', ylab='d2H (permill, VSMOW)', xlab='timesteps')	
		# Plot the signature fluctuations at the stem base
	
	lines( tstud/tF, D2Htree, col='blue')	
		# Plot the signature fluctuations at the height defined by the user

	legend( 'bottomleft', c('at stembase',paste0('at ',hom,' m')), lty=c(2,1), 
			col=c('grey', 'blue'), bty='n', cex=0.7)	
		# Add a legend 
	 


## Thanks
Special thanks to Marco D. Visser for all his appreciated help in creating this SWIFT package. I also thank all other co-authors for their help in develloping the SWIFT model.

