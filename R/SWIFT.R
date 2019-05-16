#'##############################################################################
#'##############################################################################
#'
#' @title SWIFT model (Stable Water Isotope Fluctuation within Trees)
#' 
#' 
#' @param ARtot Description: The total fine root area area of the focal plant
#'                           [in m2];  
#'              Structure:   One value, representative of the studied plant. 
#' @param Ax    Description: The total lumen area [in m2];  
#'              Structure:   One value, representative of the studied plant.
#' @param B     Description: The comunity plant Root length density at every 
#'                           discrete soil layer 
#'                           - Results from the function Bprep - [in m m-3]; 
#'              Structure:   A discrete vector of n elements, with n the number 
#'                           of soil layers.
#' @param beta  Description: Factor in the root length density distribution 
#'                           function for the focal plant; 
#'              Structure:   One value, derived from Jackson et al. (1996), 
#'                           table 1.
#' @param betaCom  Description: Factor in the root length density distribution 
#'                           function of the plant community 
#'                           (required to distribute the root length over all 
#'                           soil layers) 
#'              Structure:   One value, derived from Jackson et al. (1996), 
#'                           table 1.
#' @param D2Hxylem  Description: The deuterium signature at stem base - derived 
#'                               from the function SWIFT_SB [in permil, VSMOW];  
#'                  Structure:   A discrete vector of n elements, with n the 
#'                               number of timesteps sampled.
#' @param D2Hsoil   Description: Vector representing the stable water isotopic 
#'                               siganture profile with depth 
#'                               [in permil, VSMOW]; 
#'                  Structure:   A discrete vector of n elements, with n the 
#'                               number of soil layers.
#' @param dZ    Description: Thickness of the soil layer [in m];  
#'              Structure:   One value, Samples should be taken at discrete 
#'                           depths which have the same thickness.
#' @param hom   Description: A vector of studied measurement heights [in m]; 
#'              Structure:   The size of the vector is defined by the interest 
#'                           of the user. 
#' @param k     Description: Plant specific total soil-to-root conductance at
#'                           each specific specific soil layer [in s-1];  
#'              Structure:   A vector - results from the function SoilRootCond
#'                           - and has n elements corresponding to the number of
#'                           soil layers.  
#' @param kr    Description: The root membrane permeability [s-1];  
#'              Structure:   One value, representative of the studied plant.
#' @param PSIs  Description: Soil water potential at the each specific soil 
#'                           layer [in m]; 
#'              Structure:   A vector of n element, where n corresponds to the 
#'                           number of soil layers;  
#'              Note:        The water potentials per soil layer is currently 
#'                           considered static in time.
#' @param R0    Description: Factor in the root lenght density distribution per 
#'                           unit of soil, for entire plant community;
#'              NOTE:        Provide a Positive value!!  
#'              Structure:   One value,representative of the studied plant, here
#'                           derived from Huang et al, 2017.
#' @param SF    Description: Instantaneous sap flow over time [in m^3 s-1];  
#'              Structure:   A vector of n elements, where n corresponds to the 
#'                           number of timesteps.
#' @param Soiltype  Description:  Type of soil layer (Sand; Loamy Sand; 
#'                                Sandy Loam; Silt Loam; Loam; Sandy Clay Loam;
#'                                Silty Clay Loam; Clay Loam;Sandy Clay;
#'                                Silty Clay; Clay); 
#'                  Structure:    A Character string. 
#' @param t     Description:   Cumulative time steps [in s]; 
#'              Structure:     A vector.
#' @param tstud Description:   By the user defined timesteps to which an output 
#'                             is required [in second after initialisation]; 
#'              Structure:     A vector with n elements defined by the user;  
#'              Note:          When an output is needed for all values, then 
#'                             tstud = t.
#' @param tF    Description:   Measurement frequence per hour [in measurement 
#'                             points per hour];  
#'              Structure:     One value, defined by the user; 
#'              Note:          Lower numbers will increase the speed of the 
#'                             model, but for graphical display, higher values 
#'                             are prefered.
#' @param Z     Description:   A vector of Soil depth [in m]; 
#'              Structure:     A discrete vector  of n elements, where n 
#'                             corresponds to the number of soil layers.
#'
#'
#'   
#' @author Hannes De Deurwaerder, Marco D. Visser, Matteo Detto, Pascal Boeckx, 
#' Felicien Meunier and Hans Verbeeck.
#' 
#' @examples
#' \dontrun{
#' 
#' # Initialisation of parameters
#'   n   <- 20            # Multiple number of days studied, needed for spin up 
#'                        # of the model.
#'   tF  <- 60            # Time frequence of measurements per hour 
#'                        # [in measurments per h].
#'   t   <- seq(0,24*n,length.out = 24*tF*n)     # Discrete time vector [in h].
#'   dZ  <- 0.001         # Thickness of sampled layer [in m].
#'   L   <- 1             # maximum soil depth [in m].
#'   Z   <- seq(dZ,L,dZ)  # Discrete depth vector centered [in m].
#'   kr  <- 10*10^(-10)   # root membrane permeability [in s^(-1)] 
#'                        # (source: Leuschner et al, 2004).
#'   DBH <- 0.213         # Diameter at breast height [in m].
#'   LA  <- 0.136         # Lumen Fraction, i.e. lumen area/sapwood area [%]
#'                        # (derived from Zanne et al, 2010) [table 2, F-value].            
#'   Ax  <- LA * ((1.582*((DBH*100)^1.764)) /10^4)       
#'                        # Total lumen area of the studied tree [m^2],i.e. the 
#'                        # lumena area fraction multiplied by sapwood area 
#'                        # estimated from the DBH (Meinzer et al, 2001) . 
#'   ARtot <- exp(0.88*log(pi*2*(100*DBH/2)^2)-2)  
#'                        # ARtot in [m^2] via DBH[in cm] (Cermak et al, 2006)
#'   R0  <- 438688       # Derived from Huang et al, 2017.
#'   betaCom <- 0.976     # Derived from Jackson et al., 1996 [table 1, 
#'                        # beta-term].
#'   beta <-  0.976       # Derived from Jackson et al., 1996 [table 1, 
#'                        # beta-term].
#' # Source the SWIFT model
#'   require("SWIFT")
#'
#' 
#' # Prepare additional input data of the model
#'   data(SFday)     
#'    # SF for one day, time frequency: every min, expressed in [kg h-1].
#'    # This data is derived from Huang et al, 2017.
#'     uch <- (60/tF)*60*1000              # unit conversion: h to sec. and kg to m3
#'     SF  <- c(rep(SFday,n))/(uch) # repetition of SF day over n prefered days.
#'    
#' # Soil water potential curve with depth, from Meissner et al, 2012
#'    CTpsi <- 101.97      # Conversion factor between MPa and m H2O.
#'    Apsi  <- 198.4455;  Bpsi <- 448.9092;  Cpsi <- 255.5937; 
#'      # Apsi, Bpsi and Cpsi are values derived from  Meissner et al 2012 data
#'      # from optimizing the function defined in the paper of De Deurwaerder 
#'      # et al.  (In Review).
#'    PSIs  <- ((Apsi + Bpsi *log(Z) - Cpsi*Z^2) /10000) * CTpsi
#'  
#' # isotopic water isotopic signature with depth from Meissner et al, 2012
#'    lsoil   <- -73.98008;  msoil <- 0.148735;  soildev <- 0.001  
#'      # lsoil and msoil are values derived from the Meissner et al 2012 data
#'      # from optimizing the function defined in the paper of De Deurwaerder 
#'      # et al. (In Review).
#'    D2Hsoil <- lsoil * (Z+soildev)^msoil
#'    
#'  # NOTE that the data of Meissner et al 2012 is accessible and can be sourced
#'  # as provided below:
#'    data(MeissnerData)
#'         
#'
#'  # Using the SWIFT model functions
#'   
#'        # A. Root density distribution 
#'        B <- Bprep(betaCom, R0, Z) # output is in [m m^(-3)].
#'        
#'        # B. the soil-root conductance
#'        k <- SoilRootCond(B, kr, PSIs, Z, 'Silt Loam')
#'        
#'        # C. the ARi
#'        devio=beta^(100*Z)*(1-beta^(100*dZ))
#'        ARi=ARtot*devio/sum(devio)
#'        
#'        # D. Isotopic signature fluctuation at stem base
#'        StemBase <- SWIFT_SB(ARi, D2Hsoil,  k,  PSIs,  SF, t, Z)
#'        
#'        # E. Isotopic signature at specific height and time
#'        tstud <- seq(2*(24*tF), 3*(24*tF),1)    
#'           # Provides the data of day 3 of the modeled tree. 
#'           # (We selected day 3 to assure proper spin up of the model).
#'        hom   <- 1.3  # measured at standard coring height, 1.3 m.
#'        D2Htree <- SWIFT_H( Ax, StemBase,  hom, SF, tstud, tF)
#'        
#'        # E. The Water potential at stem base
#'        PSI0vec <- PSI0calc(ARi,  k,  PSIs,  Sf, t, Z)/CTpsi  
#'           # CTpsi to convert from m to MPa.     
#'        
#'        
#'  # Field data corroborating the hypothesis of strong isotopic fluctuations 
#'  # along the length of a stem:
#'    data(Fielddata)
#'  
#'  
#'  # Making a simple plot of the SWIFT output
#'
#'   	# MAKING A SIMPLE PLOT OF THE SWIFT OUTPUT
#'       ylabel=expression(paste(delta,""^"2","H"," [","\211",", vsmow]"))
#'      plot(t,StemBase, xlim= c( min(tstud)/tF, max(tstud)/tF), type='l', 
#'      lty=2, col='grey', ylab=ylabel, xlab='timesteps', mgp = c(2, 0.8, 0))	
#'          # Plot the signature fluctuations at the stem base
#'      lines( tstud/tF, D2Htree, col='blue')	
#'          # Plot the signature fluctuations at the height defined by the user
#'      legend( 'bottomleft', c('at stembase',paste0('at ',hom,' m')), 
#'      lty=c(2,1), col=c('grey', 'blue'), bty='n', cex=0.7)	
#'                    
#' load()
#'	}
#'	
#' @return Provided by function 
#'         - SWIFT_SB:      Vector of isotopic signatures at the stem base for 
#'                          every timestep [in permil, VSMOW].
#'         - SWIFT_H:       Vector of isotopic signatures at by the user defined
#'                          height and time [in permil, VSMOW].
#'         - SoilRootCond:  Vector of soil to root conductivity for every 
#'                          defined soil layer [in s-1].
#'         - PSI0calc:      Vector of water potential at stem base for every 
#'                          timestep  [in m H2O].
#'         - Bprep:         Vector of root length distribution for every defined
#'                          soil layer [in m m-3].
#'         
#' @export
#' 
#'##############################################################################
#'##############################################################################


SWIFT_SB<-function(ARi=NULL, D2Hsoil=NULL, k=NULL, PSIs=NULL,  SF=NULL,  t=NULL, Z=NULL){
  

#===============================================================================
#                              SWIFT MODEL PART A 
#===============================================================================

  #    STABLE WATER ISOTOPES CALUCLATED AT STEM BASE
  #-------------------------------------------------------
  # Function description: This function calculates the isotopic signature at the
  # stem base at every timestep.
    
      # Declare empty vectors
      D2Hxylem=rep(NA,length(t))    
      D2Hvec <- fi <- RWU <- rep(NA,length(Z))        
      
    
      for (a in 1:length(t)){
        # Water potential at stem base
        PSI0 <- (sum(k*ARi*(PSIs-Z)) - SF[a])/ sum(k*ARi)
    
        if(a==1){D2Hxylem[a]<-NA}     # first value is NA due to model spin-up.
        if((SF[a]=0 & a!=1) | (is.nan(SF[a]) & a!=1) ){
          D2Hxylem[a]<-D2Hxylem[a-1] 
          # When SF=0, water is stagnant --> signature at t=0 equals signature 
          # at t=-1.
        
        }else{
          
          # Relative contribution of every soil layer
          PSIdelta <- PSI0-(PSIs-Z) 
          
          RWU <- k*ARi*PSIdelta 
          RWU[is.na(PSIdelta) | PSIdelta>=0] <- NA
          
          
          # Isotopic signature at every soil layer
          fi <- RWU/sum(RWU, na.rm=TRUE)
          D2Hvec <-  fi*D2Hsoil
          
          
          # Summation over all soil layers
          D2Hxylem[a] <- sum(D2Hvec, na.rm=TRUE)
          rm(D2Hvec);
        }    
      }
      return(D2Hxylem)
}

      
################################################################################
################################################################################

   
SWIFT_H<-function(Ax=NULL, D2Hxylem=NULL, hom=NULL, SF=NULL,  tstud=NULL,  
                  tF=NULL){
      
#===============================================================================
#                              SWIFT MODEL PART B 
#===============================================================================
      
      # STABLE WATER ISOTOPES CALUCLATED AT BY THE USER DEFINED HEIGHT AND TIME
      #-------------------------------------------------------------------------
      # Function description: This function calculates the isotopic signature at
      # by the user defined height of measurement and timing of sampling.
      
      hts <-  length( as.vector(hom) )   # vector of user defined heights.
      tms <-  length( as.vector(tstud))  # vector of user defined sample times.
      
      D2Hxylem_hom <- matrix(NA,hts,tms)
      uc <-((60/tF)*60) # unit conversion.
      
      # Cumulative SF vector
      cumSF=cumsum(SF) # Cumulative sapflow 
      cumH=uc*cumSF/Ax # Cumulative Height reached by the sapflow
      cumH[which(cumSF==0)]<-NaN
      # No sapflow yet, this is the spinup, should be set at NA
      
      
      for (ab in 1:tms){
        a=tstud[ab] # By user defined timesteps. 
        
      for (w in 1:hts){
      
        
        # First value will be NA, for spin
        if (a==1){
          D2Hxylem_hom[w,ab] <- NaN
          next
        }
        
        # for signatures at stem base, i.e. h=0
        if (hom[w]==0 && a!=1){D2Hxylem_hom[w,ab] <- D2Hxylem[a]
        next
        }
        
         
        # for all other cases
            # for the spinup values 
            if(is.nan(cumH[a])){D2Hxylem_hom[w,ab] <- D2Hxylem[a-1]
            next
            
            }else{
              # not spinup values
              CumHt <- cumH -(cumH[a])+ hom[w]
              Minind <-  which(abs(CumHt)==min(abs(CumHt), na.rm=TRUE)) 
              # Total delay for the sapflow to reach hom (time)
              
              # Delay extends the datarange?
              if (length(Minind)==0){
                D2Hxylem_hom[w,ab] <- NaN  # only a problem during model spin up
              }else{
                delay=a-Minind[1]
                D2Hxylem_hom[w,ab] <- D2Hxylem[a-delay]
              }  
              
              rm(delay); rm(Minind) 
              
            }
            
        
        }
      }
      return(D2Hxylem_hom)
      
}    


################################################################################
################################################################################

SoilRootCond <- function(B=NULL, kr=NULL, PSIs=NULL, Z=NULL, Soiltype=NULL ){

#===============================================================================
#                PLANT SPECIFIC TOTAL SOIL TO ROOT CONDUCTANCE  
#===============================================================================
# Function description: This function calculates the soil to root conductance 
# while considering the soiltype.


          # Soiltype conform the soil texture triangle, with corresponding 
          # values from clapp and hornberger (1978).
          if(Soiltype == 'Sand'){
            ksmax=1.056/(100*60);   PSIsat=-0.121;   b=4.05;    sigmasat=0.395}
          if(Soiltype == 'Loamy Sand'){
            ksmax=0.938/(100*60);   PSIsat=-0.090;   b=4.38;    sigmasat=0.410}
          if(Soiltype == 'Sandy Loam'){
            ksmax=0.208/(100*60);   PSIsat=-0.218;   b=4.90;    sigmasat=0.435}
          if(Soiltype == 'Silt Loam'){
            ksmax=0.0432/(100*60);  PSIsat=-0.786;   b=5.30;    sigmasat=0.485}
          if(Soiltype == 'Loam'){
            ksmax=0.0417/(100*60);  PSIsat=-0.478;   b=5.39;    sigmasat=0.451}
          if(Soiltype == 'Sandy Clay Loam'){
            ksmax=0.0378/(100*60);  PSIsat=-0.299;   b=7.12;    sigmasat=0.420}
          if(Soiltype == 'Silty Clay Loam'){
            ksmax=0.0102/(100*60);  PSIsat=-0.356;   b=7.75;    sigmasat=0.477}
          if(Soiltype == 'Clay Loam'){
            ksmax=0.0147/(100*60);  PSIsat=-0.630;   b=8.52;    sigmasat=0.476}
          if(Soiltype == 'Sandy Clay'){
            ksmax=0.0130/(100*60);  PSIsat=-0.153;   b=10.40;   sigmasat=0.426}
          if(Soiltype == 'Silty Clay'){
            ksmax=0.0062/(100*60);  PSIsat=-0.490;   b=10.40;   sigmasat=0.492}
          if(Soiltype == 'Clay'){
            ksmax=0.0077/(100*60);  PSIsat=-0.405;   b=11.40;   sigmasat=0.482}
      
        
       # calculate k using the Clapp and Hornberger (1978) equation. 
        k=ks<-rep(NA,length(Z))
    
        for (h in 1:length(Z)){
          ks[h] <-  ksmax * (PSIsat/PSIs[h])^(2+3/b)    
          k[h]<- kr * ks[h]*sqrt(pi*B[h]) / (0.53*kr + ks[h]*sqrt(pi*B[h]) )
    
        }
        
        return(k)
        
    }


################################################################################
################################################################################


PSI0calc <- function(ARi=NULL, k=NULL,  PSIs=NULL, SF=NULL, t=NULL, Z=NULL){
  
#===============================================================================
#                          WATER POTENTIAL AT STEM BASE  
#===============================================================================
# Function description: This function calculates the water potential at stem 
# base at every timestep, indirectly derived from the sap flow. 

      PSI0vec=rep(NA,length(t))
        for (a in 1:length(t)){
          PSI0vec[a]<-(sum(k*ARi*(PSIs-Z)) - SF[a])/ sum(k*ARi)
        }
        return(PSI0vec)
}



################################################################################
################################################################################


Bprep<-function( beta=NULL, R0=NULL, Z=NULL){
  
  
#===============================================================================
#                       ROOT LENGTH DISTRIBUTION FUNCTION  
#===============================================================================
# Function description: This function calculates root length distribution 
# cf. Huang et al, 2016.


      B <- rep(NA,length(Z)) 
      UnCo <- 100 # unit conversion from m to cm.
      
        for (h in 1:length(Z)){
          B[h]= -R0* beta^(UnCo*Z[h])*log(beta, base=exp(1))    
        }
      
      return(B)

}





################################################################################
################################################################################




#' The Laussat Trees and Liana isotopic signature along the lenght of 
#' the plant dataset.
#'
#'
#' General Description:
#' This dataset contains both sample extraction information, i.e.
#' timing and height of sampling, the diameter, growth form type and species 
#' name, as well as the isotopic deuterium and oxygen 18 signature of the 
#' sampled plant
#'
#'
#' This dataset may not be
#' used without previous written permission from the
#' owner.
#'
#'
#' Ownership: Hannes De Deurwaerder  <Hannes_de_deurwaerder@hotmail.com>
#' 
#' Research Location:
#' Within the research forest of LAUSSAT, French Guyana.
#' Sampled species all are located on a white sandy subsoil.
#'
#' Dataset abstract: 
#' These data are being used to evaluate  hypothesized mechanisms of
#' xylem water isotopic variation along the length of a plant 
#' 
#' 
#' Methods:
#' The methods of sampling are in much detail described in 
#' De Deurwaerder et al. (2018) Tree Physiology. 38(7). p.1071-1083
#' 
#' 
#' The dataset contains the following labels (columns):
#' \itemize{
#' \item Number (Integer) Unique identifier per growth form.
#' \item Type (Character) Indication of the growth form, i.e. liana or tree.
#' \item Family (Character) Identification of the scientific plant family.
#' \item Family (Character) Species (Character) identification of the scientific 
#' species name.
#' \item Date (Character) Date of sampling in the year 2017.
#' \item Time (Character) Time of sampling.
#' \item Hour (Integer) Hour of measurement.
#' \item Minute (Integer) Minute of measurement.
#' \item Height (Numeric) Height of measurement [in m].
#' \item H2 (Numeric) Xylem deuterium isotopic signature [in permil, VSMOW]. 
#' \item stdev_H2 (Numeric) Standard deviation on the xylem deuterium isotopic 
#' signature [in permil, VSMOW], resulting from the Picarro measurement.
#' \item O18 (Numeric) Xylem O18 isotopic signature [in permil, VSMOW]. 
#' \item stdev_O18 (Numeric)Standard deviation on the xylem O18 isotopic 
#' signature [in permil, VSMOW], resulting from the Picarro measurement.
#' \item Diam (Numeric) Diameter of the growth form [in cm]. 
#' \D1O
#' }
#'
#'
#' @docType data
#' @keywords datasets
#' @name Fielddata
#' @usage data(Fielddata)
#' @format 
#' 
NULL


################################################################################
################################################################################


#' The Meissner et al 2012 of  water isotopic signatures
#' and the water potential in the soil. 
#'
#' General Description:
#' This data, providing an overview of both deuterium and oxygen 18
#' signatures as well as the water potential gradient discussed in the paper 
#' of meissner et al, 2012 papers entitled: (i)Partitioning of soil water among 
#' canopy trees during a soil desiccation period in a temperate mixed forest' 
#' and (ii) Soil water uptake by trees using water stable isotopes 
#' (d2H and d18O)- a method testregarding soil moisture, texture and carbonate
#'
#'
#' This dataset is freely accessible from the papers of Meissner et al 2012.
#' However, users should feel enclined to cite Meissner et al 2012 when using 
#' the data
#'
#'  
#' Ownership: M. Meissner <mmeissn3@gwdg.de>
#' 
#' 
#' Research Location:
#' Soil samples were collected at the Hainich National Park close to the vilage 
#' of Weberstedt in Central Germany. Height of sampling was approx. 350 m a.s.l.
#' and the region experiences a subatlatic climate with a mean anual 
#' precipitation around 600 mm and st.dev. around 50mm. 
#' More details are provided in both articles of Meissner et al, 2012.
#' 
#' 
#' Dataset abstract: 
#' These data is being used to test the hypothesized mechanisms of xylem water 
#' isotopic variation along the length of a plant.
#'     
#' Methods:
#' The methods of sampling are in much detail described in both papers of 
#' Meissner et al, 2012. More specifically, this data forms the driver data of 
#' the SWIFT model (1) "Soil water uptake by trees using water stable isotopes 
#' (d2H and d18O) - a method testregarding soil moisture, texture and carbonate
#' (2) Partitioning of soil water among canopy trees during a soil
#' desiccation period in a temperate mixed forest"
#' 
#' 
#' The dataset contains the following labels (columns):
#' \itemize{
#' \item Meissner_depth (Numeric) The average depth of the samples soil layer 
#' [in m].
#' \item Meissner_d2h (Numeric) The deuterium isotopic siganture of the bulk 
#' soil water at the sampled depths [in permil, VMSOW].
#' \item Meissner_sdd2h (Numeric) The standard deviations on the deuterium 
#' isotopic siganture of the bulk soil 
#' water at the sampled depths [in permil, VMSOW].
#' \item Meissner_d18o(Numeric) The oxygen 18 isotopic siganture of the bulk 
#' soil water at the sampled depths [in permil, VMSOW].
#' \item Meissner_sdd18o (Numeric) The standard deviations on the oxygen 18 
#' isotopic siganture of the bulk soil 
#' water at the sampled depths [in permil, VMSOW].
#' \item Meissner_psi (Numeric) The water potential at the sampled soil depths 
#' [in hPa].
#' \item Meissner_sdpsi (Numeric) The standard deviations on the water potential
#'  at the sampled soil depths [in hPa].
#' \D1O
#' }
#'
#'
#' @docType data
#' @keywords datasets
#' @name MeissnerData
#' @usage data(MeissnerData)
#' @format 
#' 
NULL



################################################################################
################################################################################


#' The daily sap flow derived from Huang et al 2017.
#'
#'
#' General Description:
#' The daily sap flow flux [in kg h-1] derived from the Huang et al 2017 paper, 
#' day 11 of scenario 6. Data is interpollated to derive a timefrequency per 
#' minute.
#'
#'
#' This dataset is freely accessible from the paper of Huang et al. 2017
#' However, Users should feel enclined to cite Huang et al. 2017 when using the 
#' data. paper title: The effect of plant water storage on water fluxes within 
#' the coupled soil-plant system.
#'
#'  
#' Ownership: Cheng-Wei Huang <chengweihuang1206@gmail.com>
#'
#'  
#' Research Location:
#' Details on the research location and scenario description
#' is in detail described in Huang et al, 2017, New Phytologist paper
#' entitled: The effect of plant water storage on water fluxes within
#' the coupled soil-plant system
#'
#'
#' Dataset abstract: 
#' These data are being used to evaluate  hypothesized mechanisms of
#' xylem water isotopic variation along the length of a plant.
#' This data is used to drive the SWIFT model, as presented by De Deurwaerder
#' et al. (IN REVIEW) 
#' 
#' 
#' Methods:
#' Details on the sap flow data and scenario description
#' is in detail described in Huang et al, 2017, New Phytologist paper
#' entitled: The effect of plant water storage on water fluxes within
#' the coupled soil-plant system
#' 
#' 
#' The dataset contains the following labels (columns):
#' \itemize{
#'
#' \item SFday (Numeric) sap flow [kg h-1] per time step frequency of 1 minute.
#' \item timestep (Numeric) Timesteps throughout the duration of an entire day, 
#' i.e. 24h
#' \D1O
#' }
#'
#'
#' @docType data
#' @keywords datasets
#' @name SFday
#' @usage data(SFday)
#' @format 
#' 
NULL
      
