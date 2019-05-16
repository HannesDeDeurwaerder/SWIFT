\name{SWIFT_H}
\alias{SWIFT_H}
\title{
Calculate the isotopic signatures at the users defined height and time}
\description{
This function derives the isotopic signature at by the user defined height of 
measurement and timing of sampling. It is second part of the SWIFT model, and 
should be fed with the output of the SWIFT_SB output}
\usage{	SWIFT_H(Ax = NULL, D2Hxylem = NULL, hom = NULL, SF = NULL, 
tstud = NULL, tF = NULL)}
\arguments{
  \item{Ax}{	
  Description: The Xylem lumen area [in m2];  
  Structure: One value, representative of the studied plant}
  \item{D2Hxylem}{	
  Description: The deuterium signature at stem base - derived from the function
  SWIFT_SB [in permil, VSMOW]; 
  Structure: A discrete vector of n elements, with n the number of timesteps 
  sampled}
  \item{hom}{
  Description: A vector of studied measurement heights [in m]; 
  Structure: The size of the vector is defined by the interest of the user}
  \item{SF}{ 	
  Description: Instantaneous sap flow over time [in m3 s-1];  
  Structure: A vector of n elements, where n corresponds to the number of 
  timesteps}
  \item{tstud}{ 
  Description: By the user defined timesteps to which an output is required 
  [in seconds after initialisation]; 
  Structure: A vector with n elements defined by the user;  
  Note: When an output is needed for all values, then tstud = t}
  \item{tF}{
  Description: Measurement frequence per hour [in measurement points per hour];  
  Structure: One value, defined by the user; 
  Note: Lower numbers will increase the speed of the model, but for graphical 
  display, higher values are prefered}
}
\details{}
\value{
A vector or matrix of isotopic signatures at by the user defined height(s) and 
time(s) [in permil, VSMOW]}
\references{}
\author{ Hannes De Deurwaerder, Marco D. Visser, Matteo Detto, Pascal Boeckx, 
Felicien Meunier and Hans Verbeeck}
\note{}
\seealso{}
\examples{
	## The function is currently defined as
		function (Ax = NULL, D2Hxylem = NULL, hom = NULL, SF = NULL, tstud = NULL, 
		tF = NULL){
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
}
