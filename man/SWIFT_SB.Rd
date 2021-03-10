\name{SWIFT_SB}
\alias{SWIFT_SB}
\title{
Calculate the isotopic signature at the stem base of a plant through time}
\description{	
This function derives the isotopic signature at by at the stem base of a 
plant through time. This function needs to be fed with data from the 
functions SoilRootCond and Bprep}
\usage{	SWIFT_SB(As = NULL, B = NULL, D2Hsoil = NULL, dZ = NULL, k = NULL, 
PSIs = NULL, r = NULL, rho = NULL, t = NULL, Z = NULL)}
\arguments{
  \item{ARi}{
  Description: the fine root surface area distributed in soil layer i [in m2];  
  Structure: A vector of values with length number of soil layers, 
  representative of the studied plant}
  \item{D2Hsoil}{	
  Description: Vector representing the stable water isotopic siganture profile 
  with depth [in permil, VSMOW]; 
  Structure: A discrete vector of n elements, with n the number of soil layers}
  \item{k}{ 
  Description: Plant specific total soil-to-root conductance at each specific
  specific soil layer [in s-1];  
  Structure: A vector - results from the function SoilRootCond - and has n 
  elements corresponding to the number of soil layers}
  \item{PSIs}{ 	
  Description: Soil water potential at the each specific soil layer [in m]; 
  Structure: A vector of n element, where n corresponds to the number of soil 
  layers;  
  Note: The water potentials per soil layer is currently considered static in
  time}
  \item{SF}{ 	
  Description: Instantaneous sap flow over time [in m3 s-1];  
  Structure: A vector of n elements, where n corresponds to the number of 
  timesteps}
  \item{t}{ 	
  Description: Cumulative time steps [in s];  
  Structure: A vector}
  \item{Z}{ 	
  Description: A vector of Soil depth [in m];  
  Structure: A discrete vector  of n elements, where n corresponds to the number
  of soil layers}
}
\details{}
\value{	A vector of isotopic signatures at at stem base for every timestep 
[in permil, VSMOW]}
\references{}
\author{Hannes De Deurwaerder, Marco D. Visser, Matteo Detto, Pascal Boeckx,
Felicien Meunier and Hans Verbeeck}
\note{}
\seealso{}
\examples{
## The function is currently defined as
		function (ARi = NULL, D2Hsoil = NULL, k = NULL, 
		PSIs = NULL, SF = NULL, t = NULL, Z = NULL){
      # Declare empty vectors
      D2Hxylem=rep(NaN,length(t))    
      D2Hvec <- fi <- RWU <- rep(NaN,length(Z))        
      
      for (a in 2:length(t)){
        # first value is NaN due to model spin-up.
        
        if((SF[a]==0 & a!=1) | (is.nan(SF[a]))){
              D2Hxylem[a]<-ifelse(is.nan(D2Hxylem[a-1]),
                            NaN,D2Hxylem[a-1]) 
          # When SF=0, water is stagnant --> signature at t=0
          # equals signature at t=-1.
        }else{
            # Water potential at stem base
            PSI0 <- (sum(k*ARi*(PSIs-Z)) - SF[a])/ sum(k*ARi)
            RWU0 <- k*ARi*(PSIs-Z - PSI0)
                not.use <- which(RWU0<0)
                RWU <- RWU0
                
                while (sum(not.use)>0){
                  not.use <- which(RWU<0)
                  RWU[not.use] <- 0
                  use <- which(RWU>0)
          
                  ARiuse = ARi[use]
                  kuse = k[use]
                  Zuse = Z[use]
                  PSIsuse = PSIs[use]
          
                  PSI0 <- (sum(kuse*ARiuse*(PSIsuse-Zuse)) - SF[a]) / sum(kuse*ARiuse)
                  RWU[use] <- kuse*ARiuse*(PSIsuse-Zuse - PSI0)
          
                }
          fi <- RWU/sum(RWU, na.rm=TRUE)
          D2Hxylem[a] <- sum(fi*D2Hsoil, na.rm=TRUE)
        }    
      }
      return(D2Hxylem)
}
}