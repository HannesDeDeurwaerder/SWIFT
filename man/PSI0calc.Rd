\name{PSI0calc}
\alias{PSI0calc}
\title{
Calculate the water potential at the stem base of a plant through time}
\description{
This function constructs a vector of water potentials at stem base for every 
timestep [in m H2O]}
\usage{ PSI0calc(ARi = NULL, k = NULL, PSIs = NULL, SF=NULL, 
t = NULL, Z = NULL)}
\arguments{
  \item{ARi}{
  Description: the fine root surface area distributed in soil layer i [in m2];  
  Structure: A vector of values with length number of soil layers, 
  representative of the studied plant}
  \item{k}{ 
  Description: Plant specific total soil-to-root conductance at each specific 
  soil layer [in s-1];  
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
  Structure: A discrete vector  of n elements, where n corresponds to the 
  number of soil layers}
}
\details{}
\value{ A vector of water potential at stem base for every timestep [in m H2O]}
\references{}
\author{
Hannes De Deurwaerder, Marco D. Visser, Matteo Detto, Pascal Boeckx, 
Felicien Meunier and Hans Verbeeck}
\note{}
\seealso{}
\examples{
## The function is currently defined as
		  function (ARi = NULL,  k = NULL, PSIs = NULL, SF = NULL, 
		  t = NULL, Z = NUL){
    	  PSI0vec=rep(NA,length(t))
        for (a in 1:length(t)){
          # Water potential at stem base
          PSI0 <- (sum(k*ARi*(PSIs-Z)) - SF[a])/ sum(k*ARi)
          RWU0 <- k*ARi*(PSIs-Z - PSI0)
          # make sure no negative RWU flows occurs
          use <- which(RWU0<0)
          RWU <- RWU0
          while (sum(use)>0){
              index <- which(RWU == min(RWU))
              RWU[index] <- 0
              use <- which(RWU>0)
              PSI0 <- (sum(k[use]*ARi[use]*(PSIs[use]-Z[use])) 
                        - SF[a]) / sum(k[use]*ARi[use])
              RWU[use] <- k[use]*ARi[use]*(PSIs[use]-Z[use] - PSI0)
              use <- which(RWU<0)
          }
          # assign PSI0 of one time point to the vector of PSI 
          # at stem base
          PSI0vec[a] <- PSI0
        }
      return(PSI0vec)
}
}