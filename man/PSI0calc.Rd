\name{PSI0calc}
\alias{PSI0calc}
\title{
Function to calculate water potential at the stem base of a plant through time}
\description{
This function constructs a vector of water potentials at stem base for every 
timestep [in m H2O]}
\usage{ PSI0calc(As = NULL, B = NULL, dZ = NULL, k = NULL, PSIs = NULL, 
r = NULL, rho = NULL, t = NULL, Z = NULL)}
\arguments{
  \item{As}{
  Description: The soil surface area covering the roots [in m2] 
		Structure: One value, representative of the studied plant}
  \item{B}{
  Description: Root length density at every discrete soil layer - Results from 
  the function Bprep - [in m m-3]
	Structure: A discrete vector of n elements, with n the number of soil layers} 
  \item{dZ}{
  Description: Thickness of the soil layer [in m]  
	Structure: One value, Samples should be taken at discrete depths which have
	the same thickness}
  \item{k}{ 
  Description: Plant specific total soil-to-root conductance at each specific 
  soil layer [in s-1] 
	Structure: A vector - results from the function SoilRootCond - and has n 
	elements corresponding to the number of soil layers}
  \item{PSIs}{
  Description: Soil water potential at the each specific soil layer [in m] 
	Structure: A vector of n element, where n corresponds to the number of soil 
	layers 
	Note: The water potentials per soil layer is currently considered static in 
	time}
  \item{r}{
  Description: Effective fine root radius [in m] 
	Structure: One value,representative of the studied plant}
  \item{rho}{
  Description: Density of water [kg m-3] 
	Structure: One constant value, i.e. 1000 kg m-3.}
  \item{t}{ 
  Description: Cumulative time steps [in s] 
	Structure: A vector}
  \item{Z}{ 
  Description: A vector of Soil depth [in m] 
	Structure: A discrete vector  of n elements, where n corresponds to the 
	number of soil layers}
}

\details{}
\value{ A vector of water potential at stem base for every timestep  [in m H2O]}
\references{}
\author{
Hannes De Deurwaerder, Marco D. Visser, Matteo Detto, Pascal Boeckx, 
Felicien Meunier and Hans Verbeeck}
\note{}
\seealso{}
\examples{
## The function is currently defined as
		function (As = NULL, B = NULL, dZ = NULL, k = NULL, PSIs = NULL, r = NULL, 
		rho = NULL, t = NULL, Z = NULL){
    			PSI0vec = rep(NA, length(t))
    			for (a in 1:length(t)) {
        		   PSI0vec[a] <-((sum(k * (Z + PSIs) * B)) - 
        		                ((SF[a])/(2 * pi * rho * r * As * dZ)))/sum(k * B)
    			}
    		return(PSI0vec)
  		}
}
