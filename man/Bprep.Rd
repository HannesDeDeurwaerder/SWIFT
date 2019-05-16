\name{Bprep}
\alias{Bprep}
\title{Calculate the root length distribution with depth}
\description{ 
This function constructs a vector containing the root length distribution of 
a plant at every defined soil layer [in m m-3]}
\usage{ Bprep(beta = NULL, R0 = NULL, Z = NULL)}
\arguments{
  \item{beta}{
  Description: Factor in the root length density distribution function; 
  Structure: One value, derived from Jackson et al. (1996), table 1}
  \item{R0}{
  Description: Factor in the root lenght density distribution function; 
  Structure: One value, representative of the studied plant, here derived from 
  Huang et al, [2017}. 
  NOTE: Provide a positive value !
  \item{Z}{	
  Description: A vector of soil depths [in m]; 		
  Structure: A discrete vector  of n elements, where n corresponds to the number
  of soil layers}
}
\details{}
\value{
A discrete vector of n elements, with n the number of soil layers, containing 
the root lenght distribution of a studied plant}
\references{
References for R0 and beta values: 
(a) Huang et al (2017). The effect of plant water storage on water fluxes within 
the coupled soil plant system. New Phytologist, 213(3), 1093-1106.
(b) Jackson et al (1996). A global analysis of root distributions for 
terrestrial biomes. Oecologia, 108(3), 389-411.
}
\author{ Hannes De Deurwaerder, Marco D. Visser, Matteo Detto, Pascal Boeckx, 
Felicien Meunier and Hans Verbeeck}
\note{}
\seealso{}
\examples{
## The function is currently defined as
function (beta = NULL, R0 = NULL, Z = NULL) {
    		  B <- rep(NA, length(Z))
    		  UnCo <- 100
    		  for (h in 1:length(Z)) {
    	    	    B[h] = R0 * beta^(UnCo * Z[h]) * log(beta, base=exp(1))
    		  }
    		return(B)
  		}
}