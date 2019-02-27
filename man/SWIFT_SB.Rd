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
  \item{As}{
  Description: The soil surface area covering the roots [in m2];  
  Structure: One value, representative of the studied plant}
  \item{B}{ 
  Description: Root length density at every discrete soil layer - Results from 
  the function Bprep - [in m m-3]; 
	Structure: A discrete vector of n elements, with n the number of soil layers} 
  \item{D2Hsoil}{	
  Description: Vector representing the stable water isotopic siganture profile 
  with depth [in permil, VSMOW]; 
	Structure: A discrete vector of n elements, with n the number of soil layers}
  \item{dZ}{ 
  Description: Thickness of the soil layer [in m];   
	Structure: One value, Samples should be taken at discrete depths which have 
	the same thickness}
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
  \item{r}{ 	
  Description: Effective fine root radius [in m];  
	Structure: One value,representative of the studied plant}
  \item{rho}{ 	
  Description: Density of water [kg m-3];  
	Structure: One constant value, i.e. 1000 kg m-3.}
	\item{SF}{
  Description: Instantaneous sap flow over time [in kg s-1];  
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
		function (As = NULL, B = NULL, D2Hsoil = NULL, dZ = NULL, k = NULL, 
		PSIs = NULL, r = NULL, rho = NULL, SF = NULL, t = NULL, Z = NULL){
    			D2Hxylem = rep(NA, length(t))
    			D2Hvec <- fi <- Qi <- rep(NA, length(Z))
   
			for (a in 1:length(t)) {
        			PSI0 <- ((sum(k * (Z + PSIs) * B)) - ((SF[a])/
        			(2 * pi * rho * r * As * dZ)))/sum(k * B)
        	
				if (a == 1) { D2Hxylem[a] <- NA }

        			if ((SF[a] = 0 & a != 1) | (is.nan(SF[a]) & a != 1)) {
            				D2Hxylem[a] <- D2Hxylem[a - 1] }
        			else {
            				PSIdelta <- PSI0 - PSIs
            				Qi <- PSIdelta * B * k
            				Qi[is.na(PSIdelta) | PSIdelta >= 0] <- NA

            				fi <- Qi/sum(Qi, na.rm = TRUE)
            				D2Hvec <- fi * D2Hsoil

            				D2Hxylem[a] <- sum(D2Hvec, na.rm = TRUE)
            				rm(D2Hvec)
            				Qz <- NA
        			}
    			}
    		return(D2Hxylem)
  		}
}