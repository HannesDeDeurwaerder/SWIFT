\name{SoilRootCond}
\alias{SoilRootCond}
\title{
Calculate the soil to root conductivity for every defined soil layer}
\description{
This function constructs a vector of soil to root conductivity for every defined
soil layer [in s-1], while considering the soiltype}
\usage{ SoilRootCond(B = NULL, kr = NULL, PSIs = NULL, Z = NULL, 
Soiltype = NULL)}
\arguments{
  \item{B}{
  Description: Root length density at every discrete soil layer - Results from 
  the function Bprep - [in m m-3]; 
  Structure: A discrete vector of n elements, with n the number of soil layers} 
  \item{kr}{
  Description: The root membrane permeability [s-1]; 
  Structure: One value, representative of the studied plant}
  \item{PSIs}{
  Description: Soil water potential at the each specific soil layer [in m];  
  Structure: A vector of n element, where n corresponds to the number of soil 
  layers;  
  Note: The water potentials per soil layer is currently considered static in 
  time}
  \item{Z}{
  Description: A vector of Soil depth [in m];  
  Structure: A discrete vector of n elements, where n corresponds to the number 
  of soil layers}
  \item{Soiltype}{  
  Description: Soil Tecture class studied (Sand; Loamy Sand; Sandy Loam; 
  Silt Loam; Loam; Sandy Clay Loam; Silty Clay Loam; Clay Loam;Sandy Clay; 
  Silty Clay; Clay); 
  Structure: A Character string}
}
\details{}
\value{	
A vector of soil to root conductivity for every defined soil layer [in s-1]}
\references{
Soil texture characteristics are from:
(a) Clapp and Hornberger (1978). Empirical equations for some soil 
hydraulic propert}
\author{	
Hannes De Deurwaerder, Marco D. Visser, Matteo Detto, Pascal Boeckx, 
Felicien Meunier and Hans Verbeeck}
\note{}
\seealso{}
\examples{
		## The function is currently defined as
		function (B = NULL, kr = NULL, PSIs = NULL, Z = NULL, Soiltype = NULL){
    			if (Soiltype == "Sand") {
        			ksmax = 1.056/(100 * 60)
        			PSIsat = -0.121
        			b = 4.05
        			sigmasat = 0.395}

    			if (Soiltype == "Loamy Sand") {
        			ksmax = 0.938/(100 * 60)
        			PSIsat = -0.09
        			b = 4.38
        			sigmasat = 0.41}

    			if (Soiltype == "Sandy Loam") {
        			ksmax = 0.208/(100 * 60)
        			PSIsat = -0.218
        			b = 4.9
        			sigmasat = 0.435}

    			if (Soiltype == "Silt Loam") {
 			       ksmax = 0.0432/(100 * 60)
 			       PSIsat = -0.786
  			      	b = 5.3
   			     	sigmasat = 0.485}

    			if (Soiltype == "Loam") {
       			 	ksmax = 0.0417/(100 * 60)
       			 	PSIsat = -0.478
        			b = 5.39
        			sigmasat = 0.451}

    			if (Soiltype == "Sandy Clay Loam") {
        			ksmax = 0.0378/(100 * 60)
        			PSIsat = -0.299
        			b = 7.12
        			sigmasat = 0.42}

    			if (Soiltype == "Silty Clay Loam") {
        			ksmax = 0.0102/(100 * 60)
        			PSIsat = -0.356
        			b = 7.75
        			sigmasat = 0.477}

    			if (Soiltype == "Clay Loam") {
        			ksmax = 0.0147/(100 * 60)
        			PSIsat = -0.63
        			b = 8.52
        			sigmasat = 0.476}
    
			    if (Soiltype == "Sandy Clay") {
        			ksmax = 0.013/(100 * 60)
        			PSIsat = -0.153
        			b = 10.4
        			sigmasat = 0.426}

    	    if (Soiltype == "Silty Clay") {
        			ksmax = 0.0062/(100 * 60)
        			PSIsat = -0.49
        			b = 10.4
        			sigmasat = 0.492}

			    if (Soiltype == "Clay") {
        			ksmax = 0.0077/(100 * 60)
        			PSIsat = -0.405
        			b = 11.4
        			sigmasat = 0.482}

    			k = ks <- rep(NA, length(Z))
            for (h in 1:length(Z)){
                ks[h] <-  ksmax * (PSIsat/PSIs[h])^(2+3/b)    
                k[h]<- kr * ks[h]*sqrt(pi*B[h]) / (0.53*kr + ks[h]*sqrt(pi*B[h]) )
          
              }
    			return(k)
  			}
}
