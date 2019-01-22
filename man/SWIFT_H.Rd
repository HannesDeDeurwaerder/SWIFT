\name{SWIFT_H}
\alias{SWIFT_H}
\title{
Function to calculate the isotopic signatures at the users defined height and 
time}
\description{
This function calculates the isotopic signature at by the user defined height of measurement and timing of sampling. It is second part of the SWIFT model, and 
should be fed with the output of the SWIFT_SB output}
\usage{	SWIFT_H(Ax = NULL, D2Hxylem = NULL, hom = NULL, SF = NULL, 
tstud = NULL, tF = NULL)}
\arguments{
  \item{Ax}{	
  Description: The soil surface area covering the roots [in m2] 
	Structure: One value, representative of the studied plant}
  \item{D2Hxylem}{	
  Description: The deuterium signature at stem base - derived from the function
  SWIFT_SB [in permil, VSMOW]
	Structure: A discrete vector of n elements, with n the number of timesteps 
	sampled}
  \item{hom}{
  Description: A vector of studied measurement heights [in m]
	Structure: The size of the vector is defined by the interest of the user}
  \item{SF}{
  Description: Instantaneous sap flow over time [in kg s-1] 
	Structure: A vector of n elements, where n corresponds to the number of 
	timesteps}
  \item{tstud}{ 
  Description: By the user defined timesteps to which an output is required 
  [in second after initialisation]
	Structure: A vector with n elements defined by the user. 
	Note: When an output is needed for all values, then tstud = t}
  \item{tF}{
  Description: Measurement frequence per hour [in measurement points per hour] 
	Structure: One value, defined by the user
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
    			hts <- length(as.vector(hom))
    			tms <- length(as.vector(tstud))
    			D2Hxylem_hom <- matrix(NA, hts, tms)
    			uc <- ((60/tF) * 60)/1000
    			for (ab in 1:tms) {
        			a = tstud[ab]
       				   for (w in 1:hts) {
            				if (a == 1) {
						D2Hxylem_hom[w, ab] <- NA
                				next}

            				if (hom[w] == 0 && a != 1) {
                				D2Hxylem_hom[w, ab] <- D2Hxylem[a]
                				next}
            				else {
                				delay = 0
               					hloop <- hom[w]
                				if (SF[a - delay] * uc >= (hloop * Ax)) {
                  					D2Hxylem_hom[w, ab] <- D2Hxylem[a - delay]
                  					next}

                					while ((SF[a - delay] * uc) < (hloop * Ax) && 
                					(a - delay) > 0) {
						        	      hloop <- hloop - ((SF[a - delay] * uc)/(Ax))
                  					delay <- delay + 1
                					}

					                if (a - delay <= 0) {
						                 D2Hxylem_hom[w, ab] <- NA}
                					else {
						                 D2Hxylem_hom[w, ab] <- D2Hxylem[a - delay]
                					}
            					}
        				}
    			}
    			return(D2Hxylem_hom)
  		}
}