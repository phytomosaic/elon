######################################################################
# Wavelets analysis to calculate coherence of two timeseries --
#   examine effects of changing phase, frequency, amplitude
# Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 02 Dec 2019
#   CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)
######################################################################

### package load
install.packages('biwavelet')
require(biwavelet)
citation('biwavelet')
?wtc  # wavelet coherence (analogous to 'correlation')
?pwtc # partial analysis also available

### generate timeseries data,
###     to examine effects of changing phase, frequency, amplitude
set.seed(23)
n  <- 1000         # n sampling points
z  <- 1:n/10       # time vector
z1 <- z[1:250]     # break response vector into four segments
z2 <- z[251:500]
z3 <- z[501:750]
z4 <- z[751:1000]
xn <- rnorm(n,0,1) / 5 # noise vector
yn <- rnorm(n,0,1) / 5 # noise vector

### dual plot function: 2 timeseries above, their coherence below
`dual` <- function(x, y, w, ...){
     par(mfrow=c(2,1),oma=c(0,0,0,0),mar=c(2,4,0.5,1),las=1,bty='L')
     plot(z, x, type='l', lwd=0.2, col=1, xaxs='i')
     lines(z, y, type='l', lwd=0.2, col=2)
     plot(w, plot.phase=TRUE, arrow.len=0.1)
}

### effect of phase-shift
x <- sin(z) + xn
y <- c(sin(z1),sin(z2-pi/4),sin(z3-pi/2),sin(z4-pi)) + yn
w <- wtc(cbind(z,x), cbind(z,y), nrands=9) # nrands >100 in practice
dual(x,y,w)

### effect of frequency
x <- c(sin(1*z1), sin(2*z2), sin(4*z3), sin(8*z4)) + xn
y <- c(sin(1*z1), sin(2*z2), sin(4*z3), sin(8*z4)) + yn
w <- wtc(cbind(z,x), cbind(z,y), nrands=9) # nrands >100 in practice
dual(x,y,w)

### effect of amplitude
x <- c(sin(1*z1), sin(2*z2), sin(4*z3), sin(8*z4)) + xn
y <- c(sin(1*z1), sin(2*z2)/2, sin(4*z3)/3, sin(8*z4)/4) + yn
w <- wtc(cbind(z,x), cbind(z,y), nrands=9) # nrands >100 in practice
dual(x,y,w)

####   END   ####
