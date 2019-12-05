###  Messing around to examine cor btwn two timeseries of
###       lichens and climate at HJA
###  Rob Smith, smithr2@oregonstate.edu, 16 Dec 2017
###

require(biwavelet)
require(zoo)

z1 <- enviro.data  # from the biwavelet pkg
z1 <- zoo(z1[,4:6], z1[,3]) # convert to zoo
plot(z1, screen=1, col=c(1,2,3))
t1 <- cbind(index(z1), as.numeric(z1[,1])) # convert back to matrix
t2 <- cbind(index(z1), as.numeric(z1[,2]))
t3 <- cbind(index(z1), as.numeric(z1[,3]))

t1 <- t3
# # two timeseries
# t1 <- cbind(1:100, rnorm(100))
# t2 <- cbind(1:100, rnorm(100))
# z1 <- zoo(cbind(t1[,2],t2[,2]), t1[,1])
# plot(z1, screen=1, col=c(1,3))

# Continuous wavelet transform: compute wavelet spectra
wt.t1 <- wt(t1)
wt.t2 <- wt(t2)

# Compare non-corrected vs. corrected wavelet spectrum
op <- par(mfrow=c(1,2))
plot(wt.t1, type="power.corr.norm", main="Bias-corrected")
plot(wt.t1, type="power.norm", main="Not-corrected")
par(op)

# Plot power
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
plot(wt.t1, plot.cb=TRUE, plot.phase=FALSE)

# Compute cross-wavelet
xwt.t1t2 <- xwt(t1, t2)

# Plot cross wavelet power and phase difference (arrows)
plot(xwt.t1t2, plot.cb=TRUE)

# Wavelet coherence; nrands should be large (>= 1000)
wtc.t1t2 <- wtc(t1, t2, nrands = 9)

# Plot wavelet coherence
par(oma=c(0, 0, 0, 1), mar=c(5, 4, 4, 5) + 0.1)
plot(wtc.t1t2, plot.cb=TRUE, las=1)

## Plot wavelet coherence and phase difference (arrows)
par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.1)
plot(wtc.t1t2, plot.cb = TRUE, plot.phase = TRUE)

t1 <- cbind(1:100, sin(seq(0, 10 * 2 * pi, length.out = 100)))
t2 <- cbind(1:100, sin(seq(0, 10 * 2 * pi, length.out = 100) + 0.1 * pi))



# Compute dissimilarity
wdist(wt.t1$wave, wt.t2$wave)
