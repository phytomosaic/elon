###################################################################
# DDMs per species in ELON-like quadrats (Keogh):
#   only asexual vegetative GROWTH (no fecundity involved...)
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 04 Dec 2019
#  CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)
#    see: https://cmerow.github.io/RDataScience/21_Intro_IPMs.html
###################################################################

### preamble
rm(list=ls())
require(ecole)
require(viridis)
require(raster)
cc <- viridis::viridis(99)
setwd('~/_prj/7_elon/elon/_data_raw/data_census/keogh/')
d  <- read.csv('./d_CENSUS_TIDY.csv', row.names=1, stringsAsFactors=F)

### EXPERIMENTAL: emulate measurement error in climate data
set.seed(88)
d$mat  <- d$mat  + rnorm(NROW(d), 0, 0.1)
set.seed(2112)
d$pcum <- d$pcum + rnorm(NROW(d), 0, 5)

### get the env raster (3 variables, average for 1970-2000)
# r <- raster::getData(name='worldclim',var='bio',res=0.5,lon=-105.881,lat=46.383)
# r <- dropLayer(r,c(2:11,13:19)) # only need 1=MAT and 12=MAP
# gain(r[[1]]) <- 0.1 # temperature stored as integer, need up-convert
# (r <- aggregate(r, fact = 20)) # coarsen
# writeRaster(r, filename='env_raster.tif', format='GTiff',
#             options=c('COMPRESS=LZW'), overwrite=TRUE)
r <- stack('env_raster.tif')
names(r) <- c('mat','pcum')
r # 32400 cells for prediction!
plot(r, col=cc)

###    YEAR AND SPECIES SELECTION    #################################
d_full <- d
p_sp   <- c('Bouteloua_gracilis')
d      <- d_full[d_full$species %in% p_sp,]
######################################################################


######################################################################
# A. data exploration plotting
set_par(4)
plot(d$area,jitter(d$surv),main='Survival', col='#00000010')
plot(d$area,d$areanext,main='Growth/Shrinkage/Stasis',col='#00000010')
abline(0,1)
# plot(d$size,d$fec.seed,main='Seeds')
hist(d$areanext[!is.na(d$area)], breaks=22,
     col='grey', main='Y2 size distribution, Y1 genets')
hist(d$areanext[ is.na(d$area)], breaks=22,
     col='grey', main='Y2 size distribution, Y2 recruits')


######################################################################
# B. regressions for vital rate functions
params <- data.frame(
     surv_int=NA,
     surv_slope=NA,
     growth_int=NA,
     growth_slope=NA,
     growth_sd=NA,
     # seed_int=NA,
     # seed_slope=NA,
     recruitareamean=NA,
     recruitareasd=NA
     # establ_prob=NA
)
# 1. survival regression
# surv_reg <- glm(surv~area,data=d,family=binomial(link='logit'))
fmla_s   <- surv ~ area + mat + pcum
surv_reg <- glm(fmla_s, data=d, family=binomial(link='logit'))
summary(surv_reg)
params$surv_int        <- coefficients(surv_reg)[1]
params$surv_slope      <- coefficients(surv_reg)[2]
params$surv_slope_mat  <- coefficients(surv_reg)[3]
params$surv_slope_pcum <- coefficients(surv_reg)[4]
# 2. growth regression
fmla_g     <- areanext ~ area + mat + pcum
growth_reg <- lm(fmla_g, data=d)
summary(growth_reg)
params$growth_int        <- coefficients(growth_reg)[1]
params$growth_slope      <- coefficients(growth_reg)[2]
params$growth_slope_mat  <- coefficients(growth_reg)[3]
params$growth_slope_pcum <- coefficients(growth_reg)[4]
params$growth_sd         <- sd(resid(growth_reg))

# # 2. GAM ALTERNATIVE: growth regression
# require(mgcv)
# fmla_g     <- areanext ~ s(area,m=2) + s(mat,m=3) + s(pcum,m=3)
# growth_reg <- mgcv::gam(fmla_g, data=d)
# summary(growth_reg)
# params$growth_sd <- sqrt(sum(growth_reg$residuals^2) /
#                            growth_reg$df.residual)
# g_yx <- function(y, x, mat, pcum, params) {
#   mu <- predict(growth_reg,
#                 newdata=data.frame(area=x,
#                                    mat=mean(d$mat),
#                                    pcum=mean(d$pcum)),
#                 type='response')
#   return(dnorm(y, mean=mu, sd=params$growth_sd))
# }

# # 3. seeds regression
# seed_reg <- glm(fec.seed~size,data=d,family=poisson())
# summary(seed_reg)
# 4. area distribution of recruits
params$recruitareamean <- mean(d$areanext[!is.na(d$area)], na.rm=T)
params$recruitareasd   <-   sd(d$areanext[!is.na(d$area)], na.rm=T)
# # 5. establishment probability
# params$establishment.prob=sum(is.na(d$size))/sum(d$fec.seed,na.rm=T)


######################################################################
# C. build vital rate functions
# 1. probability of surviving
s_x <- function(x, mat, pcum, params) {
  u <- exp(params$surv_int +
             params$surv_slope * x +
             params$surv_slope_mat * mat +
             params$surv_slope_pcum * pcum)
  return(u / (1 + u))
}
# 2. growth function
g_yx <- function(y, x, mat, pcum, params) {
  dnorm(y,
        mean = params$growth_int +
          params$growth_slope * x +
          params$growth_slope_mat * mat +
          params$growth_slope_pcum * pcum,
        sd = params$growth_sd)
}


# 6. plot the models over the data - figure 2
set_par(3)
xx <- seq(min(d$areanext,na.rm=T),max(d$areanext,na.rm=T), len=44)
plot(d$area,d$areanext,main='Growth/Shrinkage/Stasis', col='#00000010')
lines(xx,predict(growth_reg, data.frame(area=xx,
                                        mat=mean(d$mat),
                                        pcum=mean(d$pcum))),col=2,lwd=3)
abline(0,1, col='grey70')
plot(d$area, jitter(d$surv), main='Survival', col='#00000010')
lines(xx,predict(surv_reg,
                 data.frame(area=xx,
                            mat=mean(d$mat),
                            pcum=mean(d$pcum)),
                 type='response'), col=2, lwd=3)
# plot(d$area,d$fec,main='Seeds')
# lines(xx,predict(seed_reg,data.frame(area=xx),type='response'),col=2,lwd=3)
hist(d$areanext[!is.na(d$area)], breaks=22, freq=F, col='grey',
     main='Recruit area')
lines(xx, dnorm(xx, params$recruitareamean, params$recruitareasd),
      col=2, lwd=3)




### map predicted vital rates
`vitalmap` <- function(reg, d, r, zlims=c(0,4), areaq=0.01,
                       isgrowth=FALSE, ...){
     newd <- data.frame(area = quantile(d$area, areaq, na.rm=T),
                        na.omit(values(r)),
                        row.names = NULL)
     preds <- predict(reg, newd)
     if (isgrowth) {
          preds <- preds - mean(newd$area) # to growth *increment*
          lab <- 'Individual growth'
     } else {
          lab <- 'Survival'
     }
     preds[preds < zlims[1]] <- zlims[1]
     preds[preds > zlims[2]] <- zlims[2]
     rpred  <- r[[1]]  # copy
     nona   <- complete.cases(values(rpred))
     values(rpred)[nona] <- preds
     plot(rpred, col=cc, xaxt='n', yaxt='n', main=lab)
}
set_par(6)
vitalmap(growth_reg, d, r, zlims=c(0.5,1.1), areaq=0.05, isgrowth=T)
vitalmap(growth_reg, d, r, zlims=c(0.5,1.1), areaq=0.25, isgrowth=T)
vitalmap(growth_reg, d, r, zlims=c(0,  1.1), areaq=0.50, isgrowth=T)
vitalmap(surv_reg, d, r, zlims=c(-2,2), areaq=0.05)
vitalmap(surv_reg, d, r, zlims=c(-2,2), areaq=0.25)
vitalmap(surv_reg, d, r, zlims=c(-2,2), areaq=0.50)


# ######################################################################
# # C. build vital rate functions
# # 1. probability of surviving
# s_x <- function(x,params) {
#   u <- exp(params$surv_int+params$surv_slope*x)
#   return(u / (1 + u))
# }
# # 2. growth function
# g_yx <- function(xp,x,params) {
#   dnorm(xp, mean = params$growth_int + params$growth_slope * x,
#         sd = params$growth_sd)
# }
# # # 3. reproduction function
# # f_yx <- function(xp,x,params) {
# #      params$establ_prob *
# #           dnorm(xp,mean=params$recruitareamean,
# #                 sd=params$recruitareasd)*
# #           exp(params$seed_int+params$seed_slope*x)
# # }


# ######################################################################
######################################################################
# D. make kernels ITERATIVELY for each cell of a raster

### kernel matrix setup
minsz <- 0.9 * min(c(d$area,d$areanext),na.rm=T) # expand ranges
maxsz <- 1.1 * max(c(d$area,d$areanext),na.rm=T)
n <- 100                  # number of cells in the discretized kernel
b <- minsz+c(0:n)*(maxsz-minsz)/n    # boundary points b = cell edges
y <- 0.5*(b[1:n]+b[2:(n+1)])         # mesh points y = cell midpoints
h <- y[2]-y[1]                       # step area h = cell width

# ### for a SINGLE kernel:
mat  <- mean(d$mat)
pcum <- mean(d$pcum)
G <- h*outer(y,y,g_yx,params=params, mat, pcum) # growth kernel
# # growth correction??? rerouting growth to sizes outside the
# # allowed range to the extreme sizes (avoiding eviction):
# for(i in 1:(n/2))   G[1,i] <- G[1,i] + 1 - sum(G[,i])
# for(i in (n/2+1):n) G[n,i] <- G[n,i] + 1 - sum(G[,i])
S <- s_x(y,params=params, mat, pcum) 	     # survival kernel
P <- G 				     # initialize growth/survival
for(i in 1:n) P[,i] <- G[,i]*S[i]    # growth/survival kernel
# FF <- h*outer(y,y,f_yx,params=params) # reproduction kernel
K  <- P # + FF 	                     # complete FULL kernel
image(K, col=cc)


### to ITERATE for each raster cell (or 'xval' can be data.frame)
`calc_lam` <- function(ii) {
  if (ii%%100 == 0) { cat(ii, '... ') } # progress
  mat  <- xval[ii,'mat']
  pcum <- xval[ii,'pcum']
  G <- h*outer(y,y,g_yx,params=params, mat, pcum) # growth kernel
  S <- s_x(y,params=params, mat, pcum) 	     # survival kernel
  P <- G 				             # initialize
  for(i in 1:n) P[,i] <- G[,i]*S[i]   # growth/survival kernel
  # FF <- h*outer(y,y,f_yx,params=params) # reproduction kernel
  K   <- P # + FF 	                # complete FULL kernel
  lam <- tryCatch({Re(eigen(K)$values[1])},# lambda, pop gro rate
                  error=function(cond) {
                    message(paste("error-cell number",ii))
                    return(NA)})
  return(lam)
}
xval   <- values(r)
head(xval)
isna   <- is.na(xval[,'mat']) | is.na(xval[,'pcum'])
xval   <- xval[!isna,]
sum(isna) / ncell(r) # what proportion are NA?
lams   <- sapply(1:NROW(xval), calc_lam) # ! ! ! ! TIMEWARN ! ! ! !
rlam   <- r[[1]]
rlam[] <- NA
values(rlam)[!isna] <- lams
### the final LAMBDA MAP ! ! !
set_par(1)
plot(rlam, col=cc)
# writeRaster(rlam, filename='lambda_raster.tif', format='GTiff',
#             options=c('COMPRESS=LZW'), overwrite=TRUE)



######################################################################
# E. basic analyses

### 1.  Use popbio package
require(popbio)
eig <- popbio::eigen.analysis(K)
fm  <- popbio::fundamental.matrix(K)
fm$N[fm$N > 0.8] <- NA
image(fm$N, col=cc) # mean of the time spent in each stage class


### 2.  can manually calculate
(lam       <- Re(eigen(K)$values[1]))     # lambda, popn growth rate
w_eigen    <- Re(eigen(K)$vectors[,1])    # right eigenvector
v_eigen    <- Re(eigen(t(K))$vectors[,1]) # left eigenvector
stable_dist<- w_eigen/sum(w_eigen)        # stable size distribution
repro_val  <- v_eigen/v_eigen[1]          # reproductive value
v_dot_w    <- sum(stable_dist * repro_val) * h
# SENSITIVITY = change in lam given a small change in a matrix element
sens  <- outer(repro_val, stable_dist) / v_dot_w
# ELASTICITY = 'proportional' effect, i.e., the effect that a change
elas  <- matrix(as.vector(sens) * as.vector(K) / lam, nrow=n)
# MLE = mean and variance of life expectancy for every starting size
# mle   <- colSums(MASS::ginv(diag(NROW(K)) - K))
le <- function(K, method = c('mean','var')){
     tmp <- MASS::ginv(diag(NROW(K))-K)
     if (method == 'mean') {
          return(colSums(tmp))
     }
     if (method == 'var') {
          return(colSums(2*(tmp%*%tmp)-tmp)-colSums(tmp)*colSums(tmp))
     }
}
mle <- le(K, 'mean')
vle <- le(K, 'var')

### collect results
res <- list(lam = lam,
            stable_dist = stable_dist,
            repro_val = repro_val,
            sens = sens,
            elas = elas,
            mle = mle,
            vle = vle,
            # psg = psg,
            y = y,
            G = G,
            P = P,
            K = K)
# return(res)


### plot results
xl <- 'area (t)'
yl <- 'area (t+1)'
# tiff('./fig/fig_00_ipm_kernel.tif', wid=8, hei=6,
#      units='in', compr='lzw+p', bg='transparent', res=500)
par(mfrow=c(2,4), oma=c(0,0,0,0), mar=c(4,4,2,0)+.1, mgp=c(2.0,0.5,0),
    las=1, bty='l')
image(y,y,t(G), main='Growth kernel', col=cc, xlab=xl, ylab=yl)
contour(y,y,t(G), add=T, drawlabels=T) ; abline(0,1,lwd=3)
image(y,y,t(P), main='Survival/growth kernel', col=cc, xlab=xl,ylab=yl)
contour(y,y,t(P), add=T, drawlabels=T) ; abline(0,1,lwd=3)
# image(y,y,t(FF),main='fecundity kernel', col=cc, xlab=xl,ylab=yl)
# contour(y,y,t(FF), add=T, drawlabels=T) ; abline(0,1,lwd=3)
# plot(0, type='n', axes=F, main='Fecundity kernel', xlab='',ylab='')
image(y,y,t(K)^0.3,main='Full kernel', col=cc, xlab=xl, ylab=yl)
contour(y,y,t(K), add=T, drawlabels=T) ; abline(0,1,lwd=3)
image(y,y,t(elas), xlab=xl, ylab=yl, main='Elasticity', col=cc)
image(y,y,t(sens), xlab=xl, ylab='', main='Sensitivity', col=cc)
plot(y, stable_dist, type='l', xlab='area', ylab='',
     main='Stable area distribution')
plot(y, repro_val, xlab='area', type='l', main='Reproductive values')
plot(y,mle,xlab='area',type='l',main='Life expectancy\n+/- 1 SE')
lines(y, mle + (sqrt(vle)/sqrt(NROW(K))), col=2)
lines(y, mle - (sqrt(vle)/sqrt(NROW(K))), col=2)
# dev.off()

### perspective plot
set_par(1)
persp(y, y, t(K), col=surfcol(t(K), ngrid=length(y), pal=cc),
      border=NA, xlab=xl, ylab=yl, main='Full kernel',
      expand=0.6, phi=30, theta=-30)






plot(r, col=cc)











####    UNUSED FUNCTIONS    ####

# # PASSAGE TIME: how many years to get to a specified size?
# passagetime <- function(K, chosenSize, meshpoints=y) {
#      ### modified from IPMpack::survivorship
#      loc <- which(abs(chosenSize - meshpoints) ==
#                        min(abs(chosenSize - meshpoints)),
#                   arr.ind=TRUE)[1]
#      nr <- NROW(K)
#      Tprime        <- K
#      Tprime[,loc]  <- 0
#      Mprime        <- 1 - colSums(K)
#      Mprime[loc]   <- 0
#      Mprime        <- rbind(Mprime, rep(0,nr))
#      Mprime[2,loc] <- 1
#      Bprime        <- Mprime %*% MASS::ginv(diag(nr) - Tprime)
#      Bprime[2,][Bprime[2,]==0] <- 1
#      diagBprime    <- diag(Bprime[2,])
#      Tc   <- diagBprime %*% Tprime %*% MASS::ginv(diagBprime)
#      eta1 <- MASS::ginv(diag(nr) - Tc)
#      time_to_absorb <- colSums(eta1)
#      time_to_absorb[loc:length(time_to_absorb)] <- 0
#      return(time_to_absorb)
# }
# psg <- passagetime(K, max(d$area, na.rm=T))
# set_par(1)
# plot(NA, ylab = "Passage time", xlab = "Continuous size stage",
#      ylim=c(0,max(psg)), xlim=range(y))
# trgsizes <- quantile(d$area, probs=c(seq(0.05,0.95,len=10),0.9999), na.rm=T)
# cc       <- viridis::inferno(length(trgsizes), end=0.8)
# for (i in seq_along(trgsizes)) {
#      lines(y, passagetime(K, trgsizes[i]), col=cc[i])
#      abline(v=trgsizes[i], col=cc[i])
# }
# rm(i, cc, psg, trgsizes, passagetime)

# # SURVIVORSHIP
# survivorship <- function(K, loc, maxage=50) {
#      ### modified from IPMpack::survivorship
#      nr <- NROW(K)
#      n  <- nr
#      A1 <- tmp <-  K
#      stage_agesurv <- matrix(NA,n,maxage)
#      surv_curv     <- rep(NA,maxage)
#      popvec        <- matrix(0,n,1)
#      popvec[floor(loc),1] <- 1
#      for (a in 1:maxage) {
#           surv_curv[a] <- sum(A1[,loc])
#           stage_agesurv[c(1:n),a] <- A1[,] %*% popvec
#           A1 <- A1 %*% tmp
#      }
#      mortality <- -log(surv_curv[2:length(surv_curv)] /
#                             surv_curv[1:(length(surv_curv)-1)])
#      return(list(surv_curv = surv_curv,   # survivorship up to max age
#                  stage_agesurv = stage_agesurv, # sv by stage
#                  mortality = mortality))    # mortality per age
# }
# sv <- survivorship(K, loc=1, maxage=20)
# set_par(3)
# plot(sv$surv_curv)
# plot(sv$mortality)
# image(sv$stage_agesurv, col=cc)


####   END   ####