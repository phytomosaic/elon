###################################################################
# IPMs per species in ELON-like quadrats (Keogh):
#   only asexual vegetative GROWTH (no fecundity involved...)
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 01 Dec 2019
#  CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)
#    see: https://cmerow.github.io/RDataScience/21_Intro_IPMs.html
###################################################################

### preamble
rm(list=ls())
require(ecole)
require(viridis)
cv <- viridis::viridis(99)
setwd('~/_prj/7_elon/data/data_census/keogh/')
d  <- read.csv('./label_genets/d_clean.csv', stringsAsFactors=F)
d  <- d[,!colnames(d) %in% c('X.1','X')]
d$species <- ecole::clean_text(as.character(d$species))
names(d)[names(d) == 'objectid'] <- 'genet'
names(d)[names(d) == 'year']     <- 'yr'
d <- d[,c('genet','quad','yr','species','area','tmax','tmin','pcum')]
d <- d[order(d$quad, d$yr),]
a <- aggregate(area ~ genet+yr+quad, data=d, FUN=sum) # SUM duplicates
d <- cbind(a, d[match(interaction(a$genet,a$quad,a$yr, sep='_'),
                      interaction(d$genet,d$quad,d$yr, sep='_')),
                !names(d) %in% c('genet','quad','yr','area')]) ; rm(a)
row.names(d) <- NULL

### transformations
d$area     <- log10(d$area)     # units = log10 m^2
d$areanext <- log10(d$areanext) # units = 1og10 m^2

### get NEXT area value (look-ahead 1) for each quad+genet combn
d <- d[order(d$quad, d$genet, d$yr),]  # order critically important!
d$areanext <- NA
`gnext` <- function(f) {
     nr <- dim(f)[1]
     if (nr == 1) { f$areanext <- NA }
     else if (nr >= 2) { f$areanext <- c(f$area[2:nr], NA) }
     return(f)
}
d <- do.call(rbind, lapply(split(d,list(d$quad,d$genet)), FUN=gnext))
d <- d[order(d$quad, d$genet, d$yr),]         # order again
d <- d[!(is.na(d$area) & is.na(d$areanext)),] # omit NA rows
row.names(d) <- NULL
d_full <- d
rm(d)

###    SELECTION    ##################################################
p_sp <- c('Bouteloua_gracilis')
p_yr <- c('32','33')
######################################################################


### parse 1 species, 1 time-step:
d <- d_full[d_full$species %in% p_sp,] # keep only 1 species
d <- d[d$yr %in% p_yr,]                # keep only 2 observation years

### enforce proper structure: areanext = always Y2, area = always Y1
d$areanext[d$yr == p_yr[2]] <- d$area[d$yr == p_yr[2]]
d$area[d$yr == p_yr[2]]     <- NA
d$surv <- ifelse(is.na(d$areanext),0,1) # survived if Y1 reappeared Y2

######################################################################
# A. data exploration plotting
set_par(4)
plot(d$area,jitter(d$surv),main='Survival')
plot(d$area,d$areanext,main='Growth/Shrinkage/Stasis')
abline(0,1)
# plot(d$size,d$fec.seed,main='Seeds')
hist(d$areanext[d$areanext < 0.02 & !is.na(d$area)], breaks=22, col='grey',
     main='Y2 size distribution, Y1 genets')
hist(d$areanext[d$areanext < 0.02 & is.na(d$area)], breaks=22, col='grey',
     main='Y2 size distribution, Y2 recruits')


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
surv_reg <- glm(surv~area,data=d,family=binomial())
summary(surv_reg)
params$surv_int   <- coefficients(surv_reg)[1]
params$surv_slope <- coefficients(surv_reg)[2]
# 2. growth regression
growth_reg <- lm(areanext~area, data=d)
summary(growth_reg)
params$growth_int   <- coefficients(growth_reg)[1]
params$growth_slope <- coefficients(growth_reg)[2]
params$growth_sd    <- sd(resid(growth_reg))
# # 3. seeds regression
# seed_reg <- glm(fec.seed~size,data=d,family=poisson())
# summary(seed_reg)
# 4. area distribution of recruits
params$recruitareamean <- mean(d$areanext[!is.na(d$area)], na.rm=T)
params$recruitareasd   <- sd(d$areanext[!is.na(d$area)],   na.rm=T)
# # 5. establishment probability
# params$establishment.prob=sum(is.na(d$size))/sum(d$fec.seed,na.rm=T)
# 6. plot the models over the data - figure 2
set_par(4)
xx <- seq(min(d$areanext,na.rm=T),max(d$areanext,na.rm=T), len=44)
plot(d$area,d$areanext,main='Growth/Shrinkage/Stasis')
lines(xx,predict(growth_reg,data.frame(area=xx)),col='red',lwd=3)
abline(0,1, col='grey70')
plot(d$area,jitter(d$surv),main='Survival')
lines(xx,predict(surv_reg,data.frame(area=xx),type='response'), col='red',lwd=3)
# plot(d$area,d$fec,main='Seeds')
# lines(xx,predict(seed_reg,data.frame(area=xx),type='response'),col='red',lwd=3)
hist(d$areanext[!is.na(d$area)], breaks=22, freq=F, col='grey',
     main='Recruit area')
lines(xx, dnorm(xx, params$recruitareamean, params$recruitareasd),
      col='red', lwd=3)

######################################################################
# C. build vital rate functions
# 1. probability of surviving
s_x <- function(x,params) {
     u <- exp(params$surv_int+params$surv_slope*x)
     return(u / (1 + u))
}
# 2. growth function
g_yx <- function(xp,x,params) {
     dnorm(xp, mean = params$growth_int + params$growth_slope * x,
           sd = params$growth_sd)
}
# # 3. reproduction function
# f_yx <- function(xp,x,params) {
#      params$establ_prob *
#           dnorm(xp,mean=params$recruitareamean,
#                 sd=params$recruitareasd)*
#           exp(params$seed_int+params$seed_slope*x)
# }

######################################################################
# D. make kernels
minsz <- 0.9 * min(c(d$area,d$areanext),na.rm=T) # expand ranges
maxsz <- 1.1 * max(c(d$area,d$areanext),na.rm=T)
n <- 100                  # number of cells in the discretized kernel
b <- minsz+c(0:n)*(maxsz-minsz)/n    # boundary points b = cell edges
y <- 0.5*(b[1:n]+b[2:(n+1)])         # mesh points y = cell midpoints
h <- y[2]-y[1]                       # step area h = cell width
G <- h*outer(y,y,g_yx,params=params) # growth kernel
S <- s_x(y,params=params) 	       # survival kernel
P <- G 						  # initialize growth/survival
for(i in 1:n) P[,i] <- G[,i]*S[i]    # growth/survival kernel
# FF <- h*outer(y,y,f_yx,params=params) 	# reproduction kernel
K  <- P # + FF 	                 # complete FULL kernel

######################################################################
# E. basic analyses
(lam    <- Re(eigen(K)$values[1]))     # lambda, popn growth rate
w_eigen <- Re(eigen(K)$vectors[,1])    # right eigenvector
v_eigen <- Re(eigen(t(K))$vectors[,1]) # left eigenvector
stable_dist<- w_eigen/sum(w_eigen)     # stable size distribution
repro_val  <- v_eigen/v_eigen[1]       # reproductive value
v_dot_w <- sum(stable_dist * repro_val) * h
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
# image(sv$stage_agesurv, col=cv)

### collect results
res <- list(lam = lam,
            stable_dist = stable_dist,
            repro_val = repro_val,
            sens = sens,
            elas = elas,
            mle = mle,
            vle = vle,
            psg = psg,
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
image(y,y,t(G), main='Growth kernel', col=cv, xlab=xl, ylab=yl)
contour(y,y,t(G), add=T, drawlabels=T) ; abline(0,1,lwd=3)
image(y,y,t(P), main='Survival/growth kernel', col=cv, xlab=xl,ylab=yl)
contour(y,y,t(P), add=T, drawlabels=T) ; abline(0,1,lwd=3)
# image(y,y,t(FF),main='fecundity kernel', col=cv, xlab=xl,ylab=yl)
# contour(y,y,t(FF), add=T, drawlabels=T) ; abline(0,1,lwd=3)
# plot(0, type='n', axes=F, main='Fecundity kernel', xlab='',ylab='')
image(y,y,t(K)^0.3,main='Full kernel', col=cv, xlab=xl, ylab=yl)
contour(y,y,t(K), add=T, drawlabels=T) ; abline(0,1,lwd=3)
image(y,y,t(elas), xlab=xl, ylab=yl, main='Elasticity', col=cv)
image(y,y,t(sens), xlab=xl, ylab='', main='Sensitivity', col=cv)
plot(y, stable_dist, type='l', xlab='area', ylab='',
     main='Stable area distribution')
plot(y, repro_val, xlab='area', type='l', main='Reproductive values')
plot(y,mle,xlab='area',type='l',main='Life expectancy\n+/- 1 SE')
lines(y, mle + (sqrt(vle)/sqrt(NROW(K))), col=2)
lines(y, mle - (sqrt(vle)/sqrt(NROW(K))), col=2)
# dev.off()

####   END   ####