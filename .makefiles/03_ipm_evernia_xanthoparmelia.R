###################################################################
#
# IPMs for survival and growth (no fecundity) of:
#
#   1) Evernia mesomorpha at SPRUCE experiment in N. Minnesota
#   2) Xanthoparmelia at North Cemetery, Petersham MA (Pringle)
#
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 09 Dec 2019
#  CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)
#
###################################################################


### preamble
rm(list=ls())
require(ecole)
require(viridis)
require(raster)
require(mgcv)
cc <- viridis::viridis(99)


###################################################################
###################################################################
###################################################################
###   1) Evernia mesomorpha at SPRUCE experiment in N. Minnesota
###################################################################
###################################################################
###################################################################


### load and format data
f <- 'C:/Users/Rob/Documents/_prj/5_spruce/data/'
f <- paste0(f, '2019_data/mass_transplants_26Aug2019.csv')
x <- read.csv(f, header=T, stringsAsFactors=F)
x <- x[,c('encl','temptrmt','co2trmt','transplant','corrected2013',
          'corrected2014','corrected2015new','corrected2016',
          'corrected2017', 'corrected2019')]
names(x)[names(x) == 'transplant'] <- 'genet'
names(x)[names(x) == 'corrected2015new'] <- 'corrected2015'

### filter trouble cases
tmp <- x[,grep('corrected', names(x))]
tmp[tmp <= 0] <- NA
x[,grep('corrected', names(x))] <- tmp  ;  rm(tmp)

### get sizenext from next year's size:
w <- x[,!names(x) %in% c('encl','temptrmt','co2trmt')]
# if one year was skipped, assign MEAN of adjoining years
`infill` <- function(x) { # returns vector with single NAs smoothed
  s       <- 2:(length(x)-1) # ignore first/last column
  isna    <- s[which(is.na(x[s]))]
  x[isna] <- (x[isna-1] + x[isna+1]) / 2 # mean of neighbors
  x   # ! ! ! ! TIMEWARN ! ! ! ! need to vectorize
}
w <- data.frame(genet=w[,1], t(apply(w[,-1], 1, infill)))
nr <- NROW(w)       # number of unique genets
s  <- 2:(NCOL(w)-1) # column sequence for years
length(s)           # number of census periods to consider
m <- data.frame(    # collector data.frame
  matrix(NA, nrow=nr*length(s), ncol=4,
         dimnames=list(NULL,c('genet','yr','size','sizenext'))),
  stringsAsFactors=F)
for (j in s) {
  yr    <- gsub('corrected','',dimnames(w)[[2]][j])
  i     <- (1:nr) + (nr*(j-2))
  m[i,] <- data.frame(w$genet, yr, w[,j], w[,j+1],
                      stringsAsFactors=F)
}
m <- m[!(is.na(m$size) & is.na(m$sizenext)),] # omit if both yrs empty
d <- m
row.names(d) <- NULL

### match enclosure id
d$encl      <- x[match(d$genet, x$genet),'encl']
d$temptrmt  <- x[match(d$genet, x$genet),'temptrmt']
rm(nr,w,s,i,j,m,x,f,yr,infill)

### match climate actually experienced (cf prior timeseries analysis)
x <- data.frame(matrix(c(
  4,6.5673562, -23.8178038,
  5,0.0000000,   0.0000000,
  6,2.2187861,  -3.2270304,
  8,8.3749160, -28.4479678,
  10,10.2387452, -36.2656537,
  11,4.3896525, -16.4994324,
  13,6.6832634, -23.5557445,
  14,0.1693923,   0.4703195,
  16,8.4354980, -29.7021902,
  17,10.5573296, -36.3673755,
  19,1.9066682,  -3.6680314,
  20,4.4479378, -16.2528410),
  nrow=12, ncol=3, byrow=T,dimnames=list(NULL,c('encl','temp','rh'))))
d <- cbind(d, x[match(d$encl, x$encl),-1])  ; rm(x)

### assign survival
d$surv <- ifelse(is.na(d$sizenext),0,1) # survived if Y1 reappeared Y2

### explore
table(d$yr) # every year 2005-2011
table(d$genet) # n observations per individual
table(d$temptrmt)
table(d$temp)
table(d$encl)
# xtabs(~ genet+yr, data=d)

### data exploration plots
set_par(3)
plot(d$size,d$sizenext,main='Growth/Shrinkage/Stasis',
     col=colvec(d$temp), pch=16, cex=0.9) ; abline(0,1)
plot(d$size,jitter(d$surv),main='Survival',
     col=colvec(d$temp), pch=16, cex=0.9)
hist(d$size, breaks=33, col='grey', main='Y1 size distribution')

### SUBSET only FIRST treatment year? ? ?
d <- d[d$yr == 2015,]

######################################################################
# B. regressions for vital rate functions
params <- data.frame(
  surv_int=NA,
  surv_slope=NA,
  recruitareamean=NA,
  recruitareasd=NA
)
# 1. survival regression
fmla_s   <- surv ~ size
surv_reg <- glm(fmla_s, data=d, family=binomial(link='logit'))
summary(surv_reg)
params$surv_int        <- coefficients(surv_reg)[1]
params$surv_slope      <- coefficients(surv_reg)[2]

# 2. growth regression - GAM with constant variance
growth_reg <- mgcv::gam(sizenext ~ te(size,temp, bs='tp', m=2),
                        data=d, family='gaussian')
summary(growth_reg)

# 3. size distribution of recruits
params$recruitareamean <- mean(d$sizenext[!is.na(d$size)], na.rm=T)
params$recruitareasd   <-   sd(d$sizenext[!is.na(d$size)], na.rm=T)

### plot models over the data
set_par(3)
xx <- seq(min(d$sizenext,na.rm=T),max(d$sizenext,na.rm=T), len=44)
plot(d$size,d$sizenext,main='Growth/Shrinkage/Stasis',col='#000000')
u_temp <- sort(unique(d$temp))
cu     <- colvec(u_temp)
for (u in 1:length(u_temp)) {
  lines(xx, predict(growth_reg, data.frame(size=xx,temp=u_temp[u])),
        col=cu[u], lwd=3)
}
abline(0,1, col='grey70')
plot(d$size, jitter(d$surv), main='Survival', col='#000000')
lines(xx,predict(surv_reg,data.frame(size=xx),type='response'),
      col=2,lwd=3)
hist(d$sizenext[!is.na(d$size)], breaks=22, freq=F, col='grey',
     main='Recruit size')
lines(xx, dnorm(xx, params$recruitareamean, params$recruitareasd),
      col=2, lwd=3)


######################################################################
# C. build vital rate functions
# 1. probability of surviving
s_x <- function(x, params) {
  u <- exp(params$surv_int +
             params$surv_slope * x)
  return(u / (1 + u))
}
# 2. growth function - GAM model
g_yx <- function(y, x, reg, temp) {
  p   <- mgcv:::predict.gam(reg, data.frame(size = x,
                                            temp = temp), se.fit=T)
  mu  <- p$fit
  sig <- p$se.fit * 44
  dnorm(y, mean = mu, sd= abs(sig))
}

######################################################################
# D.  kernel matrix setup
minsz <- 0.9 * min(c(d$size,d$sizenext),na.rm=T) # expand ranges
maxsz <- 1.1 * max(c(d$size,d$sizenext),na.rm=T)
n <- 100                  # number of cells in the discretized kernel
b <- minsz+c(0:n)*(maxsz-minsz)/n    # boundary points b = cell edges
y <- 0.5*(b[1:n]+b[2:(n+1)])         # mesh points y = cell midpoints
h <- y[2]-y[1]                       # step size h = cell width


######################################################################
# E. basic population analyses for a SINGLE kernel:
temp <- mean(d$temp) # some value of interest
G <- h*outer(y,y,g_yx,reg=growth_reg, temp) # GAM growth kernel
S <- s_x(y,params=params) 	         # survival kernel
P <- G 				                       # initialize growth/survival
for(i in 1:n) P[,i] <- G[,i]*S[i]    # growth/survival kernel
K  <- P # + FF 	                     # complete FULL kernel
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
            mle  = mle,
            vle  = vle,
            y = y,
            G = G,
            P = P,
            K = K)

### ITERATE kernel for each TEMPERATURE value
`plot_k` <- function(ii) {
  temp <- xval[ii,'temp']
  G <- h*outer(y,y,g_yx,reg=growth_reg, temp) # GAM growth kernel
  S <- s_x(y,params=params) 	         # survival kernel
  P <- G 				                       # initialize growth/survival
  for(i in 1:n) P[,i] <- G[,i]*S[i]    # growth/survival kernel
  K  <- P 	                           # complete FULL kernel
  lam <- tryCatch({Re(eigen(K)$values[1])},# lambda, pop gro rate
                  error=function(cond) {
                    message(paste("error-cell number",ii))
                    return(NA)})
  image(y,y,t(K),main=paste0('T = ',round(temp, 2)), col=cc,
        xlab=xl, ylab=yl)
  contour(y,y,t(K), add=T, drawlabels=F) ; abline(0,1,lwd=2)
  ij <- which(t(K) == max(t(K)), arr.ind = TRUE)
  points(y[ij[1]], y[ij[2]], col=2, pch=16, cex=1)
  add_text(0.1,0.90,labels=round(lam,3),bold=T,col='white')
}
xval   <- data.frame(temp=u_temp)
set_par(NROW(xval))
sapply(1:NROW(xval), plot_k)

### TODO: consider EACH yearly census, and yearly climate







###################################################################
###################################################################
###################################################################
###   2) Xanthoparmelia at North Cemetery, Petersham MA (Pringle)
###################################################################
###################################################################
###################################################################




rm(list=ls())

### load and format data
u <- paste0(
     'https://portal.edirepository.org/nis/dataviewer?packageid=knb-',
     'lter-hfr.217.4&entityid=1b5854c071148a026825a11f5891e234')
x <- read.csv(url(u), header=T, stringsAsFactors=F)
head(x)
x$grave <- NULL # unneeded column
x$individual <- as.numeric(as.factor(x$individual))
x$alive <- as.factor(x$alive)
x$merge <- as.factor(x$merge)
xtabs( ~ individual + year, x)
table(x$year) # every year 2005-2011
table(x$individual) # seven observations per individual

### plotting
ecole::set_par(2)
plot(sort(log10(x$size + 1))) # rank-size curve
hist(x$size, breaks=44)       # size = squareCentimeters = Area

### convert Area to Radius
x$r <- sqrt(x$size / pi)
x$r <- x$r * 10 # convert from cm to mm

### how many births? deaths? mergers?

### filter out thalli that merged/fragmented (14 of 108 genets)?
rms <- unique(x$individual[x$merge %in% c('1','3')])
x <- x[!x$individual %in% rms,]  ;  rm(rms)

### get size and sizeNEXT in proper order for IPMs
d <- x[,!names(x) %in% c('size')]
head(d)
names(d)[names(d) == 'individual'] <- 'genet'
names(d)[names(d) == 'year']       <- 'yr'
names(d)[names(d) == 'r']          <- 'size' # use RADIUS as size
d <- d[order(d$genet, d$yr),]
xtabs(~ genet+yr, data=d)
d$size <- log10(d$size + 1)                 # units = log10 mm
row.names(d) <- NULL

### get sizenext from next year's size:
# reshape wide, where each column = size in one year
w <- reshape(d[,c('genet','yr','size')],v.names='size',idvar='genet',
             timevar='yr', direction='wide', sep='_')
head(w)
ecole::plot_heatmap(w[,-1])

### then iteratively append year-pairs
nr <- NROW(w)       # number of unique genets
s  <- 2:(NCOL(w)-1) # column sequence for years
length(s)           # number of census periods to consider
m <- data.frame(    # collector data.frame
     matrix(NA, nrow=nr*length(s), ncol=4,
            dimnames=list(NULL,c('genet','yr','size','sizenext'))),
     stringsAsFactors=F)
for (j in s) {
     yr    <- gsub('^[^_]*_','',dimnames(w)[[2]][j])
     i     <- (1:nr) + (nr*(j-2))
     m[i,] <- data.frame(w$genet, yr, w[,j], w[,j+1],
                         stringsAsFactors=F)
}
rm(nr,w,s,i,j)
m <- m[!(is.na(m$size) & is.na(m$sizenext)),] # omit if both yrs empty
d <- m
rm(m,x)

### assign survival
d$surv <- ifelse(is.na(d$sizenext),0,1) # survived if Y1 reappeared Y2


######################################################################
# A. data exploration plotting
set_par(4)
plot(d$size,d$sizenext,main='Growth/Shrinkage/Stasis',col='#000000')
plot(d$size,jitter(d$surv),main='Survival', col='#000000')
abline(0,1)
hist(d$sizenext[!is.na(d$size)], breaks=33,
     col='grey', main='Y2 size distribution, Y1 genets')
hist(d$sizenext[ is.na(d$size)], breaks=33,
     col='grey', main='Y2 size distribution, Y2 recruits')

######################################################################
# B. regressions for vital rate functions
params <- data.frame(
     surv_int=NA,
     surv_slope=NA,
     growth_int=NA,
     growth_slope=NA,
     growth_sd=NA,
     recruitareamean=NA,
     recruitareasd=NA
)
# 1. survival regression
fmla_s   <- surv ~ size
surv_reg <- glm(fmla_s, data=d, family=binomial(link='logit'))
summary(surv_reg)
params$surv_int        <- coefficients(surv_reg)[1]
params$surv_slope      <- coefficients(surv_reg)[2]

# 2. growth regression - LINEAR
fmla_g     <- sizenext ~ size
growth_reg <- lm(fmla_g, data=d)
summary(growth_reg)
params$growth_int        <- coefficients(growth_reg)[1]
params$growth_slope      <- coefficients(growth_reg)[2]
params$growth_sd         <- sd(resid(growth_reg))

# 4. size distribution of recruits
params$recruitareamean <- mean(d$sizenext[!is.na(d$size)], na.rm=T)
params$recruitareasd   <-   sd(d$sizenext[!is.na(d$size)], na.rm=T)

### plot models over the data
set_par(3)
xx <- seq(min(d$sizenext,na.rm=T),max(d$sizenext,na.rm=T), len=44)
plot(d$size,d$sizenext,main='Growth/Shrinkage/Stasis',col='#000000')
p <- predict(growth_reg, data.frame(size=xx))
lines(xx, p,  col=2, lwd=3)
abline(0,1, col='grey70')
plot(d$size, jitter(d$surv), main='Survival', col='#000000')
lines(xx,predict(surv_reg,data.frame(size=xx),type='response'),
      col=2,lwd=3)
hist(d$sizenext[!is.na(d$size)], breaks=22, freq=F, col='grey',
     main='Recruit size')
lines(xx, dnorm(xx, params$recruitareamean, params$recruitareasd),
      col=2, lwd=3)

######################################################################
# C. build vital rate functions
# 1. probability of surviving
s_x <- function(x, params) {
     u <- exp(params$surv_int +
                   params$surv_slope * x)
     return(u / (1 + u))
}
# 2. growth function - LINEAR model
g_yx <- function(y, x, params) {
  dnorm(y,
        mean = params$growth_int +
          params$growth_slope * x,
        sd = params$growth_sd)
}

######################################################################
# D. kernel matrix setup
minsz <- 0.9 * min(c(d$size,d$sizenext),na.rm=T) # expand ranges
maxsz <- 1.1 * max(c(d$size,d$sizenext),na.rm=T)
n <- 100                  # number of cells in the discretized kernel
b <- minsz+c(0:n)*(maxsz-minsz)/n    # boundary points b = cell edges
y <- 0.5*(b[1:n]+b[2:(n+1)])         # mesh points y = cell midpoints
h <- y[2]-y[1]                       # step size h = cell width


######################################################################
# E. basic population analyses

### for a SINGLE kernel:
G <- h*outer(y,y,g_yx,params=params) # LINEAR growth kernel
S <- s_x(y,params=params) 	     # survival kernel
P <- G 				     # initialize growth/survival
for(i in 1:n) P[,i] <- G[,i]*S[i]    # growth/survival kernel
# FF <- h*outer(y,y,f_yx,params=params) # reproduction kernel
K  <- P # + FF 	                     # complete FULL kernel

### 1.  manually calculate
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
            mle  = mle,
            vle  = vle,
            # psg = psg,
            y = y,
            G = G,
            P = P,
            K = K)
# return(res)

### Kernel plotting
xl <- 'size (t)'
yl <- 'size (t+1)'
# tiff('./fig/fig_00_ipm_kernel.tif', wid=8, hei=6,
#      units='in', compr='lzw+p', bg='transparent', res=500)
par(mfrow=c(2,4), oma=c(0,0,0,0), mar=c(4,4,2,0)+.1, mgp=c(2.0,0.5,0),
    las=1, bty='l')
image(y,y,t(G), main='Growth kernel', col=cc, xlab=xl, ylab=yl)
contour(y,y,t(G), add=T, drawlabels=T) ; abline(0,1,lwd=3)
image(y,y,t(P), main='Survival/growth kernel', col=cc,xlab=xl,ylab=yl)
contour(y,y,t(P), add=T, drawlabels=T) ; abline(0,1,lwd=3)
image(y,y,t(K),main='Full kernel', col=cc, xlab=xl, ylab=yl)
contour(y,y,t(K), add=T, drawlabels=T) ; abline(0,1,lwd=3)
image(y,y,t(elas), xlab=xl, ylab=yl, main='Elasticity', col=cc)
image(y,y,t(sens), xlab=xl, ylab='', main='Sensitivity', col=cc)
plot(y, stable_dist, type='l', xlab='size', ylab='',
     main='Stable size distribution')
plot(y, repro_val, xlab='size', type='l', main='Reproductive values')
plot(y,mle,xlab='size',type='l',main='Life expectancy\n+/- 1 SE')
lines(y, mle + (sqrt(vle)/sqrt(NROW(K))), col=2)
lines(y, mle - (sqrt(vle)/sqrt(NROW(K))), col=2)
# dev.off()

### perspective plot
set_par(1)
persp(y, y, t(K), col=surfcol(t(K), ngrid=length(y), pal=cc),
      border=NA, xlab=xl, ylab=yl, main='Full kernel',
      expand=0.6, phi=30, theta=-30)


####    END    ####