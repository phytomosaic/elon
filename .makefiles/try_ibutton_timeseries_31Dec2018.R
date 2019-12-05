######################################################################
# first ELON-Portland iButton data
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 31 Dec 2018
##  CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

###   preamble
setwd('~/_prj/7_elon/fig/')
rm(list=ls())
pkg <- c('viridis','zoo','plantecophys','biwavelet')
has <- pkg %in% rownames(installed.packages())
if(any(!has))install.packages(pkg[!has])
lapply(pkg, require, character.only=TRUE)
rm(pkg, has)

###   load data, acquired 25 Sept - 26 Dec 2018
cols <- inferno(7, alpha=0.7)
p <- 'C:/Users/Rob/Documents/_prj/7_elon/data/trials_ibutton_26Dec2018'
f  <- '/trials_ibutton_26Dec2018.csv'
d  <- read.csv(file=paste0(p, f), header=T, stringsAsFactors=F)
dt <- as.POSIXct(d[,1],format='%m/%d/%Y %H:%M',tz='America/Los_Angeles')
z  <- zoo(d[,c(2,4)], dt) # pick out a single ibutton T and RH
colnames(z) <- gsub('0','',tolower(colnames(z)))
colnames(z) <- c('t','rh')
rm(dt, p, f, d)

### TODO: quality checks
head(z)
anyNA(z)
any(z[,'rh'] < 0)
z[which(z[,'rh'] < 0),]

### indices for daily aggregates
bys <- as.POSIXct(as.character(index(z)),
                  format='%Y-%m-%d') # by DAILY index
# bys <- as.Date(index(z))             # by DAILY index
# bys <- as.POSIXct(d[,1],format='%m/%d/%Y %H',
#                   tz='America/Los_Angeles') # by hourly index
# bys <- as.yearmon(index(z))          # by monthly index

### function for basic summary statistics (applied DAILY)
`f` <- function(x, prefix = 't', ...){
     x[!is.finite(x)] <- NA
     m <- cbind(
          mean(x, na.rm=TRUE),
          # sd(x, na.rm=TRUE), # almost always in-phase w range
          max(x, na.rm=TRUE),
          min(x, na.rm=TRUE),
          diff(range(x, finite=TRUE))
     )
     m[!is.finite(m)] <- NA
     dimnames(m)[[2]] <-
          paste0(prefix,c('_mean','_max','_min','_rng'))
     out <- zoo(m)
     out
}

### calc DAILY statistics
daily <- cbind(aggregate(z[,'t'], by=bys, FUN=f, prefix='t'),
               aggregate(z[,'rh'], by=bys, FUN=f, prefix='rh'))

### calc daily VPD (units = kPa) ?RHtoVPD
vpd <- plantecophys::RHtoVPD(RH    = daily[,'rh_mean'],
                             TdegC = daily[,'t_mean'],
                             Pa    = 101)
daily <- cbind(daily, vpd=vpd)

### calc daily ET (Hargreaves)
`hargreaves` <- function (z, tmax = 't_max', tmin = 't_min',
                          lat = 45, elev = 45, ...){
     ### per `Evapotranspiration::ET.HargreavesSamani`
     if(!inherits(z,'zoo'))
          stop('`z` must be `zoo` object')
     if (is.null(tmax) | is.null(tmin))
          stop('Required tmax,tmin are missing')
     ii    <- index(z)
     Tmax  <- z[,tmax]
     Tmin  <- z[,tmin]
     J     <- as.POSIXlt(ii)$yday # julian day
     lambda<- 2.45  # latent heat of evaporation = 2.45 MJ kg^-1 @20C
     Gsc   <- 0.082 # solar constant = 0.0820 MJ m^-2 min^-1
     Ta    <- (Tmax + Tmin) / 2
     P     <- 101.3 * ((293 - 0.0065 * elev)/293)^5.26
     delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta + 237.3))) /
          ((Ta + 237.3)^2)
     gamma <- 0.00163 * P / lambda
     d_r2  <- 1 + 0.033 * cos(2 * pi/365 * J)
     delta2<- 0.409 * sin(2 * pi/365 * J - 1.39)
     latrad<- lat * pi / 180 # latitude in radians
     w_s   <- acos(-tan(latrad) * tan(delta2))
     N     <- 24/pi * w_s
     R_a   <- (1440/pi) * d_r2 * Gsc * (
          w_s * sin(latrad) * sin(delta2) + cos(latrad) *
               cos(delta2) * sin(w_s)
     )
     Tdiff <- Tmax - Tmin
     C_HS  <- 0.00185 * Tdiff^2 - 0.0433 * Tdiff + 0.4023
     ET    <- 0.0135 * C_HS * R_a/lambda * Tdiff^0.5 * (Ta + 17.8)
     message('*** Hargreaves-Samani Reference Crop ET ***')
     message('Times: ',min(ii), ' to ', max(ii))
     message('Mean: ', round(mean(ET, na.rm = T),
                             digits = 2))
     message('Max: ', round(max(ET, na.rm = T), digits = 2))
     message('Min: ', round(min(ET, na.rm = T), digits = 2))
     return(zoo(ET, ii))
}
et <- hargreaves(daily, lat=45.464880, elev=30)
daily <- cbind(daily, et=et)
rm(f, vpd, et)
plot(daily)

# ### OPTIONAL: match daily to intra-daily
# tmp <- coredata(daily[match(bys, index(daily)),])
# z <- cbind(z, tmp)
# rm(tmp)
# head(z)
# plot(z)

### PLOT daily
par(bty='L', mar=c(0,3,0,0))
plot(daily, screens=c(1,1,1,3,2,2,2,4,5,6), las=1, mgp=c(2.5,.7,0),
     col=c(1,2,2,1,1,2,2,1,1,1), oma=c(2,0,1,0), xaxs='i', nc=1)

### PLOT subset, pretty
# tiff('C:/Users/Rob/Desktop/fig_00_cli.tif', hei=7, wid=7,
#      unit='in', res=500, compr='lzw+p')
par(bty='L')
tmp <- daily[,!names(daily) %in% c('vpd','t_rng','rh_rng')]
names(tmp)
plot(tmp, main='Daily climate, Portland Oregon USA',
     screens=c(1,1,1,2,2,2,3), las=1, mgp=c(2.5,.7,0),
     col=c(1,2,2,1,2,2,1), oma=c(2,.5,1,.1), mar=c(.1,4,.5,.1),
     ylab=c('T (°C)','RH (%)','ET (mm)'),
     xlab='', nc=1, type='l')
rm(tmp)
# dev.off()

### plot boxes around threshold values
`plot_thresh` <- function(z, yvar, yval, dir='greater',
                          vline=FALSE, pcex=NA, lwd=2, ...){
     if(!inherits(z,'zoo')) stop('`z` must be `zoo` object')
     xx <- index(z)
     switch(as.character(dir),
            `greater` = ii <- match(xx[which(z[,yvar] > yval)], xx),
            `less`    = ii <- match(xx[which(z[,yvar] < yval)], xx),
            `equal`   = ii <- match(xx[which(z[,yvar] ==yval)], xx))
     ii <- ii[!ii %in% c(1,length(xx))] # cant include edges
     `p` <- function(x, y, ...) {
          points(x, y, cex=pcex)
          lines(x, y, lwd=lwd)
          if(vline) abline(v = xx[ii], col = 3)
          rect(xleft   = xx[ii+1],
               ybottom = min(y, na.rm=TRUE)-99,
               xright  = xx[ii-1],
               ytop    = max(y, na.rm=TRUE)+99,
               col = '#00000010', border = NA)
     }
     op <- par(bty='L')
     plot(z, panel = p, xaxs = 'i', las=1, ...)
     par(op)
}
tmp <- daily[, c('rh_mean','vpd','et','t_mean')]
yl  <- c('RH','VPD','ET','T')
plot_thresh(tmp, yvar='rh_mean', yval=70, 'less', ylab=yl)
rm(tmp, yl)

### TODO: cumulative time above 70% RH may approximate active growth time
tmp <- daily[, c('rh_mean','t_mean')]
active <- (tmp$rh_mean > 70 & tmp$t_mean > 5 & tmp$t_mean < 20)
rx  <- zoo(cumsum(active), index(tmp))
ry  <- zoo(active, index(tmp))
tmp <- cbind(tmp, growth = rx, activeday = ry)
plot_thresh(tmp, yvar='rh_mean', yval=70, 'less')







####################################################################
### wavelets coherence plot with each time-series
`plot_wc` <- function(z, x1, x2, ylaba=x1, ylabb=x2,
                      nrands=99, legend=TRUE, ...){
     if(!all(is.character(x1),is.character(x2)))
          stop('x1,x2 must be character')
     ii <- as.POSIXlt(index(z))$yday
     t1 <- cbind(ii, as.numeric(z[,x1]))
     t2 <- cbind(ii, as.numeric(z[,x2]))
     w  <- wtc(t1, t2, nrands=nrands, sig.level=0.95)
     layout(matrix(1:3,3,1),hei=c(1,1,2))
     op <- par(oma=c(0,0,0,0), mar=c(0,4.1,0,.1), las=1, bty='l')
     plot(z[,x1], ylab=ylaba, xaxs='i', xaxt='n')
     plot(z[,x2], ylab=ylabb, xaxs='i', xaxt='n', col=2)
     par(mar=c(4,4.1,.1,.1))
     plot(w, fill.cols=inferno(99, begin=0.2), plot.cb=F,
          plot.phase=T, lwd.sig=.5, col.sig='#00000070',
          arrow.len=.07, xlab='Julian day', ...)
     if(legend){
          txt <- c('Arrows:
                   Right = a,b in phase
                   Left = a,b in anti-phase
                   Up = a leads b
                   Down = b leads a')
          text(x = graphics::grconvertX(0.01, from='npc', to='user'),
               y = graphics::grconvertY(0.23, from='npc', to='user'),
               labels=txt, adj=0, font=2, col='#00000099')
     }
     par(op)
}
# tiff('C:/Users/Rob/Desktop/fig_01_wavelets.tif', hei=5, wid=5,
#      unit='in', res=500, compr='lzw+p')
plot_wc(daily, 't_mean', 'rh_mean',
        ylaba='Temperature (\u00B0C)', ylabb='Rel humidity (%)')
# dev.off()

plot_wc(daily, 't_mean', 'rh_mean')
plot_wc(daily, 'vpd', 'rh_mean')
plot_wc(daily, 'vpd', 't_mean')
plot_wc(daily, 'vpd', 'et')
plot_wc(daily, 't_mean', 'et')

####   END   ####