######################################################################
# HJA vertical microclimate and lichen transplants
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 14 Dec 2017
##  CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

# ###   preamble
setwd('~/_prj/7_elon/fig/')
rm(list=ls())
pkg <- c('zoo','maps','viridis','akima','ggplot2','gridExtra','grid')
has <- pkg %in% rownames(installed.packages())
if(any(!has))install.packages(pkg[!has])
lapply(pkg, require, character.only=TRUE)
# lapply(pkg, citation)  # credit is due these authors
rm(pkg, has)

###   load data
# Data downloaded on 14 Dec 2017 from:
#    https://andrewsforest.oregonstate.edu/sites/default/files/lter/
#         data/weather/portal/MISC/DSCMET/data/dscmet_comps_2017.csv
# Description at:
#    https://andrewsforest.oregonstate.edu/sites/default/files/lter/
#         data/weather/portal/MISC/DSCMET/data/dscmet_comps_2017.html
lbs <- c('150','1000','2000','3000','4000','5000')
cols <- inferno(7, alpha=0.7)
d <- read.csv(file='~/_prj/7_eon/data/dscmet_comps_2017.csv',
              header=T, na.strings='NaN', stringsAsFactors=F)[,-1]
d2018 <- read.csv(file='~/_prj/7_eon/data/dscmet_comps_2018.csv',
                  header=T,na.strings='NaN',stringsAsFactors=F)[-1,-1]
d <- rbind(d, d2018) ; rm(d2018)
dt <- as.POSIXct(d[,1], format='%Y-%m-%d %H:%M:%S', tz='GMT')
dz <- zoo(d[grepl('DSCMET_AIRTEMP_MEAN_', names(d))], dt)
colnames(dz) <- gsub('DSCMET_AIRTEMP_MEAN_','',colnames(dz))
plot(dz, screens=1, col=inferno(6, alpha=0.4), ylim=c(-10, 35),
     las=1, bty='l')
# correct some wonky values during Aug 29 - Sep 08:
rms <- 95600:98500
dz[rms,1] <- NA
# plot(dz[rms,1])
# lines(dz[rms,2], col=2)
# lines(dz[rms,3], col=3)
# lines(dz[rms,4], col=4)
# lines(dz[rms,5], col=5)
# lines(dz[rms,6], col=6)

# # daily temperature variance(SD)
# bys <- as.Date(index(dz))
# # bys <- as.yearmon(index(dz))
# f1 <- function(x, ...){ sd(x, na.rm=TRUE) }
# f2 <- function(x, ...){ aggregate(x, by=bys, FUN=f1) }
# f3 <- function(x, ...){ ifelse(!is.finite(x), NA, x)}
# dailysd <- zoo(sapply(dz,       FUN=f2), unique(bys))
# dailysd <- zoo(sapply(dailysd, FUN=f3), time(dailysd))
# # tiff('./fig_0_daily_t_sd.tif',wid=7.5,hei=3.5,unit='in',res=400,
# #      compr='lzw')
# layout(matrix(c(1,1,1,1,2,3), nrow=1))
# par(mar=c(4,4,0,0))
# plot(dailysd, screens=1, col=cols, ylim=c(0,10), las=1, xlab='Month',
#      ylab='Diurnal temperature SD (°C)')
# tt <- time(dailysd) ;  tx <- seq(tt[1], tt[length(tt)], by='months')
# axis(1, at=tx, labels=as.POSIXlt(tx)$mon, tcl=-0.2, las=2); rm(tt,tx)
# legend('topleft', inset=c(0,0), legend=lbs, lty=1,
#        bty='n', col=cols, cex=0.7)
# par(mar=c(4,0,0,.01))
# plot(1,type='n',xlab='EPDF',ylab='',
#      xlim=c(0,.35), ylim=c(0,10), las=1, yaxt='n')
# for(i in 1:ncol(dailysd)){
#      den <- density(dailysd[,i],na.rm=T, from=0)
#      lines(den$y, den$x, lwd=2, col=cols[i])
# }
# plot(1,type='n',xlab='ECDF',ylab='',
#      xlim=c(0,1), ylim=c(0,10), las=1, yaxt='n')
# for(i in 1:ncol(dailysd)){
#      e <- ecdf(as.numeric(dailysd[,i]))
#      lines(environment(e)$y, environment(e)$x, lwd=2, col=cols[i])
# }
# # dev.off()

# daily temperature range (not really continentality)
bys <- as.Date(index(dz))
# bys <- as.yearmon(index(dz))
f1 <- function(x, ...){ diff(range(x, na.rm=TRUE)) }
f2 <- function(x, ...){ aggregate(x, by=bys, FUN=f1) }
f3 <- function(x, ...){ ifelse(!is.finite(x), NA, x)}
dtr <- zoo(sapply(dz,       FUN=f2), unique(bys))
dtr <- zoo(sapply(dtr, FUN=f3), time(dtr))
# tiff('./fig_1_daily_t_rng.tif',wid=7.5,hei=3.5,unit='in',res=400,
#      compr='lzw')
layout(matrix(c(1,1,1,1,2,3), nrow=1))
par(mar=c(4,4,0,0))
plot(dtr, screens=1, col=cols, ylim=c(0,30), las=1, xlab='Month',
     ylab='Diurnal temperature range (°C)')
tt <- time(dtr) ;  tx <- seq(tt[1], tt[length(tt)], by='months')
axis(1, at=tx, labels=as.POSIXlt(tx)$mon, tcl=-0.2, las=2); rm(tt,tx)
legend('topleft', inset=c(0,0), legend=lbs, lty=1,
       bty='n', col=cols, cex=0.7)
par(mar=c(4,0,0,.01))
plot(1,type='n',xlab='EPDF',ylab='',
     xlim=c(0,.12), ylim=c(0,30), las=1, yaxt='n')
for(i in 1:ncol(dtr)){
     den <- density(dtr[,i],na.rm=T, from=0)
     lines(den$y, den$x, lwd=2, col=cols[i])
}
plot(1,type='n',xlab='ECDF',ylab='',
     xlim=c(0,1), ylim=c(0,30), las=1, yaxt='n')
for(i in 1:ncol(dtr)){
     e <- ecdf(as.numeric(dtr[,i]))
     lines(environment(e)$y, environment(e)$x, lwd=2, col=cols[i])
}
# dev.off()

# *monthly* averaged daily temperature range
bymon <- as.yearmon(index(dtr))
f4 <- function(x, ...){ aggregate(x, by=bymon, FUN=mean, na.rm=TRUE)}
mondtr <- zoo(sapply(dtr,  FUN=f4), unique(bymon))
# tiff('./fig_2_monthlyavg_daily_t_rng.tif',wid=7.5,hei=3.5,unit='in',
#      res=400, compr='lzw')
layout(matrix(c(1,1,1,1,2,3), nrow=1))
par(mar=c(4,4,0,0), mgp=c(2, 1, 0))
plot(mondtr, screens=1, col=cols, ylim=c(0,20), las=1, xlab='Month',
     ylab='Monthly-averaged\nDiurnal temperature range (°C)')
legend('topleft', inset=c(0,0), legend=lbs, lty=1,
       bty='n', col=cols, cex=0.7)
par(mar=c(4,0,0,.01))
plot(1,type='n',xlab='EPDF',ylab='',
     xlim=c(0,.1), ylim=c(0,20), las=1, yaxt='n')
for(i in 1:ncol(mondtr)){
     den <- density(mondtr[,i],na.rm=T, from=0)
     lines(den$y, den$x, lwd=2, col=cols[i])
}
plot(1,type='n',xlab='ECDF',ylab='',
     xlim=c(0,1), ylim=c(0,20), las=1, yaxt='n')
for(i in 1:ncol(mondtr)){
     e <- ecdf(as.numeric(mondtr[,i]))
     lines(environment(e)$y, environment(e)$x, lwd=2, col=cols[i])
}
# dev.off()

# # diff relative to 10 m understory (cant do 1.5 m because too many NA)
# f5 <- function(x){x - dtr$`1000_0_02`}
# diffs  <- zoo(sapply(dtr, FUN=f5), time(dtr)) # ~60s
# layout(1)
# plot(diffs, screens=1, col=cols, ylim=c(0, 6), las=1)
# legend('topleft', inset=c(0,0), legend=colnames(dtr),
#        lty=1, bty='n', col=cols, cex=0.7)
#
# # annual mean values, 2016-2017
# (mndsd <- colMeans(dailysd, na.rm=T))
# (mndtr <- colMeans(dtr, na.rm=T))
# (mndif <- colMeans(diffs, na.rm=T))

# calculate annual continentality, etc
bymon <- as.yearmon(index(dz))
monavg <- zoo(sapply(dz,  FUN=f4), unique(bymon)) # monthly avg T
td <- as.numeric(monavg[11,]) - as.numeric(monavg[4,])
# adjust continentality by mean (sort of like a CV)
mat <- sapply(dz, FUN=mean, na.rm=T)
tdmat <- td/mat

# tiff('./fig_33_continentalityCV.tif',wid=3.5,hei=3,unit='in',
#      res=400, compr='lzw')
par( mar=c(2,4.2,0.1,0.5), cex.axis=0.5, cex.lab=0.6, mfrow=c(1,3) )
plot(td, las=1, pch=16, ylab='Continentality (Aug-Dec T)', xaxt='n')
axis(1, at=1:6, labels=colnames(dz))
plot(mat, las=1, pch=16, ylab='MAT (C)', xaxt='n')
axis(1, at=1:6, labels=colnames(dz))
plot(tdmat, las=1, pch=16, ylab='MAT (C)', xaxt='n')
axis(1, at=1:6, labels=colnames(dz))
# dev.off()

# get ClimateNA continentality data for western US
# # generate points
# pts <- expand.grid(lat=seq(40.5, 49.5, by=0.1),
#                    long=seq(-125.5, -110.0, by=0.2))
# head(pts)
# pts <- data.frame(ID1=rep('.', nrow(pts)),
#                   ID2=rep('.', nrow(pts)),
#                   pts,
#                   el=rep('.', nrow(pts)))
# write.csv(pts, file='C:/Users/Rob/Desktop/dd.csv', row.names=F)
pts <- read.csv(file='~/_prj/7_eon/data/dd_out.csv',
                header=T, na.strings='-9999',
                stringsAsFactors=F)[,c(3,4,6:42,175)]
names(pts)[1:2] <- c('latdd', 'londd')
pts <- na.omit(pts)
pts$td <- pts$Tave07 -  pts$Tave01  # calculate continentality
pts$TDMAT <- pts$TD / pts$MAT     # calculate CVish
pts$TDMAT[which(pts$TDMAT > 15)] <- NA  # troublemaker pts
plot(pts$TD, pts$td)
plot(pts$MAT, pts$TDMAT)


# map it
eb <- element_blank()
theme_set(theme_classic() +
               theme(#legend.position='none',
                     legend.position=c(.9,.3),
                     legend.background=eb,
                     axis.title=eb, axis.text=eb,
                     axis.ticks=eb, axis.line=eb,
                     panel.background=eb,
                     plot.background=eb,
                     panel.border=element_rect(fill=NA)))

# points
p1 <- ggplot(pts, aes(x=londd, y=latdd)) +
     coord_map('albers',lat0=50,lat1=70,
               xlim=c(-125,-112), ylim=c(41,49.4))  +
     geom_point(aes(colour=TD), alpha=1, size=3, shape=16) +
     scale_color_viridis(option='B', alpha=0.8)+
     borders('world', 'canada', colour='grey30') +
     borders('state', colour='white')
p2 <- ggplot(pts, aes(x=londd, y=latdd)) +
     coord_map('albers',lat0=50,lat1=70,
               xlim=c(-125,-112), ylim=c(41,49.4))  +
     geom_point(aes(colour=MAT), alpha=1, size=3, shape=16) +
     scale_color_viridis(option='B', alpha=0.8)+
     borders('world', 'canada', colour='grey30') +
     borders('state', colour='white')
p3 <- ggplot(pts, aes(x=londd, y=latdd)) +
     coord_map('albers',lat0=50,lat1=70,
               xlim=c(-125,-112), ylim=c(41,49.4))  +
     geom_point(aes(colour=log10(TDMAT)),alpha=1,size=3,shape=16) +
     scale_color_viridis(option='B', alpha=0.8)+
     borders('world', 'canada', colour='grey30') +
     borders('state', colour='white')

# contour
`cmap` <-function(data, z, ref, title='X', incr, mult=F,
                  xlb='', ylb='', direction=1, lb=lbs, ...){
     data <- na.omit(data)
     gr <- akima::interp(x=data$londd, y=data$latdd,
                         z=data[,which(colnames(data)%in%z)],
                         nx=120, ny=120)
     gr <- cbind( expand.grid(x=gr$x, y=gr$y), z=as.vector(gr$z) )
     if(mult) { brk <- ref * incr } else { brk <- ref + incr }
     if(direction==(-1)){lb<-rev(lb);reverse<-T}else{reverse<-F}
     ggplot(gr, aes(x, y, colour=z)) + ylab(ylb) + xlab(xlb) +
          coord_map('albers',lat0=50,lat1=70,
                    xlim=c(-125,-112), ylim=c(41,49.4))  +
          stat_contour(data=gr, breaks=brk, size=0.7,
                       aes(x=x,y=y,z=z,colour=factor(..level..)))+
          borders('state') + borders('usa', colour='black') +
          guides(colour=guide_legend(
               title, override.aes=list(size=2.5),
               nrow=1,keyheight=.1,keywidth=1,reverse=reverse))+
          theme(
               legend.position='none'
               # legend.position='top',
               #legend.text=element_text(size=8),
               #axis.title=eb, axis.text=eb,
                #axis.ticks=eb, axis.line=eb,
                #plot.margin=unit(c(-0.5,0,-0.5,-0.5), 'mm')
                )+
          scale_color_viridis(option='B', direction=direction,
                              begin=0,end=.9,discrete=T,labels=lb)
}
v <- sort(td)
v <- v + (v*(v - mean(v))*0.2) # expand it a bit
p4 <- cmap(pts, z='TD', ref=0.1, incr=v, title='Height (cm)')
p5 <- cmap(pts, z='MAT', ref=-0.1, incr=mat, title='Height (cm)')
p6 <- cmap(pts, z='TDMAT', ref=-0.1, incr=tdmat, title='Height (cm)')
# tiff('./fig_4_continentality_regional.tif',wid=4.5,hei=4,unit='in',
#      res=400, compr='lzw')
print(p4)
# dev.off()

# tiff('./fig_05_isolines.tif',wid=9.5,hei=6,unit='in',
#      res=400, compr='lzw')
g <- gpar(fontsize=12)
grid.arrange(nullGrob(), textGrob('Continentality (°C)',gp=g),
             textGrob('MAT (°C)', gp=g),
             textGrob('Continentality / MAT',gp=g),
             textGrob('ClimateNA points',gp=g,rot=90),p1,p2,p3,
             textGrob('Isolines of canopy ht',gp=g,rot=90),p4,p5,p6,
             nullGrob(),nullGrob(),
             layout_matrix=matrix(c(1:12,13,14,14,14),ncol=4,byrow=T),
             heights=c(.4,5,5,.3),
             widths=c(.25,5,5,5),
             padding=unit(0,'line'))
# dev.off()

# values seem biased eastward and upward, so, check nearby met station
# what is MAT at PRIMET (436 m asl)?
d2017 <- read.csv(file='~/_prj/7_eon/data/primet_226_a_5min_2017.csv',
                  header=T, na.strings=c('NaN',''), stringsAsFactors=F)
d2018 <- read.csv(file='~/_prj/7_eon/data/primet_226_a_5min_2018.csv',
                  header=T,na.strings=c('NaN',''),stringsAsFactors=F)[-1,]
dp <- rbind(d2017, d2018) ; rm(d2018)
dtp <- as.POSIXct(dp[,1], format='%m/%d/%Y %H:%M', tz='GMT')
dp <- zoo(dp$`AIRTEMP_MEAN_150`, dtp)
plot(dp)
abline(h=8.29, col=2) # MAT at PRIMET is 8.29 C
# keep <- grepl('AIRTEMP_', names(dp)) & !grepl('Flag_', names(dp))
# dz <- zoo(dp[,keep], dtp)
# plot(dz, screens=1, col=inferno(6, alpha=0.4), ylim=c(-10, 35),
#      las=1, bty='l')
# colMeans(dz, na.rm=T)
# abline(h=8.29) # MAT at PRIMET is 8.29 C
bymon <- as.yearmon(index(dp)) # calculate *monthly* avg T
f4 <- function(x, ...){ aggregate(x, by=bymon, FUN=mean, na.rm=TRUE)}
mon_t_primet   <- zoo(sapply(dp,  FUN=f4), unique(bymon))
# mat_primet <- mean(mon_t_primet)    # MAT
mat_primet <- mean(dp, na.rm=T)
td_primet    <- 19.9644 - (-1.5272) # Continentality = Aug minus Jan
tdmat_primet <- td_primet / mat_primet     # calculate CVish


# ##################################################################
# ##################################################################
# # relative humidity: Cant really do it with these data, since wrong
# #   values included (many negative values and values >100%)
# dz2 <- zoo(d[grepl('DSCMET_RELHUM_MEAN_', names(d)) ], dt)
# colnames(dz2) <- gsub('DSCMET_RELHUM_MEAN_','',colnames(dz2))
# plot(dz2, screens=1, col=inferno(3, alpha=0.4), ylim=c(25, 105),
#      las=1, bty='l')
# hist(dz2$`150_0_01`)  # STOP here, data wonky......
# # (rhvals <- colMeans(dz2, na.rm=T))
# # v <- unlist(rhvals)
# # pts <- read.csv(file='~/_prj/7_eon/data/dd_out.csv',
# #                 header=T, na.strings='-9999',
# #                 stringsAsFactors=F)[,c(3,4,163:174)]
# # names(pts)[1:2] <- c('latdd', 'londd')
# # pts <- na.omit(pts)
# # pts$rh     <- rowMeans(pts[,3:14], na.rm=T)
# # pts$rh_mod <- pts$rh+20
# # p1 <- cmap(pts, z='rh_mod', ref=0.1, incr=v, title='Height (cm)',
# #            lb=c('150', '5000'))
# # p1
# #
# # # *monthly* averaged daily temperature range
# # cols <- inferno(3, alpha=0.7)
# # bymon <- as.yearmon(index(dz2))
# # f4 <- function(x, ...){ aggregate(x, by=bymon, FUN=mean, na.rm=T)}
# # monrh <- zoo(sapply(dz2,  FUN=f4), unique(bymon))
# # # tiff('./fig_5_monthlyavg_rh.tif',wid=7.5,hei=3.5,unit='in',
# # #      res=400, compr='lzw')
# # layout(matrix(c(1,1,1,1,2,3), nrow=1))
# # par(mar=c(4,4,0,0), mgp=c(2, 1, 0))
# # plot(monrh, screens=1, col=cols, ylim=c(0,105),las=1,xlab='Month',
# #      ylab='Monthly-averaged\nRelative humidity (%)')
# # legend('topleft', inset=c(0,0), legend=colnames(monrh), lty=1,
# #        bty='n', col=cols, cex=0.7)
# # par(mar=c(4,0,0,.01))
# # plot(1,type='n',xlab='EPDF',ylab='',
# #      xlim=c(0,.1), ylim=c(0,105), las=1, yaxt='n')
# # for(i in 1:ncol(monrh)){
# #      den <- density(monrh[,i],na.rm=T, from=0)
# #      lines(den$y, den$x, lwd=2, col=cols[i])
# # }
# # plot(1,type='n',xlab='ECDF',ylab='',
# #      xlim=c(0,1), ylim=c(0,105), las=1, yaxt='n')
# # for(i in 1:ncol(monrh)){
# #      e <- ecdf(as.numeric(monrh[,i]))
# #      lines(environment(e)$y, environment(e)$x, lwd=2, col=cols[i])
# # }
# # # dev.off()
# #
# #
# # # VPD
# # dz3 <- zoo(d[grepl('DSCMET_VPD_MEAN_', names(d)) ], dt)
# # colnames(dz3) <- gsub('DSCMET_VPD_MEAN_','',colnames(dz3))
# # plot(dz3, screens=1, col=inferno(3, alpha=0.4), #ylim=c(25, 105),
# #      las=1, bty='l')
# # colMeans(dz3, na.rm=T)


###   END   ####