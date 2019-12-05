######################################################################
# Doug Glavich's R1-R3 lichen scores for east of Cascades
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 30 Mar 2018
##  CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)

###   preamble
setwd('~/_prj/7_elon/fig/')
rm(list=ls())
pkg <- c('maps','vegan','ecole')
has <- pkg %in% rownames(installed.packages())
if(any(!has))install.packages(pkg[!has])
lapply(pkg, require, character.only=TRUE) ; rm(pkg, has)

###   load data
# lbs <- c('150','1000','2000','3000','4000','5000')
# cols <- inferno(7, alpha=0.7)
d <- read.csv(file='~/_prj/7_elon/data/glavich_30Mar2018.csv',
              header=T, stringsAsFactors=F)[,-1]
names(d) <- tolower(names(d))

###   delta NMS scores
d$d_air <- (d$rd3air - d$rd1air)
d$d_cli <- (d$rd3clim - d$rd1clim)

###   plot change R3-R1
plot(d$rd1air, d$rd3air) ; abline(0,1)
plot(d$rd1clim, d$rd3clim) ; abline(0,1)

###   plot NMS successional arrows
s3 <- scores(scores(rbind(cbind(d$rd1air, d$rd1clim),
                          cbind(d$rd3air, d$rd3clim))))
ordiplot(s3, type='n')
points(s3, col=colvec(grp), pch=16, cex=0.6)
ordiarrows(s3, lev=156, repl=2, length = 0.07, angle=15)

###   map one
par(las=1, bty='l', oma=c(0,0,0,0), mar=c(4,4,0,0))
map('county', regions='oregon')
points(d$londd, d$latdd, col=colvec(d$d_cli), pch=16)

###   map all
par(mfrow=c(2,3), las=1, bty='l', oma=c(0,0,0,0), mar=c(1,1,0,0))
nm <- c('rd1air','rd3air','r1r3air','rd1clim','rd3clim','r1r3clim')
for(i in 1:length(nm)){
     map('county', regions='oregon')
     points(d$londd, d$latdd, col=colvec(d[,nm[i]]), pch=16)
     title(nm[i])
}

###  map delta
par(mfrow=c(1,2), las=1, bty='l', oma=c(0,0,0,0), mar=c(0,0,0,0))
nm <- c('d_air','d_cli')
for(i in 1:length(nm)){
     map('county', regions='oregon')
     points(d$londd, d$latdd,
            col=colvec(d[,nm[i]],
                       pal=colorRampPalette(
                            c(rgb(1,0,0,0.7),
                              rgb(1,1,1,0.1),
                              rgb(0,0,1,0.7)),
                            alpha=TRUE)(99)),
            pch=16, cex=0.7)
     title(nm[i])
}


###   map climate scores delta
# tiff('glavich.tif',hei=3.5,wid=3.5,unit='in',res=500,compr='lzw+p')
par(mfrow=c(1,1), las=1, bty='l', oma=c(0,0,0,0), mar=c(0,0,0,0))
map('county', regions='oregon')
points(d$londd, d$latdd, pch=16, cex=0.7,
       col=colvec(d$d_cli, pal=colorRampPalette(
     c(rgb(1,0,0,.7),rgb(1,1,1,.1),rgb(0,0,1,.7)),alpha=T)(99)))
title('Climate scores')
# dev.off()

###   get DEM to plot beneath
require(elevatr)
require(raster)
m <- map('county', regions='oregon')
bb <- data.frame(x=c(-124.5664, -116.4633), y=c(41.8920, 46.2938))
p <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
r <- get_elev_raster(loc=bb, prj=p, z = 7, src = "aws")
r <- crop(r, extent(bb))
# tiff('glavich_dem.tif',hei=4.5,wid=4.5,unit='in',res=500,compr='lzw+p')
plot(r^0.8, col=inferno(99,alpha=.3,dir=-1),legend=F,axes=F,box=F,
     oma=c(0,0,0,0), mar=c(0,0,0,0))
lines(m, add=T, col='grey30')
points(d$londd, d$latdd, pch=16, cex=0.8,
       col=colvec(d$d_cli, pal=colorRampPalette(
            c(rgb(1,0,0,.9),rgb(1,1,1,.1),rgb(0,0,1,.9)),alpha=T)(99)))
# dev.off()


####   END   ####
