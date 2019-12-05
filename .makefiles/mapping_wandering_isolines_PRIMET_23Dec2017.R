
# ###   preamble
setwd('~/_prj/7_eon/fig/')
rm(list=ls())
pkg <- c('zoo','maps','viridis','akima','ggplot2','gridExtra','grid',
         'magick')
has <- pkg %in% rownames(installed.packages())
if(any(!has))install.packages(pkg[!has])
lapply(pkg, require, character.only=TRUE)
rm(pkg, has)

###   load data
lbs <- c('150','1000','2000','3000','4000','5000')
cols <- inferno(7, alpha=0.7)
d <- read.csv(file='~/_prj/7_eon/data/MS00101_v8.csv',
              header=T, na.strings='NaN', stringsAsFactors=F)[,-1]
d <- d[d$HEIGHT == '150',]  # at DBH
d <- d[!is.na(d$AIRTEMP_MEAN_DAY),]
dt <- as.POSIXct(d[,'DATE'], format='%Y-%m-%d', tz='GMT')
dz <- zoo(d[grepl('AIRTEMP_MEAN_DAY', names(d))], dt)
## remove duplicated indexes by median
dz <- aggregate(dz, index(dz), median)
dt <- index(dz)
names(dz) <-'dayt'
plot(dz, screens=1, col=inferno(6, alpha=0.4), ylim=c(-15, 30),
     las=1, bty='l')

###   MAT per year at the PRIMET site
bys <- format(dt, "%Y")
mat <- aggregate(dz, by=bys, FUN=mean)
mat

###   ClimateNA data for grid points
pts <- read.csv(file='~/_prj/7_eon/data/dd_out.csv',
                header=T, na.strings='-9999',
                stringsAsFactors=F)[,c(3,4,6:42,175)]
names(pts)[1:2] <- c('latdd', 'londd')
pts <- na.omit(pts)
pts$td <- pts$Tave07 -  pts$Tave01      # calculate continentality
pts$TDMAT <- pts$TD / pts$MAT           # calculate CVish
pts$TDMAT[which(pts$TDMAT > 15)] <- NA  # troublemaker pts

# # map it
# eb <- element_blank()
# theme_set(theme_classic() +
#                theme(#legend.position='none',
#                     legend.position=c(.9,.3),
#                     legend.background=eb,
#                     axis.title=eb, axis.text=eb,
#                     axis.ticks=eb, axis.line=eb,
#                     panel.background=eb,
#                     plot.background=eb,
#                     panel.border=element_rect(fill=NA)))

###   plot: wandering isolines
lbs <- index(mat)
z <- 'MAT'
gr <- akima::interp(x=pts$londd, y=pts$latdd,
                    z=pts[,which(colnames(pts)%in%z)],
                    nx=120, ny=120)
ref <- 1.0
incr <- mat
cols <- inferno(45)
ln   <- length(cols)
brk  <- as.vector(ref + incr)
rnk  <- rank(brk, ties.method ="first")
cols <- cols[rnk]
tiff('C:/Users/Rob/Desktop/Rplot.tif',wid=6.5*3,hei=5*3,unit='in',
     res=400, compr='lzw')
layout(matrix(1:45, 5,9, byrow=TRUE))
for(i in 1:45){
     mcols <- rep('#00000000', ln)
     mcols[i] <- cols[i]
     map('state', region=c('oreg','washing','idaho'),mar=rep(0,4))
     contour(gr, col=mcols, levels=brk, drawlabels=F, add=T)
     title(lbs[i])
}
dev.off()

###   plot: save each annual plot to tmp folder
layout(1)
for(i in 1:45){
     f <- paste0('C:/Users/Rob/Desktop/tmp/p_', lbs[i], '.tif')
     tiff(f,wid=3.5,hei=3,unit='in',res=400, compr='lzw')
     mcols <- rep('#00000000', ln)
     mcols[i] <- cols[i]
     map('state', region=c('oreg','washing','idaho'),mar=rep(0,4))
     contour(gr, col=mcols, levels=brk, drawlabels=F, add=T)
     title(lbs[i])
     dev.off()
}

### animation: wandering isolines of MAT (relative to 30-y normals)
imgs <- list.files('C:/Users/Rob/Desktop/tmp/')
imgs <- paste0('C:/Users/Rob/Desktop/tmp/', imgs)
imglst <- image_read(imgs[1]) # initiate a list of images
for(i in 2:45){
     imglst <- c(imglst, image_read(imgs[i]))
}
imglst <- image_scale(imglst, "500x400")
image_info(imglst)
ianim <- image_animate(imglst, fps=1, dispose="previous")
image_write(ianim, path = "C:/Users/Rob/Desktop/MAT_PRIMET_anim.gif",
            format = "gif")

###   END   ####