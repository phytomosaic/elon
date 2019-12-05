###################################################################
#
# Process ELON-like shapefiles to get demography per species
#
#  Goal: use segmented images to assign trackable 'genet' IDs
#
#  Rob Smith, smithr2@oregonstate.edu, Oregon State Univ, 01 Dec 2019
#
#  CC-BY-SA 4.0 License (Creative Commons Attribution-ShareAlike 4.0)
#
###################################################################


### preamble
rm(list=ls())
require(rgdal)
require(rgeos)
require(raster)
# require(ecole)

### COVER = area = units are m^2

# KEOGH data, Montana
setwd('~/_prj/7_elon/elon/_data_raw/data_census/keogh/')

### load the segmented shapefiles
pth      <- './shapefiles/shapefiles'
allfiles <- list.files(pth)  # all quadrats and years
qs <- unique(sapply(strsplit(allfiles,'_'),function(x)x[1])) # quadlist

timebegin <- Sys.time()  # starting time

### loop thru each of m quads
for (m in 1:length(qs)){
        timeiter <- Sys.time()
        message('\ndoing quad ', m, ' of ', length(qs), ' at ', timeiter, '\n')
        f  <- allfiles[grep(paste0(qs[m],'_'), allfiles)] # keep only one quad
        f  <- unique(gsub('\\..*','',f))       # remove format suffix
        f  <- sort(f[grep('_C',f)])            # keep only COVER polygons
        ### read in all shapefiles for that quad
        `readshp` <- function(dsn=pth, layer=f[1]){
                s <- rgdal::readOGR(dsn=dsn, layer=layer)
                names(s@data) <- tolower(names(s@data))
                # s$species <- factor(ecole::clean_text(as.character(s$species)))
                s$objectid <- as.numeric(s$objectid)
                s$year     <- sapply(strsplit(layer, '_'), function(x) x[2])
                s$quad     <- sapply(strsplit(layer, '_'), function(x) x[1])
                return(s)
        }
        names(f) <- f
        s <- as.list(f)  # initialize list of shapefiles
        for ( i in f ) {
                s[[i]] <- readshp(layer=f[i])  # ~3 sec each ...
        }
        ### track genets across all yrs, per Chu and Adler (2014)
        `relabel` <- function(s1, s2, ...){
                # 1 - add 5-cm BUFFER (to t-1 year)
                b1 <- rgeos::gBuffer(s1, width=0.005, byid=TRUE)
                # 2 - calc OVERLAP from t-1 to t
                o <- raster::intersect(s2,b1)
                if(!is.null(o)){ # switch if some overlap exists:
                        # 3 - calc AREA of overlap
                        o@data$overlaparea <- rgeos::gArea(o, byid=TRUE)
                        # 4 - consider overlaps only for CONSPECIFICs
                        o <- o[as.character(o@data$species.1) ==
                                       as.character(o@data$species.2),]
                        # 5 - for each new, ASSIGN label of old polygon w max overlap
                        o@data <- do.call(
                                rbind, lapply(split(o@data, o@data$objectid.1), function(x) {
                                        x$kk <- x$objectid.2[which.max(x$overlaparea)]
                                        x}))
                        s2@data$objectid <- o@data$kk[
                                match(s2@data$objectid,o@data$objectid.1)]
                        # 5 - assign new incremented label to new recruits (= no overlap)
                        isnew <- is.na(s2@data$objectid)
                        s2@data$objectid[isnew] <-
                                max(s2@data$objectid,na.rm=T) + 1:sum(isnew*1)

                } else { # if none overlap (all old died, and all new are recruits):
                        s2@data$objectid <- max(s1@data$objectid,na.rm=T) + 1:length(s2)
                }
                return(s2)
        }
        for (i in 1:(length(s)-1) ) {
                # i <- 2
                cat('iter',i,'of',length(s)-1,'\n')
                s[[i+1]] <- relabel(s[[i]], s[[i+1]]) # ~5 sec each ...
        }
        # collect all years for each quadrat into final data.frame
        d <- do.call(rbind, lapply(s, function(i) {
                i@data
        }))
        write.csv(d,file=paste0('./label_genets/',unique(d$quad),'.csv'))
        rm(s,d,f)
        gc()
}

### timing
timeend <- Sys.time()
timebegin - timeend
message('###    DONE relabeling at ', timeend, '   ###\n')

### knit them all together
rm(list=ls())
setwd('./label_genets/')
lf <- list.files('.')
lf <- lf[!lf %in% c('d_clean.csv')]
d   <- Reduce(rbind, lapply(lf, read.csv))
setwd('../')

### get size and sizeNEXT in proper order for IPMs
# setwd('~/_prj/7_elon/elon/_data_raw/data_census/keogh/')
# d  <- read.csv('./label_genets/d_clean.csv', stringsAsFactors=F)
d  <- d[,!colnames(d) %in% c('X.1','X')]
d$species <- ecole::clean_text(as.character(d$species))
names(d)[names(d) == 'objectid'] <- 'genet'
names(d)[names(d) == 'year']     <- 'yr'
d <- d[,c('genet','quad','yr','species','area')]
d <- d[order(d$quad, d$yr),]
a <- aggregate(area ~ genet+yr+quad, data=d, FUN=sum) # SUM duplicates
a$species <- d$species[match(interaction(a$genet,a$quad,a$yr,sep='_'),
                             interaction(d$genet,d$quad,d$yr,sep='_'))]
d <- a  ;  rm(a)
d$area <- log10(d$area)                 # units = log10 m^2
d <- d[order(d$quad, d$genet, d$yr),]  # order critically important!
row.names(d) <- NULL
### get areanext from next year's area:
# reshape wide, where each column = area in one year
d$uid <- paste0(d$genet,'_',d$quad)  # unique genet x quad id
w <- reshape(d[,c('uid','yr','area')], v.names='area', idvar='uid',
             timevar='yr', direction='wide', sep='_')
# # if one year was skipped, assign MEAN of adjoining years ? ? ?
# `infill` <- function(x) { # returns vector with single NAs smoothed
#         s       <- 2:(length(x)-1) # omit first/last
#         isna    <- s[which(is.na(x[s]))]
#         x[isna] <- (x[isna-1] + x[isna+1]) / 2 # mean of neighbors
#         x   # ! ! ! ! TIMEWARN ! ! ! ! need to vectorize
# }
# ww <- data.frame(w[,1], apply(w[,-1], 1, infill))

### then iteratively append year-pairs
nr <- NROW(w)       # number of unique genets
s  <- 2:(NCOL(w)-1) # column sequence for years
length(s)           # number of year-pairs to consider
m <- data.frame(    # collector data.frame
        matrix(NA, nrow=nr*length(s), ncol=4,
               dimnames=list(NULL,c('uid','yr','area','areanext'))),
        stringsAsFactors=F)
for (j in s) {
        yr    <- gsub('^[^_]*_','',dimnames(w)[[2]][j])
        i     <- (1:nr) + (nr*(j-2))
        m[i,] <- data.frame(w$uid, yr, w[,j], w[,j+1],
                            stringsAsFactors=F)
}
rm(nr,w,s,i,j)
m <- data.frame(do.call(rbind, strsplit(m$uid, '_')), m) # recover genets
names(m)[1:2] <- c('genet','quad')
m <- m[!(is.na(m$area) & is.na(m$areanext)),] # omit if both yrs empty
d <- cbind(m, species=d$species[match(m$uid, d$uid)]) # match species
rm(m)

### match climate data to each record
cli <- read.csv('./daily_climate_data.csv', stringsAsFactors=F)
cli$precip[cli$precip == 't'] <- 0.01
cli$precip <- as.numeric(cli$precip)
# ... aggregate to monthly
cli <- merge(
        aggregate(cbind(tmax=cli$max, tmin=cli$min),
                  by=list(month=cli$month, year=cli$year), mean),
        aggregate(data.frame(pcum=cli$precip),
                  by=list(month=cli$month, year=cli$year), sum)
)
cli <- cli[order(cli$year, cli$month),]
# ... aggregate to annual (consider June to June census period)
cli$year[cli$month%in%c(1:6)] <- cli$year[cli$month%in%c(1:6)] - 1
cli <- cli[!cli$year%in%c(1931,1945),] # drop non-census years
cli <- merge(
        aggregate(cbind(tmax=cli$tmax, tmin=cli$tmin),
                  by=list(year=cli$year), mean),
        aggregate(data.frame(pcum=cli$pcum),
                  by=list(year=cli$year), sum)
)
cli      <- cli[order(cli$year),]
cli$year <- cli$year - 1900
cli$mat  <- (cli$tmax - cli$tmin) / 2 # approximate MAT
d <- cbind(d, cli[match(d$yr, cli$year),c('mat','pcum')])
d$surv <- ifelse(is.na(d$areanext),0,1) # survived if Y1 reappeared Y2
write.csv(d,file='./d_CENSUS_TIDY.csv')
rm(list=ls())

####    END    ####