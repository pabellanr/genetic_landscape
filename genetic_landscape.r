#################################################################################################
## Script to obtain genetic landscapes based on pairwise population genetic divergence         ##
## as measured among multiple collection locations or populations                              ##       
##                                                                                             ##
## Pedro Abellán                                                                               ##
## October 2013                                                                                ##
## Citation: Abellán P, Svenning JC. 2104. Refugia within refugia– patterns in endemism        ##
## and genetic divergence are linked Acrobatto Late Quaternary climate stability in the        ##
## Iberian Peninsula. Biological Journal of the Linnean Society 113: 13–28.                    ##       
#################################################################################################


# load required libraries
library(geosphere)
library(tripack)
library(reshape)
library(raster)
library(gstat)
library(rgeos)

# load locality data 
dat <- read.table("locs.txt", header=TRUE, row.names=1) # table with three columns: IDs, longitude and latitude

# load  genetic distances
gdist <- read.table("gdist.txt", header=TRUE, row.names=1, fill=T) # table with pairwise genetic distances (e.g. from Arlequin)
names(gdist) <- row.names(gdist)

## Obtain a network connecting all collection locations to their nearest neighbors with non-overlapping edges
tri <- tri.mesh(dat[,1],dat[,2])

## Obtain the  midpoints between each connected edge 
t.tri <- triangles(tri)
out <- data.frame()
for (i in 1:dim(t.tri)[1]){
n1 <- t.tri[i,"node1"]
n2 <- t.tri[i,"node2"]
n3 <- t.tri[i,"node3"]
pt12 <- midPoint(c(dat[t.tri[i,"node1"],1], dat[t.tri[i,"node1"],2]), c(dat[t.tri[i,"node2"],1], dat[t.tri[i,"node2"],2]) )
pt13 <- midPoint(c(dat[t.tri[i,"node1"],1], dat[t.tri[i,"node1"],2]), c(dat[t.tri[i,"node3"],1], dat[t.tri[i,"node3"],2]) )
pt23 <- midPoint(c(dat[t.tri[i,"node2"],1], dat[t.tri[i,"node2"],2]), c(dat[t.tri[i,"node3"],1], dat[t.tri[i,"node3"],2]) )
c12 <- data.frame(n1,n2,pt12); names(c12) <- c("n1","n2","x","y")
c13 <- data.frame(n1,n3,pt13); names(c13) <- c("n1","n2","x","y")
c23 <- data.frame(n2,n3,pt23); names(c23) <- c("n1","n2","x","y")
pts <- rbind(c12,c13,c23) 
out <- rbind(out,pts)
}
xx <- unique(out); row.names(xx) <- 1:dim(xx)[1]

## The midpoints between each connected edge are mapped, and values of genetic distances are attached to these points

#  get genetic distances between each population pair
g <- data.frame()
n <- dim(gdist)[1]
for (i in 1:n){
  for(j in i:n){
    result <- data.frame(j,i,gdist[j,i])
    g <- rbind(g, result)
  }
}

pivot <- as.matrix(cast(g, i ~ j))  # obtain a pivot table with pairwise genetic distances

# Attach genetic divergence values to each midpoint
g.dist <- c()
for (i in 1:dim(xx)[1]){
dd <- c(pivot[xx$n1[i], xx$n2[i]], pivot[xx$n2[i], xx$n1[i]])
gd <- mean(dd, na.rm=T)
g.dist <- c(g.dist,gd)
}

output <- cbind(xx,g.dist)

## Obtain genetic landscape: inverse distance weighted interpolation 
r <- raster("ip_mask.asc") # load raster with a mask of study area (1 values for study area and NA's for cells outside it) 
proj4string(r) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") # set CRS of our raster
mg <- gstat(formula = g.dist~1, nmax=12, locations = ~x+y, data=output, set=list(idp = 2))
z <- interpolate(r, mg)

## The genetic landscape is clipped to the extent of the boundaries of the region of analysis
values.z <- getValues(z)
values.r <- getValues(r)
values.z[which(is.na(values.r))] <- NA
z <- setValues(z,values.z)

## Obtain converx hull polygon for localities (sampling extent)
pts <- chull(dat) # returns a vector of integer coordinates
pts2 <- c(pts,pts[1]) # close the 'polygon' by appending the first to the last
hullpts <- dat[pts2, ]
matChull <- as.matrix(hullpts)
polyg = Polygon(matChull)
ps.ch = Polygons(list(polyg),1)
ch_sp0 = SpatialPolygons(list(ps.ch))
ch_sp <- gBuffer(ch_sp0, width=0.008333)
proj4string(ch_sp) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

## The genetic landscape is clipped to the extent of the original network (sampling extent)
cells.ip <- 1:ncell(z)
cells <- cellFromPolygon(z, ch_sp)[[1]]
values <- getValues(z)
ids <- which(cells.ip%in%cells == F)
values[ids] <- NA
range01 <- function(x){(x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))}
values01 <- range01(values)
z2 <- setValues(z,values01)

# Export genetic landscape raster
writeRaster(z2, paste(sp,".asc", sep=""), format="ascii")




