# P2P Method (Point to Point)

#Script to analyze two poll tracks and find the points
#closest to each scan (RSL and VL), generating the discrepancies.

### Data path

setwd("C:/Users/Ítalo/Documents/2_Pesquisa_Doutorado/1_TESE/2_TESE_AFONS
O/Capítulo3/Base de Dados/Dados de Estudo")
getwd()

### List of packages to be used in the Script

pkg <- c("rgeos", "raster", "rgdal", "nabor")
sapply(pkg, require, character.only=TRUE)

# pts_near Function

pts_near <- function (LRS,LV,lim_dist=0, b=0)
{
#RSL - Regular Survey Line  (SpatialPointsDataframe data type)
#VL - Verification Line (SpatialPointsDataframe data type)
#lim_dist - maximum distance (in the same unit as the RSL)
#where the closest point can be considered valid, the initial value will be the diagonal 
#of the rectangle surrounding the overlap between the two strips
#b - the value of an additional threshold, i.e. a buffer at the intersection of the tracks
#the default value of b==0, that is, no buffer will be applied

# Extracting the RLS range limit

LRS_points_lim <- LRS@coords[chull(LRS@coords[,1:2]),1:2]
LRS_pol <- SpatialPolygons(list(Polygons(list(Polygon(LRS_points_lim)),
ID=1)))
#Extraindo o limite da faixa LV
LV_points_lim <- LV@coords[chull(LV@coords[,1:2]),1:2]187
LV_pol <- SpatialPolygons(list(Polygons(list(Polygon(LV_points_lim)), ID=1)))

# Intersection between the RLS and VL bands

Int <- intersect(LRS_pol,LV_pol)

# Generating a buffer at the intersection

Int_B<- gBuffer(Int, width = b)
crs(Int_B)<-crs(LRS)

# Extracting the RLS and VL points that are within the intersection between the tracks

LRS_int <- LRS[Int_B,]
LV_int <- LV[Int_B,]

#If the value of the limit distance to find the closest point is not informed, 
#the diagonal of the intersection range will be assigned as the limit

if (lim_dist==0)
lim_dist <- sqrt(sum((diff(t(as.matrix(extent(LRS_int)))))^2))

#Find the closest points from the reference range - VL, 
#find the closest point in the range to be analyzed - RSL

near <- knn((LV_int@coords[,1:2]), (LRS_int@coords[,1:2]),k=1)

### Which point is with a distance above the LIMIT DISTANCE (lim_dist)

pts_lim <- which(near$nn.dists>lim_dist)

### Update the NEAR file with the values of

#Z: points quota
#dp: discrepancy between the reference (VL) and the data (RSL), i.e. dp = VL - RSL

near <- list(near,
array(LV_int@data$Z[near$nn.idx], dim = dim(near$nn.idx)), #quota of nearby points obtained in the LV file
as.numeric(LRS_int@data$Z), #RSL file point quota
coordinates(LRS_int), #Points quota
dp = as.numeric(LV_int@data$Z[near$nn.idx[,1]])-
as.numeric(LRS_int@data$Z))
names(near) <- c("pontos", "Z_in_ref","Z_in_Dados","coordenadas","dp")
LRS_int@data <- cbind(LRS_int@data,Z2=near$Z_in_ref,dz=near$dp)

### Exclude points that have distances above the DISTANCE LIMIT

LRS_int<- LRS_int[-pts_lim,]
return(LRS_int)188
} #End of the pts_near function

### Database reading (shp file)

### Input data:

#a) Regular Survey Line (RSL) in shp format
#b) Verification Line (LV) in shp format

aa <- Sys.time()

### Reading data in shp format

Faixa_dados<-readOGR(dsn=".",layer="LRS")
Faixa_Ver<-readOGR(dsn=".",layer="LV")

### Run pts_near function

limite <- 0.01
buffer <- 0
dp <- pts_near(Faixa_dados, Faixa_Ver, limite, buffer)

### Plot database

windows(8,6,title="Database Analyzed")
par(mfrow=c(1,1), family="serif")
plot(Faixa_dados, pch=1, col=1,
xlab="E (m)", ylab= "N (m)", main="Discrepancies")
points(Faixa_Ver, pch=1, col=4)
points(dp, pch=19, col=2)
legend('bottomleft',legend=c('Regular Survey Line','Verification Line',
"Discrepancies"),col=c(1, 4, 2),pch=c(1, 1, 19))

### Exporting information:

sink("Results.txt", type="output", append=T)
cat("##### Doctoral thesis #####\n Prof. Ítalo O. Ferreira\n e-mail:
italo.ferreira@ufv.br\n\n Method for Obtaining Discrepancies in Scanning Sounds","\n",
"------------------------------------------------------","\n",
"RSL:" ,length(Faixa_dados),'Depths' ,"\n",
"VL:" ,length(Faixa_Ver),'Depths' ,"\n",
"Maximum Distance:" ,limite,"meters" ,"\n",
"Buffer:" ,buffer,"meters","\n",
"Number of Discrepancies:",length(dp@data$dz),"\n",
"------------------------------------------------------","\n",
fill=F)
sink()
shell.exec("Results.txt")

### Check columns for deletion

names(dp)

### Generate dps file in txt format (X, Y, Z, dz)

write.table(dp@data[,-c(4,5)], "dp.txt", dec=",")

### Generate dps file in shp format (X, Y, Z, dz)

writeOGR(dp[,-c(4,5)], dsn=".", layer="dp", driver="ESRI Shapefile")
Sys.time() - aa
