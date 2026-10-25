library(geosphere)
library(mapview)
library(sf)

# loading data
CMRdata<-read.csv("tullia_Sumava_merged.csv")
summary(CMRdata)
dim(CMRdata)

# to reformat date and time to date values
# the time is currently stored in format "2022/06/14 08:08:40.000"
CMRdata$cas<-CMRdata$time
CMRdata$time<-as.POSIXlt(CMRdata$date)
summary(CMRdata)
CMRdata$time<-paste(CMRdata$time,round(CMRdata$cas),sep=" ")
CMRdata$time<-as.POSIXlt(CMRdata$time,format = "%Y-%m-%d %H")

#if not ordered according to dates and time
CMRdata<-CMRdata[order(CMRdata$time),]

behavior_sezona_sex<-table(CMRdata$date,CMRdata$behaviour,CMRdata$sex)
sink("behavior_vs_sezona_vs_sex.txt")
behavior_sezona_sex
sink()

# remove rows with no coordinates and remove data from Budweis
CMRdata$lat<-CMRdata$Y
CMRdata$long<-CMRdata$X
CMRdata<-CMRdata[!is.na(CMRdata$X),]
CMRdata<-CMRdata[CMRdata$X>14.5,]
dim(CMRdata)
summary(CMRdata)
CMRdata<-CMRdata[!is.na(CMRdata$time),]

# extract the captured individuals; the individual codes are stored in NAME column!

individuals<-CMRdata[!duplicated(CMRdata$NAME),]
individuals$code<-individuals$NAME
nrow(individuals)
individuals$nrecap<-0
individuals$dist<-NA
individuals$time.between.captures<-NA

# plot the points, first captures in black, recaptures in red
plot(CMRdata$long,CMRdata$lat,pch=15,col=CMRdata$recap+1)
plot(CMRdata$long,CMRdata$lat,pch=15)

CMRdata.gis<-st_as_sf(CMRdata,coords = c(11,10))
st_crs(CMRdata.gis) <- 4326
mapview(CMRdata.gis)

# to calculate the individual distances and draw the trajectories on the map
for (i in 1:nrow(individuals)) {
  one<-CMRdata[CMRdata$code==individuals$code[i],]
  if (nrow(one)>1) {

    distance<-0
    time.distance<-0
    for (y in 2:nrow(one)) {
      distance<-distance+distGeo(c(one$long[y-1],one$lat[y-1]),c(one$long[y],one$lat[y]))
      segments(one$long[y-1],one$lat[y-1],one$long[y],one$lat[y])
      time.distance<-time.distance+as.numeric(paste(difftime(one$time[y],one$time[y-1],units="days")))
    }
    print(paste("individual",individuals$code[i],"captured",nrow(one),"times,","distance=",distance))
  individuals$nrecap[i]<-nrow(one)
  individuals$dist[i]<-distance
  individuals$time.between.captures[i]<-time.distance
  }
}

# Summary
summary(individuals)

# saving the results
write.table(individuals,"CMR_recaptures_and_distances_tullia_Sumava.txt",sep="\t")

# save the trajectorises as spatial line
plot(CMRdata.gis[6])

one<-CMRdata.gis[CMRdata.gis$NAME==individuals$code[5],]
cara<-st_linestring(as.matrix(as.data.frame(one)[,1:2]))

for (i in 6:nrow(individuals)) {
  one<-CMRdata.gis[CMRdata.gis$NAME==individuals$code[i],]
  if (nrow(one)>1) {
    
   cara2<-st_linestring(as.matrix(as.data.frame(one)[,1:2]))
   cara<-c(cara,cara2)
   plot(cara,add=T)
    print(paste("Individual",individuals$code[i],"captured",nrow(one),"times,","distance=",distance))
  }
}

plot(cara)
plot(CMRdata.gis[6],add=T)
terra::north()
terra::sbar(10)

cara<-st_sfc(geometry=cara)
st_crs(cara)<-st_crs(CMRdata.gis)
mapview(CMRdata.gis)+mapview(cara)

st_write((cara),"cary.shp")
st_write((CMRdata.gis),"body.shp")

#### for different sex
unique(individuals$sex)
individualsM<-individuals[individuals$sex=="m",]
individualsF<-individuals[individuals$sex=="f",]

# histogram showing amount of dispersal distances
hist(individuals$dist,breaks=10,xlim=c(0,2550),xlab="Vzd?lenost",ylab="Frekvence",main="")
hist(individualsM$dist,col="blue",add=T,breaks=10,xlim=c(0,2550)) #males
hist(individualsF$dist,col="red",add=T,xlim=c(0,2550)) #females

# plot showing relationship between number of recaptures and crossed distances
plot(individuals$dist,individuals$nrecap,xlab="distance",ylab="no capture events")
points(individualsM$dist,individualsM$nrecap,pch=16,col="blue") # males
points(individualsF$dist,individualsF$nrecap,pch=16,col="red") # females

# plot showing relationship between crossed distances and time between captures
plot(individuals$dist,individuals$time.between.captures,xlab="distance",ylab="time between recaptures")
points(individualsM$dist,individualsM$time.between.captures,pch=16,col="blue") #males
points(individualsF$dist,individualsF$time.between.captures,pch=16,col="red") #females

test.time.dist<-lm(individuals$time.between.captures~individuals$dist)
title(paste("F=",round(anova(test.time.dist)[1,4],2),"p=",round(anova(test.time.dist)[1,5],3)))
untb.range<-range(individuals$dist[!is.na(individuals$dist)])
untb.range<-c(untb.range[1]:untb.range[2])
lines(untb.range,(untb.range*coef(test.time.dist)[2])+coef(test.time.dist)[1],lwd=2)

test.time.dist<-lm(individualsM$time.between.captures~individualsM$dist)
#title(paste("F=",round(anova(test.time.dist)[1,4],2),"p=",round(anova(test.time.dist)[1,5],3)))
untb.range<-range(individuals$dist[!is.na(individuals$dist)])
untb.range<-c(untb.range[1]:untb.range[2])
lines(untb.range,(untb.range*coef(test.time.dist)[2])+coef(test.time.dist)[1],lwd=2,col="blue")

test.time.dist<-lm(individualsF$time.between.captures~individualsF$dist)
#title(paste("F=",round(anova(test.time.dist)[1,4],2),"p=",round(anova(test.time.dist)[1,5],3)))
untb.range<-range(individuals$dist[!is.na(individuals$dist)])
untb.range<-c(untb.range[1]:untb.range[2])
lines(untb.range,(untb.range*coef(test.time.dist)[2])+coef(test.time.dist)[1],lwd=2,col="red")

# plot with wingwear and distance
plot(individuals$dist,individuals$wingwear,xlab="distance",ylab="initial wing wear",ylim=c(1,4))
points(individualsM$dist,individualsM$wingwear,pch=16,col="blue") #males
points(individualsF$dist,individualsF$wingwear,pch=16,col="red") #females

# saving the information to .RData format
save.image()
