library(rgeos)

##  First project data into a planar coordinate system (here UTM zone 32)
utmStr <- "+proj=utm +zone=%d +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
crs <- CRS(sprintf(utmStr, 32))
pUTM <- spTransform(p, crs)
ptsUTM <- spTransform(pts, crs)

## Set up container for results
n <- length(ptsUTM)
nearestCantons <- character(n)

## For each point, find name of nearest polygon (in this case, Belgian cantons)
for (i in seq_along(nearestCantons)) {
    nearestCantons[i] <- pUTM$NAME_2[which.min(gDistance(ptsUTM[i,], pUTM, byid=TRUE))]
}

## Check that it worked
nearestCantons
# [1] "Wiltz"            "Echternach"       "Remich"           "Clervaux"        
# [5] "Echternach"       "Clervaux"         "Redange"          "Remich"          
# [9] "Esch-sur-Alzette" "Remich"   

plot(pts, pch=16, col="red")
text(pts, 1:10, pos=3)
plot(p, add=TRUE)
text(p, p$NAME_2, cex=0.7)
