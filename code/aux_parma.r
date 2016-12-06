################################################################################
require(osmar)

osmwayspoly=function(lines,proj4string="+init=epsg:4326") {
  list_of_Lines <- slot(lines, "lines")
  polygon <- SpatialPolygons(lapply(list_of_Lines, function(x) {
    Polygons(list(Polygon(slot(slot(x, "Lines")[[1]], "coords"))),
             ID=slot(x, "ID"))}),proj4string=CRS(proj4string))
  return(polygon)
  
}

gClip <- function(shp, bb, keep.attribs=TRUE) {
  
  if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  
  else b_poly <- as(extent(bb), "SpatialPolygons")
  
  new.shape <- gIntersection(shp, b_poly, byid = T)
  
  if(keep.attribs){
    
    new.attribs <- data.frame(do.call(rbind,strsplit(row.names(new.shape)," ")),stringsAsFactors = FALSE)
    
    new.attrib.data <- shp[new.attribs$X1,]@data
    
    row.names(new.shape) <- row.names(new.attrib.data)
    
    return(SpatialPolygonsDataFrame(new.shape, new.attrib.data))
    
  }else
    
  {
    
    return(new.shape)
    
  }
  
}

osmrelation2sp=function(relosm,proj4string="+init=epsg:4326") {
relosm<- get_osm(relation(relosm), full = TRUE)
lines <- as_sp(relosm, what = 'lines')
lines_closed <- gLineMerge(lines)
list_of_Lines <- slot(lines_closed, "lines")
polygon <- SpatialPolygons(lapply(list_of_Lines, function(x) {
  Polygons(list(Polygon(slot(slot(x, "Lines")[[1]], "coords"))),
           ID=slot(x, "ID"))}),proj4string=CRS(proj4string))
return(polygon)

}

osmwayspoly=function(lines,proj4string="+init=epsg:4326") {
  list_of_Lines <- slot(lines, "lines")
  polygon <- SpatialPolygons(lapply(list_of_Lines, function(x) {
    Polygons(list(Polygon(slot(slot(x, "Lines")[[1]], "coords"))),
             ID=slot(x, "ID"))}),proj4string=CRS(proj4string))
  return(polygon)
  
}

osmway2sp=function(wayosm,proj4string="+init=epsg:4326") {
  relosm<- get_osm(way(wayosm), full = TRUE)
  lines <- as_sp(relosm, what = 'lines')
  lines_closed <- gLineMerge(lines)
  list_of_Lines <- slot(lines_closed, "lines")
  polygon <- SpatialPolygons(lapply(list_of_Lines, function(x) {
    Polygons(list(Polygon(slot(slot(x, "Lines")[[1]], "coords"))),
             ID=slot(x, "ID"))}),proj4string=CRS(proj4string))
  return(polygon)
  
}


osmway2sp=function(lines,proj4string="+init=epsg:4326") {
   lines_closed <- gLineMerge(lines)
  list_of_Lines <- slot(lines_closed, "lines")
  polygon <- SpatialPolygons(lapply(list_of_Lines, function(x) {
    Polygons(list(Polygon(slot(slot(x, "Lines")[[1]], "coords"))),
             ID=slot(x, "ID"))}),proj4string=CRS(proj4string))
  return(polygon)
  
}

get.osm.structs <- function() {
    dat.types <- c("building", "amenity", "waterway", "grass", "natural", "park", 
        "highway", "boundary")
    letters <- sapply(dat.types, function(x) toupper(substr(x, 1, 1)))
    matches <- sapply(letters, function(x) which(letters %in% x))
    # returns list of matches for each letter
    indx <- which(sapply(matches, function(x) length(x) > 1))
    # returns which matches have more than one element. There won't be many of
    # these, so a loop is implemented:
    if (length(indx) > 0) {
        for (i in 1:length(indx)) {
            indx.i <- matches[[indx[i]]]  # Only need to extract first element
            letters[indx.i] <- toupper(substr(dat.types[indx.i], 1, 2))
        }
    }
    dat <- data.frame(cbind(dat.types, letters))
    row.names(dat) <- NULL
    names(dat) <- c("dat.types", "letters")
    return(dat)
}

get.polys <- function(dat, type = "building") {
    datp <- xmlParse(dat)
    dato <- as_osmar(datp)
    if (type == "grass") {
        pids <- find(dato, way(tags(k == "landuse" & v == type)))
    } else if (type == "park") {
        pids <- find(dato, way(tags(v == type)))
    } else {
        pids <- find(dato, way(tags(k == type)))
    }
    pids <- find_down(dato, way(pids))
    sp <- subset(dato, ids = pids)
    if (type == "boundary" | type == "highway") {
        sp <- as_sp(sp, "lines")
    } else {
        sp <- as_sp(sp, "polygons")
    }

    return(sp)
}

get.all.polys <- function() {
    osm.structs <- get.osm.structs()

    require(data.table)
    st0 <- Sys.time()
    cat("Downloading raw data ...\n", rep("-", 50), "\n", sep = "")
    fnames <- NULL
    for (i in 1:dim(osm.structs)[1]) {
        cat("Dowloading ", toString(osm.structs$dat.types[i]), " \t... ", sep = "")
        dat <- get.raw.dat(toString(osm.structs$dat.types[i]))
        st <- timetaken(st0)
        cat(" done; time = ", st, "\n", sep = "")
        fname <- paste("dat", toString(osm.structs$letters[i]), sep = "")
        fnames <- c(fnames, fname)  # used in subsequent processing stage
        assign(fname, dat)
        save(list = c(fname), file = fname)
    }

    cat("\nProcessing raw data ...\n", rep("-", 50), "\n", sep = "")
    st0 <- Sys.time()
    for (i in 1:dim(osm.structs)[1]) {
        cat("Processing ", toString(osm.structs$dat.types[i]), " \t... ", sep = "")
        dat <- get.polys(get(fnames[i]), type = toString(osm.structs$dat.types[i]))
        st <- timetaken(st0)
        cat(" done; time = ", st, "\n", sep = "")
        fname <- paste("poly", toString(osm.structs$letters[i]), sep = "")
        assign(fname, dat)
        save(list = c(fname), file = fname)
    }
}


remove.small.polys <- function(sp = sp, limit = 0.75) {
    areas <- unlist(lapply(slot(sp, "polygons"), function(x) slot(x, "area")))
    ai <- sort(areas)[floor(length(areas) * limit)]
    indx <- which(areas >= ai)
    sp <- sp[indx, ]
    return(sp)
}

get.xylims <- function (poly, aspect1=TRUE)
{
  xrange <- range (unlist (lapply (slot (poly, "polygons"),  
                    function (x) min (slot (x, "labpt") [1]))))
  yrange <- range (unlist (lapply (slot (poly, "polygons"),  
                    function (x) min (slot (x, "labpt") [2]))))
  if (aspect1) {
    if (diff (xrange) > diff (yrange)) {
      r <- diff (xrange) / diff (yrange)
      yrange [1] <- yrange [1] - r * diff (yrange) / 2
      yrange [2] <- yrange [2] + r * diff (yrange) / 2
    } else {
      r <- diff (yrange) / diff (xrange)
      xrange [1] <- xrange [1] - r * diff (xrange) / 2
      xrange [2] <- xrange [2] + r * diff (xrange) / 2
    }
  }
  return (list (x=xrange, y=yrange))
}

plot.polys <- function(poly = poly, xylims = xylims, filedump = FALSE, col = "gray40") {
    if (filedump) {
        png(filename = "London_map.png", width = 1680, height = 945, type = "cairo-png", 
            bg = "white")
    } else {
        x11(width = 14)
    }
    par(mar = c(0, 0, 0, 0))
    plot(NULL, NULL, xlim = xylims$x, ylim = xylims$y, xaxt = "n", yaxt = "n", 
        xlab = "", ylab = "", bty = "n")
    usr <- par("usr")
    rect(usr[1], usr[3], usr[2], usr[4], border = NA, col = "gray20")
    plotfun <- function(i) {
        xy <- slot(slot(i, "Polygons")[[1]], "coords")
        polypath(xy, border = NA, col = col)
    }
    junk <- lapply(slot(poly, "polygons"), plotfun)
}

add.polys <- function(poly = poly, col = "gray40") {
    plotfun <- function(i) {
        xy <- slot(slot(i, "Polygons")[[1]], "coords")
        polypath(xy, border = NA, col = col)
    }
    junk <- lapply(slot(poly, "polygons"), plotfun)
}

add.lines <- function(lines = lines, col = "gray60") {
    plotfun <- function(i) {
        xy <- slot(slot(i, "Lines")[[1]], "coords")
        lines(xy, col = col)
    }
    junk <- lapply(slot(lines, "lines"), plotfun)
}

makemap <- function (filedump=TRUE, streets=FALSE)
{
  osm.structs <- get.osm.structs ()

  for (i in 1:dim (osm.structs) [1]) {
    fname <- paste ("poly", toString (osm.structs$letters [i]), sep="")
    load (fname)
  }

  col.grass <- rgb (100, 120, 100, maxColorValue=255)
  col.water <- rgb (77, 77, 92, maxColorValue=255)
  xylims <- get.xylims (polyBU, aspect1=FALSE)

  # Streets are plotted underneath everything, so first polygons are merely default fillers, to be re-plotted later.
  plot.polys (polyN, xylims, filedump=filedump, col=col.water)
  if (!is.na (streets)) add.lines (polyH, col=streets) # highways   
  #add.lines (polyBO, col="white") # boundaries
  # Then amenities, grass and parks:
  add.polys (polyA, col=col.grass)
  add.polys (polyG, col=col.grass)
  add.polys (polyP, col=col.grass)
  # Then buildings:
  add.polys (polyBU, col="gray40")
  # Water has to be plotted last to go over the top of grass and parks   
  add.polys (polyW, col=col.water)
  add.polys (polyN, col=col.water)
  if (filedump) { junk <- dev.off (which=dev.cur()) }
}

get.raw.dat <- function(type = "building") {
    bbox <- "[bbox=-0.21,51.484,-0.06,51.54]"
    if (type == "building") {
        dat <- getURL(paste("http://open.mapquestapi.com/xapi/api/0.6/", "way[building=*]", 
            bbox, sep = ""))
    } else if (type == "waterway") {
        dat <- getURL(paste("http://open.mapquestapi.com/xapi/api/0.6/", "way[waterway=*]", 
            bbox, sep = ""))
    } else if (type == "natural") {
        dat <- getURL(paste("http://open.mapquestapi.com/xapi/api/0.6/", "way[natural=water]", 
            bbox, sep = ""))
    } else if (type == "grass") {
        dat <- getURL(paste("http://open.mapquestapi.com/xapi/api/0.6/", "way[landuse=grass]", 
            bbox, sep = ""))
    } else if (type == "park") {
        dat <- getURL(paste("http://open.mapquestapi.com/xapi/api/0.6/", "way[leisure=park]", 
            bbox, sep = ""))
    } else if (type == "amenity") {
        dat <- getURL(paste("http://open.mapquestapi.com/xapi/api/0.6/", "way[amenity=*]", 
            bbox, sep = ""))
    } else if (type == "shop") {
        dat <- getURL(paste("http://open.mapquestapi.com/xapi/api/0.6/", "way[shop=*]", 
            bbox, sep = ""))
    } else if (type == "boundary") {
        dat <- getURL(paste("http://open.mapquestapi.com/xapi/api/0.6/", "way[boundary=*]", 
            bbox, sep = ""))
    } else if (type == "highway") {
        dat <- getURL(paste("http://open.mapquestapi.com/xapi/api/0.6/", "way[highway=*]", 
            bbox, sep = ""))
    }
    return(dat)
}

# http://rgeomatic.hypotheses.org/category/osrm


