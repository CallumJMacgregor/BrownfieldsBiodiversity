###########################################################################
####   Script for quantifying ex-landfill coverage in survey squares   ####
###########################################################################

##### script setup ####

### Clear the current workspace (DOESN'T CLEAR LOADED LIBRARIES)
rm(list=ls())

### install if necessary and then load the libraries you need

j <- c("rstudioapi","plyr","raster","dplyr")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

### set working directory to script's saved location

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



## in QGIS, I've identified every location where there is an ex-landfill site within a structured survey square
# (I've done this for BBS and NPMS, still awaiting data from WCBS)
# and measured the area of those sites
# however, this doesn't account for the fact that there are often multiple landfill shapes in a single square,
# so all we need to do now is add up how much landfill there is in each square in total

#### BBS

# read in the data summary

BBS <- read.csv("BTO BBS/BBS_landfill_intersect_cleaned.csv")
summary(BBS)

# for example, we can see in the summary that square SU5863 has some 10 separate bits of landfill in it!

BBS_sum <- ddply(BBS, .(GRIDREF, Grnd_Tt), summarise,
                 LANDFILL_AREA = sum(Intersection_area))

summary(BBS_sum)

# now we want to calculate the percentage of each square which is brownfield
# the area of each square is 1,000,000 m2

BBS_sum$PERC <- BBS_sum$LANDFILL_AREA*100/1000000

# and count (by picking out) how many squares are > 5% ex-landfill by area
BBS_good <- BBS_sum[which(BBS_sum$PERC > 5), ]




#### NPMS

# read in the data summary

NPMS <- read.csv("NPMS/Analysis/NPMS_landfill_intersect_cleaned.csv")
summary(NPMS)

# for example, we can see in the summary that square SU5863 has some 10 separate bits of landfill in it!

NPMS_sum <- ddply(NPMS, .(NPMS_sq, RECORDS), summarise,
                 LANDFILL_AREA = sum(Intersection_area))

summary(NPMS_sum)

# now we want to calculate the percentage of each square which is brownfield
# the area of each square is 1,000,000 m2

NPMS_sum$PERC <- NPMS_sum$LANDFILL_AREA*100/1000000

# and count (by picking out) how many squares are > 5% ex-landfill by area
NPMS_good <- NPMS_sum[which(NPMS_sum$PERC > 5), ]




#### WCBS

# read in the data summary

WCBS <- read.csv("WCBS/Analysis/WCBS_landfill_intersect_cleaned.csv")
summary(WCBS)

# for example, we can see in the summary that square SU3524 has some 9 separate bits of landfill in it!

WCBS_sum <- ddply(WCBS, .(WCBS_sq, VISITS), summarise,
                  LANDFILL_AREA = sum(Intersection_area))

summary(WCBS_sum)

# now we want to calculate the percentage of each square which is brownfield
# the area of each square is 1,000,000 m2

WCBS_sum$PERC <- WCBS_sum$LANDFILL_AREA*100/1000000

# and count (by picking out) how many squares are > 5% ex-landfill by area
WCBS_good <- WCBS_sum[which(WCBS_sum$PERC > 5), ]







#### merge

# it will be interesting to see whether there are any brownfield-heavy squares that are covered by multiple recording schemes
# so the next step is to merge these dataframes above

# first we need to tidy some things up

colnames(NPMS_sum)[1] <- "GRIDREF"
colnames(BBS_sum)[2] <- "RECORDS"
colnames(WCBS_sum)[1] <- "GRIDREF"
colnames(WCBS_sum)[2] <- "RECORDS"

BBS_sum$SCHEME <- factor("BBS")
NPMS_sum$SCHEME <- factor("NPMS")
WCBS_sum$SCHEME <- factor("WCBS")

summary(BBS_sum)
summary(NPMS_sum)
summary(WCBS_sum)

# now we can bind

combined_sum <- rbind(BBS_sum,NPMS_sum,WCBS_sum)

summary(combined_sum)


# and summarise how many schemes each square is covered by, and which

combined_sum$COUNT <- 1

combined <- ddply(combined_sum, .(GRIDREF), summarise,
                  SCHEMES = sum(COUNT),
                  BBS = ifelse("BBS" %in% SCHEME, T, F),
                  NPMS = ifelse("NPMS" %in% SCHEME, T, F),
                  WCBS = ifelse("WCBS" %in% SCHEME, T, F))

summary(combined)


## sadly only one square is covered by all three schemes, and it's only got ~2% brownfield cover

## a simpler option (and one underpinned by better data) would be to do bipartite butterfly-hostplant networks
# so what's the overlap between NPMS and WCBS specifically?

combined_bipartite <- rbind(NPMS_sum, WCBS_sum)
summary(combined_bipartite)

# unfortunately only one brownfield square, so very sadly a network approach mightn't work well here
# just because this would be interesting to play around with, I'll check overlap of all squares a little later on...



##### land cover

# now we want to introduce data from the Land Cover Map 2015
# and use it to identify dominant land cover in each square,
# and then use this to identify the comparison squares for each brownfield square

# first read in the LandCoverMap

lcm2015 <- raster("Spatial data/Download_LCM2015_York_only_1505408/lcm-2015-25m_3522330/lcm2015gb25m.tif")
lcm2015

plot(lcm2015)

### BBS

# now, we want to generate a small subset of this data for every site
# but to do so we need to read in every site (i.e. even the non-brownfield ones)

BBS_all <- read.csv("BTO BBS/BBSsquares_loc.csv")

summary(BBS_all)

# we can label target brownfield squares in this file before moving on
BBS_target <- BBS_sum[which(BBS_sum$PERC > 5), ]

BBS_all$TARGET <- ifelse(BBS_all$GRIDREF %in% BBS_target$GRIDREF, T, F)


# now, for every square in this list, we need to extract a few items into a new dataframe
# that's 6574 squares so it might take a few minutes!

# we need the function that calculates the mode

getmode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}


# create a list to loop over and seed an output frame

BBS_all_squares <- levels(droplevels(BBS_all$GRIDREF))
n <- 0

BBS_landuse <- data.frame()

## and now begin the loop

for (x in BBS_all_squares){
  
  # first track progress
  n <- n + 1
  if (n == round(n, -2)){
    print(n)
  }
  
  # next, pull out square details
  square <- BBS_all[which(BBS_all$GRIDREF == x), ]
  out <- square[which(square$CORNER == "SW"), ]
  
  # create a polygon from the square's coordinates
  poly <- list(Polygons(list(Polygon(square[,5:6])), ID = x)) 
  SP <- SpatialPolygons(poly)
  
  # extract the LCM class for each 25x25m cell in the square
  cells <- extract(lcm2015, SP)
  
  # get rid of cells with value 0 (i.e. the sea)
  cells[[1]] <- cells[[1]][!cells[[1]] == 0]
  
  out$LANDUSE <- getmode(cells[[1]])
  BBS_landuse <- rbind(BBS_landuse, out)
  
}


summary(BBS_landuse)

# this took a while - do an out and in to preserve

write.csv(BBS_landuse, "BTO BBS/BBS_landuse.csv", row.names = F)

BBS_landuse <- read.csv("BTO BBS/BBS_landuse.csv")



## next we want to identify the non-brownfield neighbours to each brownfield site
# this time we only need to loop over the brownfield sites, and within each one, the non-brownfield sites

BBS_target_squares <- levels(droplevels(BBS_target$GRIDREF))
BBS_nontarget_squares <- BBS_all_squares[which(!(BBS_all_squares %in% BBS_target_squares))]

n <- 0

BBS_neighbours <- data.frame()


# do it in a loop

for (x in BBS_target_squares){
  n <- n + 1
  if (n == round(n, -1)){
    print(n)
  }
  
  square <- BBS_landuse[which(BBS_landuse$GRIDREF == x), ]
  
  pairdists <- data.frame()
  
  # now an internal loop calculates distance to every other square
  for (y in BBS_nontarget_squares){
    neighbour <- BBS_landuse[which(BBS_landuse$GRIDREF == y), ]
      
    SQUARE <- x
    NEIGHBOUR <- y
      
    SQUARE_use <- square$LANDUSE[1]
    NEIGHBOUR_use <- neighbour$LANDUSE[1]
      
    eastdist <- square$EASTING - neighbour$EASTING
    northdist <- square$NORTHING - neighbour$NORTHING
      
    PAIR_DIST <- sqrt(eastdist^2 + northdist^2)
      
    back <- cbind(SQUARE, NEIGHBOUR, SQUARE_use, NEIGHBOUR_use, PAIR_DIST)
    pairdists <- rbind(pairdists, back)
  }
  
  pairdists$PAIR_DIST <- as.numeric(as.character(pairdists$PAIR_DIST))
  
  # split into paired and different land use
  
  pairdists_matched <- pairdists[which(as.character(pairdists$SQUARE_use) == as.character(pairdists$NEIGHBOUR_use)), ]
  pairdists_other <- pairdists[which(as.character(pairdists$SQUARE_use) != as.character(pairdists$NEIGHBOUR_use)), ]
  
  # now take the top row of the matched frame and the top 3 of the other one
  
  matched_neighbour <- top_n(pairdists_matched, -1, pairdists_matched$PAIR_DIST)
  unmatched_neighbours <- top_n(pairdists_other, -3, pairdists_other$PAIR_DIST)
  
  neighbours <- rbind(matched_neighbour, unmatched_neighbours)
  
  BBS_neighbours <- rbind(BBS_neighbours, neighbours)
  
}

# this took a while - do an out and in to preserve

write.csv(BBS_neighbours, "BTO BBS/BBS_neighbours.csv", row.names = F)

BBS_neighbours <- read.csv("BTO BBS/BBS_neighbours.csv")


## finally, generate a non-redundant list of sites of interest

BBS_good$CLASS <- "BROWNFIELD"
BBS_brownfields <- BBS_good[,c(1,5)]

BBS_neighbours$COUNT <- 1
BBS_pairs <- ddply(BBS_neighbours, .(NEIGHBOUR), summarise,
                   PAIRS = sum(COUNT))

BBS_pairs$CLASS <- "PAIR"
BBS_pairs <- BBS_pairs[,c(1,3)]
colnames(BBS_pairs)[1] <- "GRIDREF"

BBS_required <- rbind(BBS_brownfields, BBS_pairs)

# write this list out

write.csv(BBS_required, "BTO BBS/BBS_requested.csv", row.names = F)




### NPMS

# now, we want to generate a small subset of this data for every site
# but to do so we need to read in every site (i.e. even the non-brownfield ones)

NPMS_all <- read.csv("NPMS/Analysis/NPMSsquares_loc.csv")

summary(NPMS_all)

# we can label target brownfield squares in this file before moving on
NPMS_target <- NPMS_sum[which(NPMS_sum$PERC > 5), ]

NPMS_all$TARGET <- ifelse(NPMS_all$NPMS_square %in% NPMS_target$GRIDREF, T, F)


# now, for every square in this list, we need to extract a few items into a new dataframe

# create a list to loop over and seed an output frame

NPMS_all_squares <- levels(droplevels(NPMS_all$NPMS_square))
n <- 0

NPMS_landuse <- data.frame()

## and now begin the loop

for (x in NPMS_all_squares){
  
  # first track progress
  n <- n + 1
  if (n == round(n, -2)){
    print(n)
  }
  
  # next, pull out square details
  square <- NPMS_all[which(NPMS_all$NPMS_square == x), ]
  out <- square[which(square$CORNER == "SW"), ]
  
  # create a polygon from the square's coordinates
  poly <- list(Polygons(list(Polygon(square[,6:7])), ID = x)) 
  SP <- SpatialPolygons(poly)
  
  # extract the LCM class for each 25x25m cell in the square
  cells <- extract(lcm2015, SP)
  
  # get rid of cells with value 0 (i.e. the sea)
  cells[[1]] <- cells[[1]][!cells[[1]] == 0]
  
  out$LANDUSE <- getmode(cells[[1]])
  NPMS_landuse <- rbind(NPMS_landuse, out)
  
}


summary(NPMS_landuse)

# this took a while - do an out and in to preserve

write.csv(NPMS_landuse, "NPMS/Analysis/NPMS_landuse.csv", row.names = F)

NPMS_landuse <- read.csv("NPMS/Analysis/NPMS_landuse.csv")



## next we want to identify the non-brownfield neighbours to each brownfield site
# this time we only need to loop over the brownfield sites, and within each one, the non-brownfield sites

NPMS_target_squares <- levels(droplevels(NPMS_target$GRIDREF))
NPMS_nontarget_squares <- NPMS_all_squares[which(!(NPMS_all_squares %in% NPMS_target_squares))]

n <- 0

NPMS_neighbours <- data.frame()


# do it in a loop

for (x in NPMS_target_squares){
  n <- n + 1
  if (n == round(n, -1)){
    print(n)
  }
  
  square <- NPMS_landuse[which(NPMS_landuse$NPMS_square == x), ]
  
  pairdists <- data.frame()
  
  # now an internal loop calculates distance to every other square
  for (y in NPMS_nontarget_squares){
    neighbour <- NPMS_landuse[which(NPMS_landuse$NPMS_square == y), ]
    
    SQUARE <- x
    NEIGHBOUR <- y
    
    SQUARE_use <- square$LANDUSE[1]
    NEIGHBOUR_use <- neighbour$LANDUSE[1]
    
    eastdist <- square$REast - neighbour$REast
    northdist <- square$RNorth - neighbour$RNorth
    
    PAIR_DIST <- sqrt(eastdist^2 + northdist^2)
    
    back <- cbind(SQUARE, NEIGHBOUR, SQUARE_use, NEIGHBOUR_use, PAIR_DIST)
    pairdists <- rbind(pairdists, back)
  }
  
  pairdists$PAIR_DIST <- as.numeric(as.character(pairdists$PAIR_DIST))
  
  # split into paired and different land use
  
  pairdists_matched <- pairdists[which(as.character(pairdists$SQUARE_use) == as.character(pairdists$NEIGHBOUR_use)), ]
  pairdists_other <- pairdists[which(as.character(pairdists$SQUARE_use) != as.character(pairdists$NEIGHBOUR_use)), ]
  
  # now take the top row of the matched frame and the top 3 of the other one
  
  matched_neighbour <- top_n(pairdists_matched, -1, pairdists_matched$PAIR_DIST)
  unmatched_neighbours <- top_n(pairdists_other, -3, pairdists_other$PAIR_DIST)
  
  neighbours <- rbind(matched_neighbour, unmatched_neighbours)
  
  NPMS_neighbours <- rbind(NPMS_neighbours, neighbours)
  
}

# this took a while - do an out and in to preserve

write.csv(NPMS_neighbours, "NPMS/Analysis/NPMS_neighbours.csv", row.names = F)

NPMS_neighbours <- read.csv("NPMS/Analysis/NPMS_neighbours.csv")


## finally, generate a non-redundant list of sites of interest

NPMS_good$CLASS <- "BROWNFIELD"
NPMS_brownfields <- NPMS_good[,c(1,5)]
colnames(NPMS_brownfields)[1] <- "GRIDREF"

NPMS_neighbours$COUNT <- 1
NPMS_pairs <- ddply(NPMS_neighbours, .(NEIGHBOUR), summarise,
                   PAIRS = sum(COUNT))

NPMS_pairs$CLASS <- "PAIR"
NPMS_pairs <- NPMS_pairs[,c(1,3)]
colnames(NPMS_pairs)[1] <- "GRIDREF"

NPMS_required <- rbind(NPMS_brownfields, NPMS_pairs)

# write this list out

write.csv(NPMS_required, "NPMS/Analysis/NPMS_requested.csv", row.names = F)



### WCBS

# now, we want to generate a small subset of this data for every site
# but to do so we need to read in every site (i.e. even the non-brownfield ones)

WCBS_all <- read.csv("WCBS/WCBSsquares_loc.csv")

summary(WCBS_all)

# we can label target brownfield squares in this file before moving on
WCBS_target <- WCBS_sum[which(WCBS_sum$PERC > 5), ]

WCBS_all$TARGET <- ifelse(WCBS_all$GRIDREF %in% WCBS_target$GRIDREF, T, F)


# now, for every square in this list, we need to extract a few items into a new dataframe


# create a list to loop over and seed an output frame

WCBS_all_squares <- levels(droplevels(WCBS_all$GRIDREF))
n <- 0

WCBS_landuse <- data.frame()

## and now begin the loop

for (x in WCBS_all_squares){
  
  # first track progress
  n <- n + 1
  if (n == round(n, -2)){
    print(n)
  }
  
  # next, pull out square details
  square <- WCBS_all[which(WCBS_all$GRIDREF == x), ]
  out <- square[which(square$CORNER == "SW"), ]
  
  # create a polygon from the square's coordinates
  poly <- list(Polygons(list(Polygon(square[,3:4])), ID = x)) 
  SP <- SpatialPolygons(poly)
  
  # extract the LCM class for each 25x25m cell in the square
  cells <- extract(lcm2015, SP)
  
  # get rid of cells with value 0 (i.e. the sea)
  cells[[1]] <- cells[[1]][!cells[[1]] == 0]
  
  out$LANDUSE <- getmode(cells[[1]])
  WCBS_landuse <- rbind(WCBS_landuse, out)
  
}


summary(WCBS_landuse)

# this took a while - do an out and in to preserve

write.csv(WCBS_landuse, "WCBS/WCBS_landuse.csv", row.names = F)

WCBS_landuse <- read.csv("WCBS/WCBS_landuse.csv")



## next we want to identify the non-brownfield neighbours to each brownfield site
# this time we only need to loop over the brownfield sites, and within each one, the non-brownfield sites

WCBS_target_squares <- levels(droplevels(WCBS_target$GRIDREF))
WCBS_nontarget_squares <- WCBS_all_squares[which(!(WCBS_all_squares %in% WCBS_target_squares))]

n <- 0

WCBS_neighbours <- data.frame()


# do it in a loop

for (x in WCBS_target_squares){
  n <- n + 1
  if (n == round(n, -1)){
    print(n)
  }
  
  square <- WCBS_landuse[which(WCBS_landuse$GRIDREF == x), ]
  
  pairdists <- data.frame()
  
  # now an internal loop calculates distance to every other square
  for (y in WCBS_nontarget_squares){
    neighbour <- WCBS_landuse[which(WCBS_landuse$GRIDREF == y), ]
    
    SQUARE <- x
    NEIGHBOUR <- y
    
    SQUARE_use <- square$LANDUSE[1]
    NEIGHBOUR_use <- neighbour$LANDUSE[1]
    
    eastdist <- square$EASTING - neighbour$EASTING
    northdist <- square$NORTHING - neighbour$NORTHING
    
    PAIR_DIST <- sqrt(eastdist^2 + northdist^2)
    
    back <- cbind(SQUARE, NEIGHBOUR, SQUARE_use, NEIGHBOUR_use, PAIR_DIST)
    pairdists <- rbind(pairdists, back)
  }
  
  pairdists$PAIR_DIST <- as.numeric(as.character(pairdists$PAIR_DIST))
  
  # split into paired and different land use
  
  pairdists_matched <- pairdists[which(as.character(pairdists$SQUARE_use) == as.character(pairdists$NEIGHBOUR_use)), ]
  pairdists_other <- pairdists[which(as.character(pairdists$SQUARE_use) != as.character(pairdists$NEIGHBOUR_use)), ]
  
  # now take the top row of the matched frame and the top 3 of the other one
  
  matched_neighbour <- top_n(pairdists_matched, -1, pairdists_matched$PAIR_DIST)
  unmatched_neighbours <- top_n(pairdists_other, -3, pairdists_other$PAIR_DIST)
  
  neighbours <- rbind(matched_neighbour, unmatched_neighbours)
  
  WCBS_neighbours <- rbind(WCBS_neighbours, neighbours)
  
}

# this took a while - do an out and in to preserve

write.csv(WCBS_neighbours, "WCBS/WCBS_neighbours.csv", row.names = F)

WCBS_neighbours <- read.csv("WCBS/WCBS_neighbours.csv")


## finally, generate a non-redundant list of sites of interest

WCBS_good$CLASS <- "BROWNFIELD"
WCBS_brownfields <- WCBS_good[,c(1,5)]
colnames(WCBS_brownfields)[1] <- "GRIDREF"

WCBS_neighbours$COUNT <- 1
WCBS_pairs <- ddply(WCBS_neighbours, .(NEIGHBOUR), summarise,
                   PAIRS = sum(COUNT))

WCBS_pairs$CLASS <- "PAIR"
WCBS_pairs <- WCBS_pairs[,c(1,3)]
colnames(WCBS_pairs)[1] <- "GRIDREF"

WCBS_required <- rbind(WCBS_brownfields, WCBS_pairs)

# write this list out

write.csv(WCBS_required, "WCBS/WCBS_requested.csv", row.names = F)




### finally check the overlap between squares, regardless of brownfield status

# first create a non-redundant list of all squares
all_squares_redundant <- append(BBS_all_squares,append(NPMS_all_squares,WCBS_all_squares))

# then loop over this for a non-redundant list
all_squares <- vector()
for (x in all_squares_redundant){
  if (!(x %in% all_squares)){
    all_squares <- append(all_squares,x)
  }
}

# now loop over this to identify coverage in the various datasets
coverage <- data.frame()

for (x in all_squares){
  SQUARE <- x
  
  BBS <- ifelse(x %in% BBS_all_squares, T, F)
  NPMS <- ifelse(x %in% NPMS_all_squares, T, F)
  WCBS <- ifelse(x %in% WCBS_all_squares, T, F)
  
  BIPARTITE <- ifelse(NPMS == T & WCBS == T, T, F)
  TRIPARTITE <- ifelse(BIPARTITE == T & BBS == T, T, F)
  
  out <- cbind(SQUARE, BBS, NPMS, WCBS, BIPARTITE, TRIPARTITE)
  coverage <- rbind(coverage, out)
}

summary(coverage)


