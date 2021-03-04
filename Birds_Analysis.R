##############################################################
####   Script for first attempt at analysing bird data   ####
##############################################################

##### script setup ####

### Clear the current workspace (DOESN'T CLEAR LOADED LIBRARIES)
rm(list=ls())

### install if necessary and then load the libraries you need

j <- c("rstudioapi","plyr","raster","dplyr","vegan","reshape2","Rarity","lme4","MASS","glmmTMB")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

### set working directory to script's saved location

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## read in some functions we might need

source("../CheckResidsFunction.R")
source("../CheckConvergenceFunction.R")

### read in the data

raw_data <- read.csv("Data/Brownfield BBS Data_All years_MACGREGOR.csv")

summary(raw_data)


# a key column here is "TOT" but it's a bit messy atm
# let's tidy it up

raw_data$TOT <- as.numeric(as.character(raw_data$TOT))

summary(raw_data$TOT)



# clear out the NAs (there aren't very many and they aren't much use), and also a few zeroes and (weirdly) negatives

good_data <- raw_data[which(raw_data$TOT > 0), ]

summary(good_data)



### now we need to attach square data to these

# lists of squares and their pairs
pairs <- read.csv("BBS_neighbours.csv", stringsAsFactors = T)
summary(pairs)

pairs$TYPE <- ifelse(pairs$SQUARE_use == pairs$NEIGHBOUR_use, "Matched","Neighbouring")


# make a list of squares we need to analyse first, and then a list of all squares (which we'll need for rarity analysis)

squares <- levels(droplevels(pairs$SQUARE))
neighbours <- levels(droplevels(pairs$NEIGHBOUR))

combined <- append(squares,neighbours)


# read in some location data

locations  <- read.csv("BBS_landuse.csv", stringsAsFactors = T)


# pick out the useful columns

locs_useful <- locations[,c(1,5:6,8:9,11)]


# merge this into the data by the matching column (which annoyingly is named different things in each)
data_locs <- merge(locs_useful, good_data)

summary(data_locs)


## now we can start to pick out actual stats on these squares
# let's start very simple:
# figure out species richness (raw and extrapolated)


all.sqs <- levels(droplevels(data_locs$GRIDREF))

richness <- data.frame()


for(x in all.sqs){
  sq.dat <- data_locs[which(data_locs$GRIDREF == x), ]
  print(x)
  
  # first crack raw spp richness
  sq.dat$COUNT <- 1
  
  sq.spp <- ddply(sq.dat, .(ENGLISH_NAME, CBC_CODE), summarise,
                  RECS <- sum(COUNT),
                  TOT <- sum(TOT),
                  MEAN <- mean(TOT))
  
  # now extrapolate spp richness across samples
  
  sq.pres <- sq.dat[,c(10,13,15)]
  
  # wrangle the data to be in wide form (one row per sample)
  
  sq.mat <- dcast(sq.pres, OBS_DT ~ ENGLISH_NAME,
                   value.var = "COUNT",
                   fun.aggregate = sum)
  
  # get rid of the sample ID column, otherwise it gets treated as an extra species!
  sq.mat <- sq.mat[,-1]
  
  out <- specpool(sq.mat)
  out$Square <- x
  
  richness <- rbind(richness, out)

}

summary(richness)



#### figure out species rarity 
# use R package "Rarity"


### need a full site-level matrix (presence/absence, one row per square)

# trim out only the columns we need
data_locs_mat <- data_locs[,c(1,10,14)]
data_locs_mat$presence <- 1


data_mat <- dcast(data_locs_mat, GRIDREF ~ ENGLISH_NAME,
                value.var = "presence",
                fun.aggregate = sum)

# make rownames from col 1

rownames(data_mat) <- data_mat[,1]

data_mat <- data_mat[,-1]

# finally transpose the matrix so that assemblages are in columns and species in rows

data_mat_t <- t(data_mat)


### and a list of species with number of square-level occurrences

# first narrow down to one row per species*square
spp_occ_sq <- ddply(data_locs_mat, .(GRIDREF, ENGLISH_NAME), summarise,
                    presence = mean(presence))

summary(spp_occ_sq)

# then summarise the data for each species
spp_occ <- ddply(spp_occ_sq, .(ENGLISH_NAME), summarise,
                 squares = sum(presence))

# need to rejig the structure of this slightly

spp_occ_mat <- spp_occ$squares
names(spp_occ_mat) <- spp_occ$ENGLISH_NAME
head(spp_occ_mat)

### then need to apply the functions rWeights and Irr

rarity.weights <- rWeights(spp_occ_mat)

square_rarity <- data.frame(Irr(assemblages = data_mat_t, W = rarity.weights))
hist(square_rarity$Irr)

# this actually causes some problems further down the line, because this index generates a lot of species with weight zero
# so I'm going to try and generate my own simpler version for comparison...
# there are 708 squares in total, so each species can be observed in 1-708 squares
# for my index I'll take 1-(obs.sqs/708)

rarity.weights$CJM <- 1-(rarity.weights$Q/708)

# normalise these between 0 and 1 (slightly boosting the rare species and taking away from the common)
minW <- min(rarity.weights$CJM)
maxW <- max(rarity.weights$CJM)

rarity.weights$CJMz <- (rarity.weights$CJM - minW)/(maxW - minW)

summary(rarity.weights$CJMz)

rarity.weights.cjm <- rarity.weights[,c(1,2,6,4)]
colnames(rarity.weights.cjm)[3] <- "W"

square_rarity.cjm <- data.frame(Irr(assemblages = data_mat_t, W = rarity.weights.cjm))
hist(square_rarity.cjm$Irr)


### now merge all these indices back together so that we have one row per square in a single dataframe
square_rarity$Square <- rownames(square_rarity)
square_rarity.cjm$Square <- rownames(square_rarity.cjm)

colnames(square_rarity.cjm) <- c("Irr.cjm", "Richness.cjm", "Square")

square_rr <- merge(square_rarity, square_rarity.cjm)

square_stats <- merge(richness, square_rr)


### now we are ready to start some basic analyses
# first pick out the stats for the squares in pairs

# several ways to tackle this
# I'll want a long-form version (still one row for every square, but labelling which group they belong to)
# and possibly a wide-form version (one row per brownfield square, columns for )

# first, add a factor label to each square

ss_label <- data.frame()

y <- 0

for (x in squares){
  print(x)
  y <- y+1
  
  sq.pair <- pairs[which(pairs$SQUARE == x), ]
  
  sq.sqs <- append(x, levels(droplevels(sq.pair$NEIGHBOUR)))
  
  sq.dat <- square_stats[which(square_stats$Square %in% sq.sqs), ]
  
  
  sq.pair.base <- sq.pair[,c('NEIGHBOUR','TYPE')]
  sq.pair.base$NEIGHBOUR <- as.character(sq.pair.base$NEIGHBOUR)
  sq.pair.base[5,] <- c(x, "Target")
  colnames(sq.pair.base) <- c("Square", "Type")
  
  
  sq.dat$Group <- y
  
  sq.type <- merge(sq.pair.base, sq.dat)
  
  
  ss_label <- rbind(ss_label, sq.type)
  
  
}

summary(ss_label)

ss_label$Type <- as.factor(ss_label$Type)
ss_label$Group <- as.factor(ss_label$Group)

ss_label$Completeness <- ss_label$Species/ss_label$chao

summary(ss_label)


# tack in the area of brownfield site

BBS_landfill <- read.csv("BBS_landfill_intersect_cleaned.csv")
summary(BBS_landfill)

# we have two versions of the end-date for each site in this frame - 
# "licence surrendered" and "last input"
# the sites variously have neither, one or both of these recorded

# for simplicity, we want to condense these into a single column, taking "last input" if available,
# and "licence surrendered" otherwise

BBS_landfill$EndDate <- ifelse(BBS_landfill$lastinput != "", BBS_landfill$lastinput, BBS_landfill$lic_surren)

# now clean this new variable up and turn it into an actual date in R

BBS_landfill$EndDate[which(BBS_landfill$EndDate == "")] <- NA

BBS_landfill$EndDate <- as.Date(BBS_landfill$EndDate, format = "%d/%m/%Y")

# now we're ready to process further
# for example, we can see in the summary that square SU5863 has some 10 separate bits of landfill in it!

BBS_sum <- ddply(BBS_landfill, .(GRIDREF), summarise,
                 LANDFILL_AREA = sum(Intersection_area),
                 EndDate = max(EndDate, na.rm = T))

summary(BBS_sum)

# and now also add a numerical "days since end date" variable
BBS_sum$TimeSinceClosure <- as.numeric(as.Date("31/12/2019", format = "%d/%m/%Y") - BBS_sum$EndDate)
BBS_sum$TimeSinceClosure[which(BBS_sum$TimeSinceClosure == "Inf")] <- NA

# and due to scaling issues, make this into years (approximately)
BBS_sum$TSCY <- BBS_sum$TimeSinceClosure/365.25



# now we want to calculate the percentage of each square which is brownfield
# the area of each square is 1,000,000 m2

BBS_sum$PERC <- BBS_sum$LANDFILL_AREA*100/1000000
colnames(BBS_sum)[1] <- "Square"

ss_landfill <- merge(BBS_sum, ss_label, all.y = T)


# write this out for easy figure-plotting

write.csv(ss_landfill, "../PreppedData/all_birds.csv")


### we are finally ready to begin interrogating the data!

# first a very simple comparison

ss_landfill$Type <- relevel(ss_landfill$Type, ref = "Target")

matched.dat <- ss_landfill[which(ss_landfill$Type != "Neighbouring"), ]
grouped.dat <- ss_landfill[which(ss_landfill$Type != "Matched"), ]


# matched

hist(matched.dat$Species)

mod1r <- glmer(Species ~ Type + (1|Group),
              family = poisson (link = "log"),
              data = matched.dat)

summary(mod1r)
drop1(mod1r, test = "Chi")

chkres(mod1r, matched.dat$Type)


# grouped 

mod2r <- glmer(Species ~ Type + (1|Group),
            family = poisson (link = "log"),
            data = grouped.dat)

summary(mod2r)
drop1(mod2r, test = "Chi")

chkres(mod2r, grouped.dat$Type)


# combined

mod3r <- glmer(Species ~ Type + (1|Group),
               family = poisson (link = "log"),
               data = ss_landfill)

summary(mod3r)
drop1(mod3r, test = "Chi")

chkres(mod3r, ss_landfill$Type)

plot(ss_landfill$Species ~ ss_landfill$Type)
# appears that target sites have higher SR than either other category...


## repeat with estimated richness

# matched

hist(matched.dat$chao)

mod1e <- glmer(round(chao) ~ Type + (1|Group),
              family = poisson (link = "log"),
              data = matched.dat)

summary(mod1e)
drop1(mod1e, test = "Chi")

chkres(mod1e, matched.dat$Type)


# grouped 

mod2e <- glmer(round(chao) ~ Type + (1|Group),
              family = poisson (link = "log"),
              data = grouped.dat)

summary(mod2e)
drop1(mod2e, test = "Chi")

chkres(mod2e, grouped.dat$Type)


# combined

mod3e <- glmer(round(chao) ~ Type + (1|Group),
               family = poisson (link = "log"),
               data = ss_landfill)

summary(mod3e)
drop1(mod3e, test = "Chi")

chkres(mod3e, ss_landfill$Type)

plot(ss_landfill$chao ~ ss_landfill$Type)

# this also suggests that target sites have HIGHER SR than either other category!

# no immediate sign of the bias here that was present in the plant data, but
# let's check sampling completeness out anyway

hist(ss_landfill$Completeness)

# matched

mod1c <- glmer(Completeness ~ Type + (1|Group),
               family = binomial (link = "logit"),
               weights = n,
              data = matched.dat)

summary(mod1c)
drop1(mod1c, test = "Chi")

chkres(mod1c, matched.dat$Type)



# grouped

mod2c <- glmer(Completeness ~ Type + (1|Group),
               family = binomial (link = "logit"),
               weights = n,
              data = grouped.dat)

summary(mod2c)
drop1(mod2c, test = "Chi")

chkres(mod2c, grouped.dat$Type)



# combined

mod3c <- glmer(Completeness ~ Type + (1|Group),
               family = binomial (link = "logit"),
               weights = n,
              data = ss_landfill)

summary(mod3c)
drop1(mod3c, test = "Chi")

chkres(mod3c)

plot(ss_landfill$Completeness ~ ss_landfill$Type)

# no cause for concern here



## species richness is significantly higher in brownfield-type habitats than other habitat types, 
# and also significantly higher in brownfield squares than non-brownfields of same type (though with a smaller effect)

# in other words, the habitats associated with brownfield sites have higher species richness than other local habitat types
# AND brownfield sites have higher SR than other sites with similar habitat...?


# repeat for Irr

# matched

hist(matched.dat$Irr.cjm)

mod1i <- glmer(Irr.cjm ~ Type + (1|Group),
            family = binomial (link = "logit"),
            weights = n,
            data = matched.dat)

summary(mod1i)
drop1(mod1i, test = "Chi")

chkres(mod1i, matched.dat$Type)


# grouped

mod2i <- glmer(Irr.cjm ~ Type + (1|Group),
              family = binomial (link = "logit"),
              weights = n,
              data=grouped.dat)

summary(mod2i)
drop1(mod2i, test = "Chi")

chkres(mod2i)


# combined

mod3i <- glmer(Irr.cjm ~ Type + (1|Group),
               family = binomial (link = "logit"),
               weights = n,
               data=ss_landfill)

summary(mod3i)
drop1(mod3i, test = "Chi")

chkres(mod3i)

plot(ss_landfill$Irr.cjm ~ ss_landfill$Type)


### AND brownfield sites have higher rarity scores than non-brownfields of either type!
# (this result was absent from plants)



## the next set of questions relates to the effect of *area of brownfield* within target sites


ss_brownfield <- ss_landfill[which(ss_landfill$PERC > 5), ]


# first, estimated SR

mod4sr <- glm(round(chao) ~ LANDFILL_AREA,
                 family = poisson (link = "log"),
                 data = ss_brownfield)

summary(mod4sr)
drop1(mod4sr, test = "Chi")

chkres(mod4sr, ss_brownfield$LANDFILL_AREA, ss_brownfield$Type)


plot(log(ss_brownfield$chao) ~ log(ss_brownfield$LANDFILL_AREA))


# next, completeness (just to check)

mod4c <-  glm(Completeness ~ LANDFILL_AREA,
              data = ss_brownfield)

summary(mod4c)
drop1(mod4c, test = "Chi")

chkres(mod4c, ss_brownfield$LANDFILL_AREA, ss_brownfield$Type)

plot(ss_brownfield$Completeness ~ log(ss_brownfield$LANDFILL_AREA))



# unlike last time, this appears to have a slight problem - sampling completeness is higher where brownfield sites are smaller
# however, the effect size is minute
# so, moving on to Irr

mod4i <-  glm(Irr.cjm ~ LANDFILL_AREA,
              family = binomial (link = "logit"),
              data = ss_brownfield)


summary(mod4i)
drop1(mod4i, test = "Chi")

chkres(mod4i)


plot(ss_brownfield$Irr.cjm ~ log(ss_brownfield$LANDFILL_AREA))


## for completeness, do these analyses with the percentage landfill by area (instead of raw area)

# first, estimated SR

mod5sr <- glm(round(chao) ~ log(PERC),
              family = poisson (link = "log"),
              data = ss_brownfield)

summary(mod5sr)
drop1(mod5sr, test = "Chi")

chkres(mod5sr)

plot(log(ss_brownfield$chao) ~ log(ss_brownfield$PERC))


# and observed

mod5so <- glm(Species ~ log(PERC),
              family = poisson (link = "log"),
              data = ss_brownfield)

summary(mod5so)
drop1(mod5so, test = "Chi")

chkres(mod5so)

plot(log(ss_brownfield$Species) ~ log(ss_brownfield$PERC))




# next, completeness (just to check)

mod5c <-  glm(Completeness ~ PERC,
              data = ss_brownfield)

summary(mod5c)
drop1(mod5c, test = "Chi")

chkres(mod5c)

plot(ss_brownfield$Completeness ~ log(ss_brownfield$PERC))


# unlike last time, this appears to be fine
# so, moving on to Irr

mod5i <-  glm(Irr.cjm ~ PERC,
              family = binomial (link = "logit"),
              data = ss_brownfield)


summary(mod5i)
drop1(mod5i, test = "Chi")

chkres(mod5i)


plot(ss_brownfield$Irr.cjm ~ log(ss_brownfield$PERC))



## and now we have a set of questions relating to the age of the brownfield sites

ss_brownfield_age <- ss_brownfield[which(!is.na(ss_brownfield$TimeSinceClosure)), ]


# first, estimated SR

mod6sr <- glm(round(chao) ~ TSCY,
              family = poisson (link = "log"),
              data = ss_brownfield_age)

summary(mod6sr)
drop1(mod6sr, test = "Chi")

chkres(mod6sr)


plot(log(ss_brownfield$chao) ~ ss_brownfield$TSCY)

# and observed

mod6so <- glm(Species ~ TSCY,
              family = poisson (link = "log"),
              data = ss_brownfield_age)

summary(mod6so)
drop1(mod6so, test = "Chi")

chkres(mod6so)

plot(log(ss_brownfield$Species) ~ ss_brownfield$TSCY)


# next, completeness (just to check)

mod6c <-  glm(Completeness ~ TSCY,
              data = ss_brownfield_age)

summary(mod6c)
drop1(mod6c, test = "Chi")

chkres(mod6c)

plot(ss_brownfield$Completeness ~ ss_brownfield$TSCY)



# no problem here
# so, moving on to Irr

mod6i <-  glm(Irr.cjm ~ TSCY,
              family = binomial (link = "logit"),
              data = ss_brownfield_age)


summary(mod6i)
drop1(mod6i, test = "Chi")

chkres(mod6i)


plot(ss_brownfield$Irr.cjm ~ ss_brownfield$TimeSinceClosure)









### this script stops here -
# it's not possible to do any within-square analyses with the BBS data

