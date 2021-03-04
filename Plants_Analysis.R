##############################################################
####   Script for first attempt at analysing plant data   ####
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

source("../../CheckResidsFunction.R")
source("../../CheckConvergenceFunction.R")

### read in the data

raw_data <- read.csv("../Data/2015to2019/occurrences_2015to2019.csv", stringsAsFactors = T)

summary(raw_data)


# a key column here is "domin" but it's a bit messy atm
# let's tidy it up

good_domin <- c("1","2","3","4","5","6","7","8","9","10")

raw_data$domin_fixed <- ifelse(raw_data$domin %in% good_domin, as.character(raw_data$domin),
                               ifelse(raw_data$domin == "NULL", NA, 
                                      ifelse(raw_data$domin == "10. 91-100%", 10,
                                             substr(raw_data$domin,1,1))))

raw_data$domin_fixed <- as.numeric(as.character(raw_data$domin_fixed))

summary(raw_data$domin_fixed)



# clear out the NAs (there aren't very many and they aren't much use)

good_data <- raw_data[!is.na(raw_data$domin_fixed), ]

# and make domin_fixed into an ordinal factor

good_data$domin_fixed <- ordered(good_data$domin_fixed, levels = 1:10)

summary(good_data$domin_fixed)


### now we need to attach square data to these

locations <- read.csv("../Data/2015to2019/sampleInfoWithLatLong_2015to2019.csv", stringsAsFactors = T)

summary(locations)

# lists of squares and their pairs
pairs <- read.csv("NPMS_neighbours.csv", stringsAsFactors = T)
summary(pairs)

pairs$TYPE <- ifelse(pairs$SQUARE_use == pairs$NEIGHBOUR_use, "Matched","Neighbouring")


# make a list of squares we need to analyse first, and then a list of all squares (which we'll need for rarity analysis)

squares <- levels(droplevels(pairs$SQUARE))
neighbours <- levels(droplevels(pairs$NEIGHBOUR))

combined <- append(squares,neighbours)


# pick out the useful columns

locs_useful <- locations[,c(1:4,9:11)]


# merge this into the data by the matching column (which annoyingly is named different things in each)
colnames(locs_useful)[1] <- "sample_id"

data_locs <- merge(locs_useful, good_data)

summary(data_locs)


## now we can start to pick out actual stats on these squares
# let's start very simple:
# figure out species richness (raw and extrapolated)


all.sqs <- levels(droplevels(data_locs$monad))

richness <- data.frame()


for(x in all.sqs){
  sq.dat <- data_locs[which(data_locs$monad == x), ]
  print(x)
  
  # first crack raw spp richness
  sq.dat$COUNT <- 1
  
  sq.spp <- ddply(sq.dat, .(taxonversionkey, preferred_taxon), summarise,
                  RECS <- sum(COUNT))
  
  # now extrapolate spp richness across samples
  
  sq.pres <- sq.dat[,c(1,10,15)]
  
  # wrangle the data to be in wide form (one row per sample)
  
  sq.mat <- dcast(sq.pres, sample_id ~ preferred_taxon,
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
data_locs_mat <- data_locs[,c(7,10,14)]
data_locs_mat$presence <- 1


data_mat <- dcast(data_locs_mat, monad ~ preferred_taxon,
                value.var = "presence",
                fun.aggregate = sum)

# make rownames from col 1

rownames(data_mat) <- data_mat[,1]

data_mat <- data_mat[,-1]

# finally transpose the matrix so that assemblages are in columns and species in rows

data_mat_t <- t(data_mat)


### and a list of species with number of square-level occurrences

# first narrow down to one row per species*square
spp_occ_sq <- ddply(data_locs_mat, .(monad, preferred_taxon), summarise,
                    presence = mean(presence))

summary(spp_occ_sq)

# then summarise the data for each species
spp_occ <- ddply(spp_occ_sq, .(preferred_taxon), summarise,
                 squares = sum(presence))

# need to rejig the structure of this slightly

spp_occ_mat <- spp_occ$squares
names(spp_occ_mat) <- spp_occ$preferred_taxon
head(spp_occ_mat)

### then need to apply the functions rWeights and Irr

rarity.weights <- rWeights(spp_occ_mat)

square_rarity <- data.frame(Irr(assemblages = data_mat_t, W = rarity.weights))
hist(square_rarity$Irr)

# this actually causes some problems further down the line, because this index generates a lot of species with weight zero
# so I'm going to try and generate my own simpler version for comparison...
# there are 996 squares in total, so each species can be observed in 1-996 squares
# for my index I'll take 1-(obs.sqs/996)

rarity.weights$CJM <- 1-(rarity.weights$Q/996)

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

NPMS_landfill <- read.csv("NPMS_landfill_intersect_cleaned.csv")
summary(NPMS_landfill)


# we have two versions of the end-date for each site in this frame - 
# "licence surrendered" and "last input"
# the sites variously have neither, one or both of these recorded

# for simplicity, we want to condense these into a single column, taking "last input" if available,
# and "licence surrendered" otherwise

NPMS_landfill$EndDate <- ifelse(NPMS_landfill$lastinput != "", NPMS_landfill$lastinput, NPMS_landfill$lic_surren)

# now clean this new variable up and turn it into an actual date in R

NPMS_landfill$EndDate[which(NPMS_landfill$EndDate == "")] <- NA

NPMS_landfill$EndDate <- as.Date(NPMS_landfill$EndDate, format = "%d/%m/%Y")

# now we're ready to process further
# for example, we can see in the summary that square SU5863 has some 10 separate bits of landfill in it!

NPMS_sum <- ddply(NPMS_landfill, .(NPMS_sq), summarise,
                  LANDFILL_AREA = sum(Intersection_area),
                  EndDate = max(EndDate, na.rm = T))

summary(NPMS_sum)

# and now also add a numerical "days since end date" variable
NPMS_sum$TimeSinceClosure <- as.numeric(as.Date("31/12/2019", format = "%d/%m/%Y") - NPMS_sum$EndDate)
NPMS_sum$TimeSinceClosure[which(NPMS_sum$TimeSinceClosure == "Inf")] <- NA

# and due to scaling issues, make this into years (approximately)
NPMS_sum$TSCY <- NPMS_sum$TimeSinceClosure/365.25


# now we want to calculate the percentage of each square which is brownfield
# the area of each square is 1,000,000 m2

NPMS_sum$PERC <- NPMS_sum$LANDFILL_AREA*100/1000000
colnames(NPMS_sum)[1] <- "Square"

ss_landfill <- merge(NPMS_sum, ss_label, all.y = T)



# write this out for easy figure-plotting

write.csv(ss_landfill, "../../PreppedData/all_plants.csv")


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
# appears that target sites have lower SR than either other category... but...


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

# this suggests that target sites have HIGHER SR than either other category!

## that's a bit strange
# raw species richness suggests there's a significant effect in one direction,
# and estimated species richness in the other

# this could indicate there's a systematic bias in sampling completeness 
# (i.e. brownfield squares tend to be less well sampled) -
# so let's check that out

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

## this confirms it - sampling completeness is significantly lower in brownfield sites
# therefore we should only trust the estimated SR, not the raw numbers


summary(mod3e)



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


# unlike last time, this appears to be fine
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

chkres(mod5sr)

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






################################


### now we want to do something more complex - 
# identify which observations within brownfield squares come from brownfield sites
# and compare those to other observations within the same squares


## three steps -
# 1. identify all locations in the "landfill" squares
# 2. tag the ones that are in the brownfield sites themselves
# 3. compare!

# first read in all sites:

survey_locs <- read.csv("NPMSsites_loc.csv")

# now read in the sites which intersect with the brownfield sites themselves

survey_landfill <- read.csv("NPMSsites_landfill.csv")

landfill_locs <- levels(droplevels(as.factor(survey_landfill$location_id)))

# now mark each row in survey_locs_bf by whether or not the location is actually on a brownfield or not

survey_locs$landfill <- ifelse(survey_locs$location_id %in% landfill_locs, T, F)


## n.b. some of these are in squares *not* previously treated as brownfield (i.e. < 5% brownfield coverage)
# but we still want to include them here

# pick out all squares with a brownfield survey in them...

survey_locs_bf <- survey_locs[which(survey_locs$landfill == T), ]

landfill_squares <- levels(droplevels(survey_locs_bf$NPMS_square))


# identify all locations within these squares

data_intra <- data_locs[which(data_locs$monad %in% landfill_squares), ]



# and figure out richness for each *location* in each square

# figure out species richness (raw and extrapolated)

all.locs <- levels(droplevels(data_intra$location_id))

richness_intra <- data.frame()


for(x in all.locs){
  loc.dat <- data_locs[which(data_locs$location_id == x), ]
  print(x)
  
  # first crack raw spp richness
  loc.dat$COUNT <- 1
  
  loc.spp <- ddply(loc.dat, .(taxonversionkey, preferred_taxon), summarise,
                  RECS <- sum(COUNT))
  
  # now extrapolate spp richness across samples
  
  loc.pres <- loc.dat[,c(1,10,15)]
  
  # wrangle the data to be in wide form (one row per sample)
  
  loc.mat <- dcast(loc.pres, sample_id ~ preferred_taxon,
                  value.var = "COUNT",
                  fun.aggregate = sum)
  
  # get rid of the sample ID column, otherwise it gets treated as an extra species!
  loc.mat <- loc.mat[,-1]
  
  out <- specpool(loc.mat)
  out$Location <- x
  
  richness_intra <- rbind(richness_intra, out)
  
}

summary(richness_intra)


## now do the same for rarity


### need a full location-level matrix (presence/absence, one row per location)

# trim out only the columns we need
data_locs_mat_intra <- data_locs[,c(4,10,14)]
data_locs_mat_intra$presence <- 1


data_mat_intra <- dcast(data_locs_mat_intra, location_id ~ preferred_taxon,
                  value.var = "presence",
                  fun.aggregate = sum)

# make rownames from col 1

rownames(data_mat_intra) <- data_mat_intra[,1]

data_mat_intra <- data_mat_intra[,-1]

# finally transpose the matrix so that assemblages are in columns and species in rows

data_mat_intra_t <- t(data_mat_intra)


### and a list of species with number of square-level occurrences

# first narrow down to one row per species*location
spp_occ_sq_intra <- ddply(data_locs_mat_intra, .(location_id, preferred_taxon), summarise,
                    presence = mean(presence))

summary(spp_occ_sq_intra)

# then summarise the data for each species
spp_occ_intra <- ddply(spp_occ_sq_intra, .(preferred_taxon), summarise,
                 squares = sum(presence))

# need to rejig the structure of this slightly

spp_occ_mat_intra <- spp_occ_intra$squares
names(spp_occ_mat_intra) <- spp_occ_intra$preferred_taxon
head(spp_occ_mat_intra)

### then need to apply the functions rWeights and Irr

rarity.weights_intra <- rWeights(spp_occ_mat_intra)

square_rarity_intra <- data.frame(Irr(assemblages = data_mat_intra_t, W = rarity.weights_intra))
hist(square_rarity$Irr)

# this actually causes some problems further down the line, because this index generates a lot of species with weight zero
# so I'm going to try and generate my own simpler version for comparison...
# there are 3899 locations in total, so each species can be observed in 1-3899 squares
# for my index I'll take 1-(obs.sqs/3899)

rarity.weights_intra$CJM <- 1-(rarity.weights_intra$Q/3899)

# normalise these between 0 and 1 (slightly boosting the rare species and taking away from the common)
minW_intra <- min(rarity.weights_intra$CJM)
maxW_intra <- max(rarity.weights_intra$CJM)

rarity.weights_intra$CJMz <- (rarity.weights_intra$CJM - minW_intra)/(maxW_intra - minW_intra)

summary(rarity.weights_intra$CJMz)

rarity.weights.cjm_intra <- rarity.weights_intra[,c(1,2,6,4)]
colnames(rarity.weights.cjm_intra)[3] <- "W"

square_rarity.cjm_intra <- data.frame(Irr(assemblages = data_mat_intra_t, W = rarity.weights.cjm_intra))
hist(square_rarity.cjm_intra$Irr)


### now merge all these indices back together so that we have one row per square in a single dataframe
square_rarity_intra$Location <- rownames(square_rarity_intra)
square_rarity.cjm_intra$Location <- rownames(square_rarity.cjm_intra)

colnames(square_rarity.cjm_intra) <- c("Irr.cjm", "Richness.cjm", "Location")

square_rr_intra <- merge(square_rarity_intra, square_rarity.cjm_intra)

square_stats_intra <- merge(richness_intra, square_rr_intra)


### now we're ready to run some tests



### now we are ready to start some basic analyses
# first pick out the stats for the squares

# first, label whether each loc is on landfill or not

square_stats_intra$landfill <- ifelse(square_stats_intra$Location %in% landfill_locs, "Landfill","Neighbour")


# re-introduce the squares data for a grouping factor

squares_intra <- ddply(data_intra, .(monad, location_id), summarise,
                       arbit = mean(LATITUDE))

colnames(squares_intra) <- c("Square", "Location","Arbit")

# merge

ss_intra <- merge(squares_intra,square_stats_intra)

# calculate sampling completeness

ss_intra$Completeness <- ss_intra$Species*100/ss_intra$chao

summary(ss_intra)


### we are finally ready to begin interrogating the data!

# first a very simple comparison

ss_intra$landfill <- relevel(as.factor(ss_intra$landfill), ref = "Landfill")


mod6r <- glmer(Species ~ landfill + (1|Square),
               family = poisson (link = "log"),
               data = ss_intra)

summary(mod6r)
drop1(mod6r, test = "Chi")

chkres(mod6r)

plot(ss_intra$Species ~ ss_intra$landfill)

# appears that target sites have lower SR than either other category... but...


## repeat with estimated richness

mod6e <- glmer(round(chao) ~ landfill + (1|Square),
               family = poisson (link = "log"),
               data = ss_intra)

summary(mod6e)
drop1(mod6e, test = "Chi")

chkres(mod6e)

plot(ss_intra$chao ~ ss_intra$landfill)


## that's a bit strange
# raw species richness suggests there's a significant effect in one direction,
# and estimated species richness no effect

# again this could indicate there's a systematic bias in sampling completeness 
# (i.e. brownfield sites tend to be less well sampled) -
# so let's check that out

hist(ss_intra$Completeness)

mod6c <- lmer(Completeness ~ landfill + (1|Square),
              data = ss_intra)

summary(mod6c)
drop1(mod6c, test = "Chi")

chkres(mod6c)

plot(ss_intra$Completeness ~ ss_intra$landfill)

## in this instance it looks ok - so we could probably trust either

summary(mod6e)




## repeat for rarity


mod6i <- glmer(Irr.cjm ~ landfill + (1|Square),
               family = binomial (link = "logit"),
               weights = Species,
               data=ss_intra)

summary(mod6i)
drop1(mod6i, test = "Chi")

chkres(mod6i)

plot(ss_intra$Irr.cjm ~ ss_intra$landfill)












