##############################################################
####   Script for first attempt at analysing insect data   ####
##############################################################

##### script setup ####

### Clear the current workspace (DOESN'T CLEAR LOADED LIBRARIES)
rm(list=ls())

### install if necessary and then load the libraries you need

j <- c("rstudioapi","plyr","raster","dplyr","vegan","reshape2","Rarity","lme4","MASS","glmmTMB","nlme","car")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

### set working directory to script's saved location

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## read in some functions we might need

source("../../CheckResidsFunction.R")
source("../../CheckConvergenceFunction.R")

### print session info

sessionInfo()


### read in the data

raw_data <- read.csv("../REFERENCE RECORDS 2006-2019.csv")

summary(raw_data)


# a key column here is "COUNT" but it's a bit messy atm
# let's tidy it up

# clear out the NAs (there aren't very many and they aren't much use), and also a few zeroes and (weirdly) negatives

good_data <- raw_data[which(raw_data$COUNT > 0), ]

summary(good_data)


### unlike the other data sets the species data is held separately - read it in:

species <- read.csv("../SPECIES UPDATED.csv")

summary(species)

## as well as butterflies, this dataset contains a lot of other insects
# and in addition, the surveys were focussed on "common butterflies" (so "other butterflies" may not have been a target for recording)

# therefore it will be useful to have a couple of factors to divide things up later

species$BUTTERFLY <- ifelse(!is.na(species$COMMON_BUTTERFLY), T, 
                            ifelse(!is.na(species$OTHER_BUTTERFLY), T, F))

species$GROUP <- ifelse(!is.na(species$OTHER_BUTTERFLY), "Other butterfly",
                        ifelse(!is.na(species$COMMON_BUTTERFLY), "Common butterfly",
                               ifelse(!is.na(species$MOTH), "Moth",
                                      ifelse(!is.na(species$DRAGONFLY), "Dragonfly", 
                                             ifelse(!is.na(species$DAMSELFLY), "Damselfly", "Other")))))

species <- species[c(1,2,8,9)]

good_data_sp <- merge(species, good_data)


### now we need to attach square data to these

# lists of squares and their pairs
pairs <- read.csv("../WCBS_neighbours.csv", stringsAsFactors = T)
summary(pairs)

pairs$TYPE <- ifelse(pairs$SQUARE_use == pairs$NEIGHBOUR_use, "Matched","Neighbouring")


# make a list of squares we need to analyse first, and then a list of all squares (which we'll need for rarity analysis)

squares <- levels(droplevels(pairs$SQUARE))
neighbours <- levels(droplevels(pairs$NEIGHBOUR))

combined <- append(squares,neighbours)


# read in some location data

locations  <- read.csv("../REFERENCE VISITS 2006-2019.csv", stringsAsFactors = T)


# pick out the useful columns

locs_useful <- locations[,c(1:4)]


# merge this into the data by the matching column (which annoyingly is named different things in each)
data_locs <- merge(locs_useful, good_data_sp)

summary(data_locs)

all.sqs <- ddply(data_locs, .(GRID_REF), summarise,
                 TOTAL = sum(COUNT))

colnames(all.sqs)[1] <- "Square"

summary(all.sqs)


#### BUTTERFLIES

# start the analyses from this point on the butterflies only, since they're the focus of the scheme
# we'll do some replicated analyses on other subsets later but let's start with this one

data_locs_but <- data_locs[which(data_locs$BUTTERFLY==T), ]


## now we can start to pick out actual stats on these squares
# let's start very simple:
# figure out species richness (raw and extrapolated)


all.sqs.but <- levels(droplevels(data_locs_but$GRID_REF))

richness_but <- data.frame()


for(x in all.sqs.but){
  sq.dat <- data_locs_but[which(data_locs_but$GRID_REF == x), ]
  print(x)
  
  # first crack raw spp richness
  sq.dat$PRESENT <- 1
  
  sq.spp <- ddply(sq.dat, .(COMMON_NAME, BMSCODE), summarise,
                  RECS = sum(PRESENT),
                  TOT = sum(COUNT),
                  MEAN = mean(COUNT))
  
  # now extrapolate spp richness across samples
  
  sq.pres <- sq.dat[,c(4,6,12)]
  
  # wrangle the data to be in wide form (one row per sample)
  
  sq.mat <- dcast(sq.pres, VISIT_DATE ~ COMMON_NAME,
                   value.var = "COUNT",
                   fun.aggregate = sum)
  
  # get rid of the sample ID column, otherwise it gets treated as an extra species!
  sq.mat <- sq.mat[,-1]
  
  out <- specpool(sq.mat)
  out$Square <- x
  
  richness_but <- rbind(richness_but, out)

}

summary(richness_but)



#### figure out species rarity 
# use R package "Rarity"


### need a full site-level matrix (presence/absence, one row per square)

# trim out only the columns we need
data_locs_mat_but <- data_locs_but[,c(2,6,12)]
data_locs_mat_but$presence <- 1


data_mat_but <- dcast(data_locs_mat_but, GRID_REF ~ COMMON_NAME,
                value.var = "presence",
                fun.aggregate = sum)

# make rownames from col 1

rownames(data_mat_but) <- data_mat_but[,1]

data_mat_but <- data_mat_but[,-1]

# finally transpose the matrix so that assemblages are in columns and species in rows

data_mat_but_t <- t(data_mat_but)


### and a list of species with number of square-level occurrences

# first narrow down to one row per species*square
spp_occ_sq_but <- ddply(data_locs_mat_but, .(GRID_REF, COMMON_NAME), summarise,
                    presence = mean(presence))

summary(spp_occ_sq_but)

# then summarise the data for each species
spp_occ_but <- ddply(spp_occ_sq_but, .(COMMON_NAME), summarise,
                 squares = sum(presence))

# need to rejig the structure of this slightly

spp_occ_mat_but <- spp_occ_but$squares
names(spp_occ_mat_but) <- spp_occ_but$COMMON_NAME
head(spp_occ_mat_but)

### then need to apply the functions rWeights and Irr

rarity.weights_but <- rWeights(spp_occ_mat_but)

square_rarity_but <- data.frame(Irr(assemblages = data_mat_but_t, W = rarity.weights_but))
hist(square_rarity_but$Irr)

# this actually causes some problems further down the line, because this index generates a lot of species with weight zero
# so I'm going to try and generate my own simpler version for comparison...
# there are 2096 squares in total, so each species can be observed in 1-2096 squares
# for my index I'll take 1-(obs.sqs/2096)

rarity.weights_but$CJM <- 1-(rarity.weights_but$Q/2096)

# normalise these between 0 and 1 (slightly boosting the rare species and taking away from the common)
minW <- min(rarity.weights_but$CJM)
maxW <- max(rarity.weights_but$CJM)

rarity.weights_but$CJMz <- (rarity.weights_but$CJM - minW)/(maxW - minW)

summary(rarity.weights_but$CJMz)

rarity.weights.cjm_but <- rarity.weights_but[,c(1,2,6,4)]
colnames(rarity.weights.cjm_but)[3] <- "W"

square_rarity.cjm_but <- data.frame(Irr(assemblages = data_mat_but_t, W = rarity.weights.cjm_but))
hist(square_rarity.cjm_but$Irr)


### now merge all these indices back together so that we have one row per square in a single dataframe
square_rarity_but$Square <- rownames(square_rarity_but)
square_rarity.cjm_but$Square <- rownames(square_rarity.cjm_but)

colnames(square_rarity.cjm_but) <- c("Irr.cjm", "Richness.cjm", "Square")

square_rr_but <- merge(square_rarity_but, square_rarity.cjm_but)

square_stats_but <- merge(richness_but, square_rr_but)

# restore squares with zero richness, and fill in blanks
square_stats_complete_but <- merge(all.sqs, square_stats_but, all = T)

square_stats_complete_but[is.na(square_stats_complete_but)] <- 0



### now we are ready to start some basic analyses
# first pick out the stats for the squares in pairs

# several ways to tackle this
# I'll want a long-form version (still one row for every square, but labelling which group they belong to)
# and possibly a wide-form version (one row per brownfield square, columns for )

# first, add a factor label to each square

ss_label_but <- data.frame()

y <- 0

for (x in squares){
  print(x)
  y <- y+1
  
  sq.pair <- pairs[which(pairs$SQUARE == x), ]
  
  sq.sqs <- append(x, levels(droplevels(sq.pair$NEIGHBOUR)))
  
  sq.dat <- square_stats_complete_but[which(square_stats_complete_but$Square %in% sq.sqs), ]
  
  
  sq.pair.base <- sq.pair[,c('NEIGHBOUR','TYPE')]
  sq.pair.base$NEIGHBOUR <- as.character(sq.pair.base$NEIGHBOUR)
  sq.pair.base[5,] <- c(x, "Target")
  colnames(sq.pair.base) <- c("Square", "Type")
  
  
  sq.dat$Group <- y
  
  sq.type <- merge(sq.pair.base, sq.dat)
  
  
  ss_label_but <- rbind(ss_label_but, sq.type)
  
  
}

summary(ss_label_but)

ss_label_but$Type <- as.factor(ss_label_but$Type)
ss_label_but$Group <- as.factor(ss_label_but$Group)

ss_label_but$Completeness <- ss_label_but$Species/ss_label_but$chao

summary(ss_label_but)


# tack in the area of brownfield site

WCBS_landfill <- read.csv("WCBS_landfill_intersect_cleaned.csv")
summary(WCBS_landfill)


# we have two versions of the end-date for each site in this frame - 
# "licence surrendered" and "last input"
# the sites variously have neither, one or both of these recorded

# for simplicity, we want to condense these into a single column, taking "last input" if available,
# and "licence surrendered" otherwise

WCBS_landfill$EndDate <- ifelse(WCBS_landfill$lastinput != "", WCBS_landfill$lastinput, WCBS_landfill$lic_surren)

# now clean this new variable up and turn it into an actual date in R

WCBS_landfill$EndDate[which(WCBS_landfill$EndDate == "")] <- NA

WCBS_landfill$EndDate <- as.Date(WCBS_landfill$EndDate, format = "%d/%m/%Y")

# now we're ready to process further
# for example, we can see in the summary that square SU3524 has some 10 separate bits of landfill in it!

WCBS_sum <- ddply(WCBS_landfill, .(WCBS_sq), summarise,
                  LANDFILL_AREA = sum(Intersection_area),
                  EndDate = max(EndDate, na.rm = T))

summary(WCBS_sum)

# and now also add a numerical "days since end date" variable
WCBS_sum$TimeSinceClosure <- as.numeric(as.Date("31/12/2019", format = "%d/%m/%Y") - WCBS_sum$EndDate)
WCBS_sum$TimeSinceClosure[which(WCBS_sum$TimeSinceClosure == "Inf")] <- NA

# and due to scaling issues, make this into years (approximately)
WCBS_sum$TSCY <- WCBS_sum$TimeSinceClosure/365.25




# now we want to calculate the percentage of each square which is brownfield
# the area of each square is 1,000,000 m2

WCBS_sum$PERC <- WCBS_sum$LANDFILL_AREA*100/1000000
colnames(WCBS_sum)[1] <- "Square"

ss_landfill_but <- merge(WCBS_sum, ss_label_but, all.y = T)


# write this out for easy figure-plotting

write.csv(ss_landfill_but, "../../PreppedData/all_butterflies.csv")


### we are finally ready to begin interrogating the data!

# first a very simple comparison

ss_landfill_but$Type <- relevel(ss_landfill_but$Type, ref = "Target")

matched.dat <- ss_landfill_but[which(ss_landfill_but$Type != "Neighbouring"), ]
grouped.dat <- ss_landfill_but[which(ss_landfill_but$Type != "Matched"), ]


# matched

hist(matched.dat$Species)

mod1rb <- glmer(Species ~ Type + (1|Group),
              family = poisson (link = "log"),
              data = matched.dat)

summary(mod1rb)
drop1(mod1rb, test = "Chi")

chkres(mod1rb, matched.dat$Type)


# grouped 

mod2rb <- glmer(Species ~ Type + (1|Group),
            family = poisson (link = "log"),
            data = grouped.dat)

summary(mod2rb)
drop1(mod2rb, test = "Chi")

chkres(mod2rb, grouped.dat$Type)


# combined

mod3rb <- glmer(Species ~ Type + (1|Group),
               family = poisson (link = "log"),
               data = ss_landfill_but)

summary(mod3rb)
drop1(mod3rb, test = "Chi")

chkres(mod3rb)

plot(ss_landfill_but$Species ~ ss_landfill_but$Type)
# appears that target sites have higher SR than either other category...


## repeat with estimated richness

# matched

hist(matched.dat$chao)

mod1eb <- glmer(round(chao) ~ Type + (1|Group),
              family = poisson (link = "log"),
              data = matched.dat)

summary(mod1eb)
drop1(mod1eb, test = "Chi")

chkres(mod1eb, matched.dat$Type)


# grouped 

mod2eb <- glmer(round(chao) ~ Type + (1|Group),
              family = poisson (link = "log"),
              data = grouped.dat)

summary(mod2eb)
drop1(mod2eb, test = "Chi")

chkres(mod2eb, grouped.dat$Type)


# combined

mod3eb <- glmer(round(chao) ~ Type + (1|Group),
               family = poisson (link = "log"),
               data = ss_landfill_but)

summary(mod3eb)
drop1(mod3eb, test = "Chi")

chkres(mod3eb)

plot(ss_landfill_but$chao ~ ss_landfill_but$Type)

# this also suggests that target sites have HIGHER SR than matched habitat, but in this case possibly LOWER than surrounding landscape

# some sign of the bias here that was present in the plant data, so
# let's check sampling completeness out

hist(ss_landfill_but$Completeness)


# matched

mod1cb <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data = matched.dat)

summary(mod1cb)
drop1(mod1cb, test = "Chi")

chkres(mod1cb, matched.dat$Type)


# grouped 

mod2cb <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data = grouped.dat)

summary(mod2cb)
drop1(mod2cb, test = "Chi")

chkres(mod2cb, grouped.dat$Type)


# combined

mod3cb <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data = ss_landfill_but)

summary(mod3cb)
drop1(mod3cb, test = "Chi")

chkres(mod3cb)

plot(ss_landfill_but$Completeness ~ ss_landfill_but$Type)

# no cause for concern here
# maybe slightly lower (but not significant) in neighbouring sites



## species richness is significantly lower in brownfield-type habitats than other habitat types, 
# but significantly higher in brownfield squares than non-brownfields of same type

# in other words, the habitats associated with brownfield sites have lower species richness than other local habitat types
# but brownfield sites have higher SR than other sites with similar habitat...?


# repeat for Irr

# matched

hist(matched.dat$Irr.cjm)

mod1ib <- glmer(Irr.cjm ~ Type + (1|Group),
            family = binomial (link = "logit"),
            weights = n,
            data = matched.dat)

summary(mod1ib)
drop1(mod1ib, test = "Chi")

chkres(mod1ib, matched.dat$Type)


# grouped

mod2ib <- glmer(Irr.cjm ~ Type + (1|Group),
              family = binomial (link = "logit"),
              weights = n,
              data=grouped.dat)

summary(mod2ib)
drop1(mod2ib, test = "Chi")

chkres(mod2ib)


# combined

mod3ib <- glmer(Irr.cjm ~ Type + (1|Group),
               family = binomial (link = "logit"),
               weights = n,
               data=ss_landfill_but)

summary(mod3ib)
drop1(mod3ib, test = "Chi")

chkres(mod3ib)

plot(ss_landfill_but$Irr.cjm ~ ss_landfill_but$Type)


# like plants, no sign of a rarity effect


## the next set of questions relates to the effect of *area of brownfield* within target sites


ss_brownfield_but <- ss_landfill_but[which(ss_landfill_but$PERC > 5), ]


# first, estimated SR

mod4srb <- glm(round(chao) ~ LANDFILL_AREA,
                 family = poisson (link = "log"),
                 data = ss_brownfield_but)

summary(mod4srb)
drop1(mod4srb, test = "Chi")

chkres(mod4srb)


plot(log(ss_brownfield_but$chao) ~ log(ss_brownfield_but$LANDFILL_AREA))


# next, completeness (just to check)

mod4cb <-  glm(Completeness ~ LANDFILL_AREA,
              data = ss_brownfield_but)

summary(mod4cb)
drop1(mod4cb, test = "Chi")

chkres(mod4cb)

plot(ss_brownfield_but$Completeness ~ log(ss_brownfield_but$LANDFILL_AREA))

# no problem here


mod4ib <-  glm(Irr.cjm ~ LANDFILL_AREA,
              family = binomial (link = "logit"),
              data = ss_brownfield_but)


summary(mod4ib)
drop1(mod4ib, test = "Chi")

chkres(mod4ib)


plot(ss_brownfield_but$Irr.cjm ~ log(ss_brownfield_but$LANDFILL_AREA))


## for completeness, do these analyses with the percentage landfill by area (instead of raw area)

# first, estimated SR

mod5srb <- glm(round(chao) ~ log(PERC),
              family = poisson (link = "log"),
              data = ss_brownfield_but)

summary(mod5srb)
drop1(mod5srb, test = "Chi")

chkres(mod5srb)

plot(log(ss_brownfield_but$chao) ~ log(ss_brownfield_but$PERC))


# and observed SR

mod5srb <- glm(Richness ~ log(PERC),
               family = poisson (link = "log"),
               data = ss_brownfield_but)

summary(mod5srb)
drop1(mod5srb, test = "Chi")

chkres(mod5srb)

plot(log(ss_brownfield_but$chao) ~ log(ss_brownfield_but$PERC))



# next, completeness (just to check)

mod5cb <-  glm(Completeness ~ log(PERC),
              data = ss_brownfield_but)

summary(mod5cb)
drop1(mod5cb, test = "Chi")

chkres(mod5cb)

plot(ss_brownfield_but$Completeness ~ log(ss_brownfield_but$PERC))


# unlike last time, this appears to be fine
# so, moving on to Irr

mod5ib <-  glm(Irr.cjm ~ log(PERC),
              family = binomial (link = "logit"),
              data = ss_brownfield_but)


summary(mod5ib)
drop1(mod5ib, test = "Chi")

chkres(mod5ib)


plot(ss_brownfield_but$Irr.cjm ~ log(ss_brownfield_but$PERC))



## and now we have a set of questions relating to the age of the brownfield sites

ss_brownfield_age_but <- ss_brownfield_but[which(!is.na(ss_brownfield_but$TimeSinceClosure)), ]


# first, estimated SR

mod6srb <- glm(round(chao) ~ TSCY,
              family = poisson (link = "log"),
              data = ss_brownfield_age_but)

summary(mod6srb)
drop1(mod6srb, test = "Chi")

chkres(mod6srb)


plot(log(ss_brownfield_age_but$chao) ~ ss_brownfield_age_but$TSCY)

# and observed

mod6sob <- glm(Species ~ TSCY,
              family = poisson (link = "log"),
              data = ss_brownfield_age_but)

summary(mod6sob)
drop1(mod6sob, test = "Chi")

chkres(mod6sob)

plot(log(ss_brownfield_age_but$Species) ~ ss_brownfield_age_but$TSCY)


# next, completeness (just to check)

mod6cb <-  glm(Completeness ~ TSCY,
              data = ss_brownfield_age_but)

summary(mod6cb)
drop1(mod6cb, test = "Chi")

chkres(mod6cb)

plot(ss_brownfield_age_but$Completeness ~ ss_brownfield_age_but$TSCY)



# no problem here
# so, moving on to Irr

mod6ib <-  glm(Irr.cjm ~ TSCY,
              family = binomial (link = "logit"),
              data = ss_brownfield_age_but)


summary(mod6ib)
drop1(mod6ib, test = "Chi")

chkres(mod6ib)


plot(ss_brownfield_age_but$Irr.cjm ~ ss_brownfield_age_but$TimeSinceClosure)



### WIDER COUNTRYSIDE ONLY
# the analysis above used data on all butterflies
# however, with butterflies we have a convenient and established division between generalists and specialists (based on the Millennium Atlas)
# which gives us the chance to investigate whether the two groups respond differently

data_locs_wcb <- data_locs[which(data_locs$GROUP=="Common butterfly"), ]


## now we can start to pick out actual stats on these squares
# let's start very simple:
# figure out species richness (raw and extrapolated)


all.sqs.wcb <- levels(droplevels(data_locs_wcb$GRID_REF))

richness_wcb <- data.frame()


for(x in all.sqs.wcb){
  sq.dat <- data_locs_wcb[which(data_locs_wcb$GRID_REF == x), ]
  print(x)
  
  # first crack raw spp richness
  sq.dat$PRESENT <- 1
  
  sq.spp <- ddply(sq.dat, .(COMMON_NAME, BMSCODE), summarise,
                  RECS = sum(PRESENT),
                  TOT = sum(COUNT),
                  MEAN = mean(COUNT))
  
  # now extrapolate spp richness across samples
  
  sq.pres <- sq.dat[,c(4,6,12)]
  
  # wrangle the data to be in wide form (one row per sample)
  
  sq.mat <- dcast(sq.pres, VISIT_DATE ~ COMMON_NAME,
                  value.var = "COUNT",
                  fun.aggregate = sum)
  
  # get rid of the sample ID column, otherwise it gets treated as an extra species!
  sq.mat <- sq.mat[,-1]
  
  out <- specpool(sq.mat)
  out$Square <- x
  
  richness_wcb <- rbind(richness_wcb, out)
  
}

summary(richness_wcb)



#### figure out species rarity 
# use R package "Rarity"


### need a full site-level matrix (presence/absence, one row per square)

# trim out only the columns we need
data_locs_mat_wcb <- data_locs_wcb[,c(2,6,12)]
data_locs_mat_wcb$presence <- 1


data_mat_wcb <- dcast(data_locs_mat_wcb, GRID_REF ~ COMMON_NAME,
                      value.var = "presence",
                      fun.aggregate = sum)

# make rownames from col 1

rownames(data_mat_wcb) <- data_mat_wcb[,1]

data_mat_wcb <- data_mat_wcb[,-1]

# finally transpose the matrix so that assemblages are in columns and species in rows

data_mat_wcb_t <- t(data_mat_wcb)


### and a list of species with number of square-level occurrences

# first narrow down to one row per species*square
spp_occ_sq_wcb <- ddply(data_locs_mat_wcb, .(GRID_REF, COMMON_NAME), summarise,
                        presence = mean(presence))

summary(spp_occ_sq_wcb)

# then summarise the data for each species
spp_occ_wcb <- ddply(spp_occ_sq_wcb, .(COMMON_NAME), summarise,
                     squares = sum(presence))

# need to rejig the structure of this slightly

spp_occ_mat_wcb <- spp_occ_wcb$squares
names(spp_occ_mat_wcb) <- spp_occ_wcb$COMMON_NAME
head(spp_occ_mat_wcb)

### then need to apply the functions rWeights and Irr

rarity.weights_wcb <- rWeights(spp_occ_mat_wcb)

square_rarity_wcb <- data.frame(Irr(assemblages = data_mat_wcb_t, W = rarity.weights_wcb))
hist(square_rarity_wcb$Irr)

# this actually causes some problems further down the line, because this index generates a lot of species with weight zero
# so I'm going to try and generate my own simpler version for comparison...
# there are 2096 squares in total, so each species can be observed in 1-2096 squares
# for my index I'll take 1-(obs.sqs/2096)

rarity.weights_wcb$CJM <- 1-(rarity.weights_wcb$Q/2096)

# normalise these between 0 and 1 (slightly boosting the rare species and taking away from the common)
minW <- min(rarity.weights_wcb$CJM)
maxW <- max(rarity.weights_wcb$CJM)

rarity.weights_wcb$CJMz <- (rarity.weights_wcb$CJM - minW)/(maxW - minW)

summary(rarity.weights_wcb$CJMz)

rarity.weights.cjm_wcb <- rarity.weights_wcb[,c(1,2,6,4)]
colnames(rarity.weights.cjm_wcb)[3] <- "W"

square_rarity.cjm_wcb <- data.frame(Irr(assemblages = data_mat_wcb_t, W = rarity.weights.cjm_wcb))
hist(square_rarity.cjm_wcb$Irr)


### now merge all these indices back together so that we have one row per square in a single dataframe
square_rarity_wcb$Square <- rownames(square_rarity_wcb)
square_rarity.cjm_wcb$Square <- rownames(square_rarity.cjm_wcb)

colnames(square_rarity.cjm_wcb) <- c("Irr.cjm", "Richness.cjm", "Square")

square_rr_wcb <- merge(square_rarity_wcb, square_rarity.cjm_wcb)

square_stats_wcb <- merge(richness_wcb, square_rr_wcb)

# restore squares with zero richness, and fill in blanks
square_stats_complete_wcb <- merge(all.sqs, square_stats_wcb, all = T)

square_stats_complete_wcb[is.na(square_stats_complete_wcb)] <- 0




### now we are ready to start some basic analyses
# first pick out the stats for the squares in pairs

# several ways to tackle this
# I'll want a long-form version (still one row for every square, wcb labelling which group they belong to)
# and possibly a wide-form version (one row per brownfield square, columns for )

# first, add a factor label to each square

ss_label_wcb <- data.frame()

y <- 0

for (x in squares){
  print(x)
  y <- y+1
  
  sq.pair <- pairs[which(pairs$SQUARE == x), ]
  
  sq.sqs <- append(x, levels(droplevels(sq.pair$NEIGHBOUR)))
  
  sq.dat <- square_stats_complete_wcb[which(square_stats_complete_wcb$Square %in% sq.sqs), ]
  
  
  sq.pair.base <- sq.pair[,c('NEIGHBOUR','TYPE')]
  sq.pair.base$NEIGHBOUR <- as.character(sq.pair.base$NEIGHBOUR)
  sq.pair.base[5,] <- c(x, "Target")
  colnames(sq.pair.base) <- c("Square", "Type")
  
  
  sq.dat$Group <- y
  
  sq.type <- merge(sq.pair.base, sq.dat)
  
  
  ss_label_wcb <- rbind(ss_label_wcb, sq.type)
  
  
}

summary(ss_label_wcb)

ss_label_wcb$Type <- as.factor(ss_label_wcb$Type)
ss_label_wcb$Group <- as.factor(ss_label_wcb$Group)

ss_label_wcb$Completeness <- ss_label_wcb$Species/ss_label_wcb$chao

summary(ss_label_wcb)


# tack in the area of brownfield site

ss_landfill_wcb <- merge(WCBS_sum, ss_label_wcb, all.y = T)



# write this out for easy figure-plotting

write.csv(ss_landfill_wcb, "../../PreppedData/all_wc_butterflies.csv")


### we are finally ready to begin interrogating the data!

# first a very simple comparison

ss_landfill_wcb$Type <- relevel(ss_landfill_wcb$Type, ref = "Target")

matched.dat <- ss_landfill_wcb[which(ss_landfill_wcb$Type != "Neighbouring"), ]
grouped.dat <- ss_landfill_wcb[which(ss_landfill_wcb$Type != "Matched"), ]


# matched

hist(matched.dat$Species)

mod1rw <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = matched.dat)

summary(mod1rw)
drop1(mod1rw, test = "Chi")

chkres(mod1rw, matched.dat$Type)


# grouped 

mod2rw <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = grouped.dat)

summary(mod2rw)
drop1(mod2rw, test = "Chi")

chkres(mod2rw, grouped.dat$Type)


# combined

mod3rw <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = ss_landfill_wcb)

summary(mod3rw)
drop1(mod3rw, test = "Chi")

chkres(mod3rw)

plot(ss_landfill_wcb$Species ~ ss_landfill_wcb$Type)
# appears that target sites have higher SR than either other category...


## repeat with estimated richness

# matched

hist(matched.dat$chao)

mod1ew <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = matched.dat)

summary(mod1ew)
drop1(mod1ew, test = "Chi")

chkres(mod1ew, matched.dat$Type)


# grouped 

mod2ew <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = grouped.dat)

summary(mod2ew)
drop1(mod2ew, test = "Chi")

chkres(mod2ew, grouped.dat$Type)


# combined

mod3ew <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = ss_landfill_wcb)

summary(mod3ew)
drop1(mod3ew, test = "Chi")

chkres(mod3ew)

plot(ss_landfill_wcb$chao ~ ss_landfill_wcb$Type)

# this also suggests that target sites have HIGHER SR than matched habitat, wcb in this case possibly LOWER than surrounding landscape

# some sign of the bias here that was present in the plant data, so
# let's check sampling completeness out

hist(ss_landfill_wcb$Completeness)

# matched

mod1cw <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
               data = matched.dat)

summary(mod1cw)
drop1(mod1cw, test = "Chi")

chkres(mod1cw, matched.dat$Type)



# grouped

mod2cw <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
               data = grouped.dat)

summary(mod2cw)
drop1(mod2cw, test = "Chi")

chkres(mod2cw, grouped.dat$Type)



# combined

mod3cw<- glmer(Completeness ~ Type + (1|Group),
               family = binomial (link = "logit"),
               weights = n,
               data = ss_landfill_wcb)

summary(mod3cw)
drop1(mod3cw, test = "Chi")

chkres(mod3cw)

plot(ss_landfill_wcb$Completeness ~ ss_landfill_wcb$Type)

# no cause for concern here
# maybe slightly lower (but not significant) in neighbouring sites



## species richness is significantly lower in brownfield-type habitats than other habitat types, 
# but significantly higher in brownfield squares than non-brownfields of same type

# in other words, the habitats associated with brownfield sites have lower species richness than other local habitat types
# but brownfield sites have higher SR than other sites with similar habitat...?


# repeat for Irr

# matched

hist(matched.dat$Irr.cjm)

mod1iw <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data = matched.dat)

summary(mod1iw)
drop1(mod1iw, test = "Chi")

chkres(mod1iw, matched.dat$Type)


# grouped

mod2iw <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data=grouped.dat)

summary(mod2iw)
drop1(mod2iw, test = "Chi")

chkres(mod2iw)


# combined

mod3iw <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data=ss_landfill_wcb)

summary(mod3iw)
drop1(mod3iw, test = "Chi")

chkres(mod3iw)

plot(ss_landfill_wcb$Irr.cjm ~ ss_landfill_wcb$Type)


# like plants, no sign of a rarity effect


## the next set of questions relates to the effect of *area of brownfield* within target sites


ss_brownfield_wcb <- ss_landfill_wcb[which(ss_landfill_wcb$PERC > 5), ]


# first, estimated SR

mod4srw <- glm(round(chao) ~ LANDFILL_AREA,
               family = poisson (link = "log"),
               data = ss_brownfield_wcb)

summary(mod4srw)
drop1(mod4srw, test = "Chi")

chkres(mod4srw)


plot(log(ss_brownfield_wcb$chao) ~ log(ss_brownfield_wcb$LANDFILL_AREA))


# next, completeness (just to check)

mod4cw <-  glm(Completeness ~ LANDFILL_AREA,
               data = ss_brownfield_wcb)

summary(mod4cw)
drop1(mod4cw, test = "Chi")

chkres(mod4cw)

plot(ss_brownfield_wcb$Completeness ~ log(ss_brownfield_wcb$LANDFILL_AREA))

# no problem here


mod4iw <-  glm(Irr.cjm ~ LANDFILL_AREA,
               family = binomial (link = "logit"),
               data = ss_brownfield_wcb)


summary(mod4iw)
drop1(mod4iw, test = "Chi")

chkres(mod4iw)


plot(ss_brownfield_wcb$Irr.cjm ~ log(ss_brownfield_wcb$LANDFILL_AREA))


## for completeness, do these analyses with the percentage landfill by area (instead of raw area)

# first, estimated SR

mod5srw <- glm(round(chao) ~ log(PERC),
               family = poisson (link = "log"),
               data = ss_brownfield_wcb)

summary(mod5srw)
drop1(mod5srw, test = "Chi")

chkres(mod5srw)

plot(log(ss_brownfield_wcb$chao) ~ log(ss_brownfield_wcb$PERC))

# and observed SR

mod5sow <- glm(Species ~ log(PERC),
               family = poisson (link = "log"),
               data = ss_brownfield_wcb)

summary(mod5sow)
drop1(mod5sow, test = "Chi")

chkres(mod5sow)

plot(log(ss_brownfield_wcb$Species) ~ log(ss_brownfield_wcb$PERC))



# next, completeness (just to check)

mod5cw <-  glm(Completeness ~ log(PERC),
               data = ss_brownfield_wcb)

summary(mod5cw)
drop1(mod5cw, test = "Chi")

chkres(mod5cw)

plot(ss_brownfield_wcb$Completeness ~ log(ss_brownfield_wcb$PERC))


# unlike last time, this appears to be fine
# so, moving on to Irr

mod5iw <-  glm(Irr.cjm ~ log(PERC),
               family = binomial (link = "logit"),
               data = ss_brownfield_wcb)


summary(mod5iw)
drop1(mod5iw, test = "Chi")

chkres(mod5iw)


plot(ss_brownfield_wcb$Irr.cjm ~ log(ss_brownfield_wcb$PERC))



## and now we have a set of questions relating to the age of the brownfield sites

ss_brownfield_age_wcb <- ss_brownfield_wcb[which(!is.na(ss_brownfield_wcb$TimeSinceClosure)), ]


# first, estimated SR

mod6srw <- glm(round(chao) ~ TSCY,
               family = poisson (link = "log"),
               data = ss_brownfield_age_wcb)

summary(mod6srw)
drop1(mod6srw, test = "Chi")

chkres(mod6srw)


plot(log(ss_brownfield_age_wcb$chao) ~ ss_brownfield_age_wcb$TSCY)

# and observed

mod6sow <- glm(Species ~ TSCY,
               family = poisson (link = "log"),
               data = ss_brownfield_age_wcb)

summary(mod6sow)
drop1(mod6sow, test = "Chi")

chkres(mod6sow)

plot(log(ss_brownfield_age_wcb$Species) ~ ss_brownfield_age_wcb$TSCY)


# next, completeness (just to check)

mod6cw <-  glm(Completeness ~ TSCY,
               data = ss_brownfield_age_wcb)

summary(mod6cw)
drop1(mod6cw, test = "Chi")

chkres(mod6cw)

plot(ss_brownfield_age_wcb$Completeness ~ ss_brownfield_age_wcb$TSCY)



# no problem here
# so, moving on to Irr

mod6iw <-  glm(Irr.cjm ~ TSCY,
               family = binomial (link = "logit"),
               data = ss_brownfield_age_wcb)


summary(mod6iw)
drop1(mod6iw, test = "Chi")

chkres(mod6iw)


plot(ss_brownfield_age_wcb$Irr.cjm ~ ss_brownfield_age_wcb$TimeSinceClosure)



### repeat for the minor groups-

# other butterflies


data_locs_oth <- data_locs[which(data_locs$GROUP=="Other butterfly"), ]


## now we can start to pick out actual stats on these squares
# let's start very simple:
# figure out species richness (raw and extrapolated)


all.sqs.oth <- levels(droplevels(data_locs_oth$GRID_REF))

richness_oth <- data.frame()


for(x in all.sqs.oth){
  sq.dat <- data_locs_oth[which(data_locs_oth$GRID_REF == x), ]
  print(x)
  
  # first crack raw spp richness
  sq.dat$PRESENT <- 1
  
  sq.spp <- ddply(sq.dat, .(COMMON_NAME, BMSCODE), summarise,
                  RECS = sum(PRESENT),
                  TOT = sum(COUNT),
                  MEAN = mean(COUNT))
  
  # now extrapolate spp richness across samples
  
  sq.pres <- sq.dat[,c(4,6,12)]
  
  # wrangle the data to be in wide form (one row per sample)
  
  sq.mat <- dcast(sq.pres, VISIT_DATE ~ COMMON_NAME,
                  value.var = "COUNT",
                  fun.aggregate = sum)
  
  # get rid of the sample ID column, otherwise it gets treated as an extra species!
  sq.mat <- sq.mat[,-1]
  
  out <- specpool(sq.mat)
  out$Square <- x
  
  richness_oth <- rbind(richness_oth, out)
  
}

summary(richness_oth)



#### figure out species rarity 
# use R package "Rarity"


### need a full site-level matrix (presence/absence, one row per square)

# trim out only the columns we need
data_locs_mat_oth <- data_locs_oth[,c(2,6,12)]
data_locs_mat_oth$presence <- 1


data_mat_oth <- dcast(data_locs_mat_oth, GRID_REF ~ COMMON_NAME,
                      value.var = "presence",
                      fun.aggregate = sum)

# make rownames from col 1

rownames(data_mat_oth) <- data_mat_oth[,1]

data_mat_oth <- data_mat_oth[,-1]

# finally transpose the matrix so that assemblages are in columns and species in rows

data_mat_oth_t <- t(data_mat_oth)


### and a list of species with number of square-level occurrences

# first narrow down to one row per species*square
spp_occ_sq_oth <- ddply(data_locs_mat_oth, .(GRID_REF, COMMON_NAME), summarise,
                        presence = mean(presence))

summary(spp_occ_sq_oth)

# then summarise the data for each species
spp_occ_oth <- ddply(spp_occ_sq_oth, .(COMMON_NAME), summarise,
                     squares = sum(presence))

# need to rejig the structure of this slightly

spp_occ_mat_oth <- spp_occ_oth$squares
names(spp_occ_mat_oth) <- spp_occ_oth$COMMON_NAME
head(spp_occ_mat_oth)

### then need to apply the functions rWeights and Irr

rarity.weights_oth <- rWeights(spp_occ_mat_oth)

square_rarity_oth <- data.frame(Irr(assemblages = data_mat_oth_t, W = rarity.weights_oth))
hist(square_rarity_oth$Irr)

# this actually causes some problems further down the line, because this index generates a lot of species with weight zero
# so I'm going to try and generate my own simpler version for comparison...
# there are 2096 squares in total, so each species can be observed in 1-2096 squares
# for my index I'll take 1-(obs.sqs/2096)

rarity.weights_oth$CJM <- 1-(rarity.weights_oth$Q/2096)

# normalise these between 0 and 1 (slightly boosting the rare species and taking away from the common)
minW <- min(rarity.weights_oth$CJM)
maxW <- max(rarity.weights_oth$CJM)

rarity.weights_oth$CJMz <- (rarity.weights_oth$CJM - minW)/(maxW - minW)

summary(rarity.weights_oth$CJMz)

rarity.weights.cjm_oth <- rarity.weights_oth[,c(1,2,6,4)]
colnames(rarity.weights.cjm_oth)[3] <- "W"

square_rarity.cjm_oth <- data.frame(Irr(assemblages = data_mat_oth_t, W = rarity.weights.cjm_oth))
hist(square_rarity.cjm_oth$Irr)


### now merge all these indices back together so that we have one row per square in a single dataframe
square_rarity_oth$Square <- rownames(square_rarity_oth)
square_rarity.cjm_oth$Square <- rownames(square_rarity.cjm_oth)

colnames(square_rarity.cjm_oth) <- c("Irr.cjm", "Richness.cjm", "Square")

square_rr_oth <- merge(square_rarity_oth, square_rarity.cjm_oth)

square_stats_oth <- merge(richness_oth, square_rr_oth)

# restore squares with zero richness, and fill in blanks
square_stats_complete_oth <- merge(all.sqs, square_stats_oth, all = T)

square_stats_complete_oth[is.na(square_stats_complete_oth)] <- 0



### now we are ready to start some basic analyses
# first pick out the stats for the squares in pairs

# several ways to tackle this
# I'll want a long-form version (still one row for every square, oth labelling which group they belong to)
# and possibly a wide-form version (one row per brownfield square, columns for )

# first, add a factor label to each square

ss_label_oth <- data.frame()

y <- 0

for (x in squares){
  print(x)
  y <- y+1
  
  sq.pair <- pairs[which(pairs$SQUARE == x), ]
  
  sq.sqs <- append(x, levels(droplevels(sq.pair$NEIGHBOUR)))
  
  sq.dat <- square_stats_complete_oth[which(square_stats_complete_oth$Square %in% sq.sqs), ]
  
  
  sq.pair.base <- sq.pair[,c('NEIGHBOUR','TYPE')]
  sq.pair.base$NEIGHBOUR <- as.character(sq.pair.base$NEIGHBOUR)
  sq.pair.base[5,] <- c(x, "Target")
  colnames(sq.pair.base) <- c("Square", "Type")
  
  
  sq.dat$Group <- y
  
  sq.type <- merge(sq.pair.base, sq.dat)
  
  
  ss_label_oth <- rbind(ss_label_oth, sq.type)
  
  
}

summary(ss_label_oth)

ss_label_oth$Type <- as.factor(ss_label_oth$Type)
ss_label_oth$Group <- as.factor(ss_label_oth$Group)

ss_label_oth$Completeness <- ss_label_oth$Species/ss_label_oth$chao

summary(ss_label_oth)


# tack in the area of brownfield site

ss_landfill_oth <- merge(WCBS_sum, ss_label_oth, all.y = T)


ss_landfill_oth$LANDFILL_AREA[is.na(ss_landfill_oth$LANDFILL_AREA)] <- 0
ss_landfill_oth$PERC[is.na(ss_landfill_oth$PERC)] <- 0

# write this out for easy figure-plotting

write.csv(ss_landfill_oth, "../../PreppedData/all_hs_butterflies.csv")




### we are finally ready to begin interrogating the data!

# first a very simple comparison

ss_landfill_oth$Type <- relevel(ss_landfill_oth$Type, ref = "Target")

matched.dat <- ss_landfill_oth[which(ss_landfill_oth$Type != "Neighbouring"), ]
grouped.dat <- ss_landfill_oth[which(ss_landfill_oth$Type != "Matched"), ]


# matched

hist(matched.dat$Species)

mod1ro <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = matched.dat)

summary(mod1ro)
drop1(mod1ro, test = "Chi")

chkres(mod1ro, matched.dat$Type)


# grouped 

mod2ro <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = grouped.dat)

summary(mod2ro)
drop1(mod2ro, test = "Chi")

chkres(mod2ro, grouped.dat$Type)


# combined

mod3ro <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = ss_landfill_oth)

summary(mod3ro)
drop1(mod3ro, test = "Chi")

chkres(mod3ro)

plot(ss_landfill_oth$Species ~ ss_landfill_oth$Type)

# appears that neighbouring sites have higher SR than either other category... but not significantly


## repeat with estimated richness

# matched

hist(matched.dat$chao)

mod1eo <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = matched.dat)

summary(mod1eo)
drop1(mod1eo, test = "Chi")

chkres(mod1eo, matched.dat$Type)


# grouped 

mod2eo <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = grouped.dat)

summary(mod2eo)
drop1(mod2eo, test = "Chi")

chkres(mod2eo, grouped.dat$Type)


# combined

mod3eo <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = ss_landfill_oth)

summary(mod3eo)
drop1(mod3eo, test = "Chi")

chkres(mod3eo)

plot(ss_landfill_oth$chao ~ ss_landfill_oth$Type)

# this also suggests that neighbouring sites have HIGHER SR than target or matched habitat

# let's check sampling completeness out

hist(ss_landfill_oth$Completeness)

# matched

mod1co <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
               data = matched.dat)

summary(mod1co)
drop1(mod1co, test = "Chi")

chkres(mod1co, matched.dat$Type)



# grouped

mod2co <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
               data = grouped.dat)

summary(mod2co)
drop1(mod2co, test = "Chi")

chkres(mod2co, grouped.dat$Type)



# combined

mod3co <- glmer(Completeness ~ Type + (1|Group),
               family = binomial (link = "logit"),
               weights = n,
               data = ss_landfill_oth)

summary(mod3co)
drop1(mod3co, test = "Chi")

chkres(mod3co)

plot(ss_landfill_oth$Completeness ~ ss_landfill_oth$Type)


## species richness is significantly lower in brownfield-type habitats than other habitat types, 
# but similar in brownfield squares and non-brownfields of same type

# in other words, the habitats associated with brownfield sites have lower species richness than other local habitat types


# repeat for Irr

# matched

hist(matched.dat$Irr.cjm)

mod1io <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data = matched.dat)

summary(mod1io)
drop1(mod1io, test = "Chi")

chkres(mod1io, matched.dat$Type)


# grouped

mod2io <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data=grouped.dat)

summary(mod2io)
drop1(mod2io, test = "Chi")

chkres(mod2io)


# combined

mod3io <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data=ss_landfill_oth)

summary(mod3io)
drop1(mod3io, test = "Chi")

chkres(mod3io)

plot(ss_landfill_oth$Irr.cjm ~ ss_landfill_oth$Type)


# a marginal signal of an effect where brownfield sites lack rare species
# but it's not statistically significant


## the next set of questions relates to the effect of *area of brownfield* within target sites


ss_brownfield_oth <- ss_landfill_oth[which(ss_landfill_oth$PERC > 5), ]


# first, estimated SR

mod4sro <- glm(round(chao) ~ LANDFILL_AREA,
               family = poisson (link = "log"),
               data = ss_brownfield_oth)

summary(mod4sro)
drop1(mod4sro, test = "Chi")

chkres(mod4sro)


plot(log(ss_brownfield_oth$chao) ~ log(ss_brownfield_oth$LANDFILL_AREA))


# next, completeness (just to check)

mod4co <-  glm(Completeness ~ LANDFILL_AREA,
               data = ss_brownfield_oth)

summary(mod4co)
drop1(mod4co, test = "Chi")

chkres(mod4co)

plot(ss_brownfield_oth$Completeness ~ log(ss_brownfield_oth$LANDFILL_AREA))

# no problem here


mod4io <-  glm(Irr.cjm ~ LANDFILL_AREA,
               family = binomial (link = "logit"),
               data = ss_brownfield_oth)


summary(mod4io)
drop1(mod4io, test = "Chi")

chkres(mod4io)


plot(ss_brownfield_oth$Irr.cjm ~ log(ss_brownfield_oth$LANDFILL_AREA))


## for completeness, do these analyses with the percentage landfill by area (instead of raw area)

# first, estimated SR

mod5sro <- glm(round(chao) ~ log(PERC),
               family = poisson (link = "log"),
               data = ss_brownfield_oth)

summary(mod5sro)
drop1(mod5sro, test = "Chi")

chkres(mod5sro)

plot(log(ss_brownfield_oth$chao) ~ log(ss_brownfield_oth$PERC))

# and observed

mod5soo <- glm(Species ~ log(PERC),
               family = poisson (link = "log"),
               data = ss_brownfield_oth)

summary(mod5soo)
drop1(mod5soo, test = "Chi")

chkres(mod5soo)

plot(log(ss_brownfield_oth$chao) ~ log(ss_brownfield_oth$PERC))




# next, completeness (just to check)

mod5co <-  glm(Completeness ~ log(PERC),
               data = ss_brownfield_oth)

summary(mod5co)
drop1(mod5co, test = "Chi")

chkres(mod5co)

plot(ss_brownfield_oth$Completeness ~ log(ss_brownfield_oth$PERC))



# so, moving on to Irr

mod5io <-  glm(Irr.cjm ~ log(PERC),
               family = binomial (link = "logit"),
               data = ss_brownfield_oth)


summary(mod5io)
drop1(mod5io, test = "Chi")

chkres(mod5io)


plot(ss_brownfield_oth$Irr.cjm ~ log(ss_brownfield_oth$PERC))


## and now we have a set of questions relating to the age of the brownfield sites

ss_brownfield_age_oth <- ss_brownfield_oth[which(!is.na(ss_brownfield_oth$TimeSinceClosure)), ]


# first, estimated SR

mod6sro <- glm(round(chao) ~ TSCY,
               family = poisson (link = "log"),
               data = ss_brownfield_age_oth)

summary(mod6sro)
drop1(mod6sro, test = "Chi")

chkres(mod6sro)


plot(log(ss_brownfield_age_oth$chao) ~ ss_brownfield_age_oth$TSCY)

# and observed

mod6soo <- glm(Species ~ TSCY,
               family = poisson (link = "log"),
               data = ss_brownfield_age_oth)

summary(mod6soo)
drop1(mod6soo, test = "Chi")

chkres(mod6soo)

plot(log(ss_brownfield_age_oth$Species) ~ ss_brownfield_age_oth$TSCY)


# next, completeness (just to check)

mod6co <-  glm(Completeness ~ TSCY,
               data = ss_brownfield_age_oth)

summary(mod6co)
drop1(mod6co, test = "Chi")

chkres(mod6co)

plot(ss_brownfield_age_oth$Completeness ~ ss_brownfield_age_oth$TSCY)



# no problem here
# so, moving on to Irr

mod6io <-  glm(Irr.cjm ~ TSCY,
               family = binomial (link = "logit"),
               data = ss_brownfield_age_oth)


summary(mod6io)
drop1(mod6io, test = "Chi")

chkres(mod6io)


plot(ss_brownfield_age_oth$Irr.cjm ~ ss_brownfield_age_oth$TimeSinceClosure)







# moths


data_locs_mot <- data_locs[which(data_locs$GROUP=="Moth"), ]


## now we can start to pick out actual stats on these squares
# let's start very simple:
# figure out species richness (raw and extrapolated)


all.sqs.mot <- levels(droplevels(data_locs_mot$GRID_REF))

richness_mot <- data.frame()


for(x in all.sqs.mot){
  sq.dat <- data_locs_mot[which(data_locs_mot$GRID_REF == x), ]
  print(x)
  
  # first crack raw spp richness
  sq.dat$PRESENT <- 1
  
  sq.spp <- ddply(sq.dat, .(COMMON_NAME, BMSCODE), summarise,
                  RECS = sum(PRESENT),
                  TOT = sum(COUNT),
                  MEAN = mean(COUNT))
  
  # now extrapolate spp richness across samples
  
  sq.pres <- sq.dat[,c(4,6,12)]
  
  # wrangle the data to be in wide form (one row per sample)
  
  sq.mat <- dcast(sq.pres, VISIT_DATE ~ COMMON_NAME,
                  value.var = "COUNT",
                  fun.aggregate = sum)
  
  # get rid of the sample ID column, otherwise it gets treated as an extra species!
  sq.mat <- sq.mat[,-1]
  
  out <- specpool(sq.mat)
  out$Square <- x
  
  richness_mot <- rbind(richness_mot, out)
  
}

summary(richness_mot)



#### figure out species rarity 
# use R package "Rarity"


### need a full site-level matrix (presence/absence, one row per square)

# trim out only the columns we need
data_locs_mat_mot <- data_locs_mot[,c(2,6,12)]
data_locs_mat_mot$presence <- 1


data_mat_mot <- dcast(data_locs_mat_mot, GRID_REF ~ COMMON_NAME,
                      value.var = "presence",
                      fun.aggregate = sum)

# make rownames from col 1

rownames(data_mat_mot) <- data_mat_mot[,1]

data_mat_mot <- data_mat_mot[,-1]

# finally transpose the matrix so that assemblages are in columns and species in rows

data_mat_mot_t <- t(data_mat_mot)


### and a list of species with number of square-level occurrences

# first narrow down to one row per species*square
spp_occ_sq_mot <- ddply(data_locs_mat_mot, .(GRID_REF, COMMON_NAME), summarise,
                        presence = mean(presence))

summary(spp_occ_sq_mot)

# then summarise the data for each species
spp_occ_mot <- ddply(spp_occ_sq_mot, .(COMMON_NAME), summarise,
                     squares = sum(presence))

# need to rejig the structure of this slightly

spp_occ_mat_mot <- spp_occ_mot$squares
names(spp_occ_mat_mot) <- spp_occ_mot$COMMON_NAME
head(spp_occ_mat_mot)

### then need to apply the functions rWeights and Irr

rarity.weights_mot <- rWeights(spp_occ_mat_mot)

square_rarity_mot <- data.frame(Irr(assemblages = data_mat_mot_t, W = rarity.weights_mot))
hist(square_rarity_mot$Irr)

# this actually causes some problems further down the line, because this index generates a lot of species with weight zero
# so I'm going to try and generate my own simpler version for comparison...
# there are 2096 squares in total, so each species can be observed in 1-2096 squares
# for my index I'll take 1-(obs.sqs/2096)

rarity.weights_mot$CJM <- 1-(rarity.weights_mot$Q/2096)

# normalise these between 0 and 1 (slightly boosting the rare species and taking away from the common)
minW <- min(rarity.weights_mot$CJM)
maxW <- max(rarity.weights_mot$CJM)

rarity.weights_mot$CJMz <- (rarity.weights_mot$CJM - minW)/(maxW - minW)

summary(rarity.weights_mot$CJMz)

rarity.weights.cjm_mot <- rarity.weights_mot[,c(1,2,6,4)]
colnames(rarity.weights.cjm_mot)[3] <- "W"

square_rarity.cjm_mot <- data.frame(Irr(assemblages = data_mat_mot_t, W = rarity.weights.cjm_mot))
hist(square_rarity.cjm_mot$Irr)


### now merge all these indices back together so that we have one row per square in a single dataframe
square_rarity_mot$Square <- rownames(square_rarity_mot)
square_rarity.cjm_mot$Square <- rownames(square_rarity.cjm_mot)

colnames(square_rarity.cjm_mot) <- c("Irr.cjm", "Richness.cjm", "Square")

square_rr_mot <- merge(square_rarity_mot, square_rarity.cjm_mot)

square_stats_mot <- merge(richness_mot, square_rr_mot)

# restore squares with zero richness, and fill in blanks
square_stats_complete_mot <- merge(all.sqs, square_stats_mot, all = T)




### now we are ready to start some basic analyses
# first pick out the stats for the squares in pairs

# several ways to tackle this
# I'll want a long-form version (still one row for every square, mot labelling which group they belong to)
# and possibly a wide-form version (one row per brownfield square, columns for )

# first, add a factor label to each square

ss_label_mot <- data.frame()

y <- 0

for (x in squares){
  print(x)
  y <- y+1
  
  sq.pair <- pairs[which(pairs$SQUARE == x), ]
  
  sq.sqs <- append(x, levels(droplevels(sq.pair$NEIGHBOUR)))
  
  sq.dat <- square_stats_complete_mot[which(square_stats_complete_mot$Square %in% sq.sqs), ]
  
  
  sq.pair.base <- sq.pair[,c('NEIGHBOUR','TYPE')]
  sq.pair.base$NEIGHBOUR <- as.character(sq.pair.base$NEIGHBOUR)
  sq.pair.base[5,] <- c(x, "Target")
  colnames(sq.pair.base) <- c("Square", "Type")
  
  
  sq.dat$Group <- y
  
  sq.type <- merge(sq.pair.base, sq.dat)
  
  
  ss_label_mot <- rbind(ss_label_mot, sq.type)
  
  
}

summary(ss_label_mot)

ss_label_mot$Type <- as.factor(ss_label_mot$Type)
ss_label_mot$Group <- as.factor(ss_label_mot$Group)

ss_label_mot$Completeness <- ss_label_mot$Species/ss_label_mot$chao

summary(ss_label_mot)


# tack in the area of brownfield site

ss_landfill_mot <- merge(WCBS_sum, ss_label_mot, all.y = T)


ss_landfill_mot$LANDFILL_AREA[is.na(ss_landfill_mot$LANDFILL_AREA)] <- 0
ss_landfill_mot$PERC[is.na(ss_landfill_mot$PERC)] <- 0

# write this out for easy figure-plotting

write.csv(ss_landfill_mot, "../../PreppedData/all_moths.csv")


### we are finally ready to begin interrogating the data!

# first a very simple comparison

ss_landfill_mot$Type <- relevel(ss_landfill_mot$Type, ref = "Target")

matched.dat <- ss_landfill_mot[which(ss_landfill_mot$Type != "Neighbouring"), ]
grouped.dat <- ss_landfill_mot[which(ss_landfill_mot$Type != "Matched"), ]


# matched

hist(matched.dat$Species)

mod1rm <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = matched.dat)

summary(mod1rm)
drop1(mod1rm, test = "Chi")

chkres(mod1rm, matched.dat$Type)


# grouped 

mod2rm <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = grouped.dat)

summary(mod2rm)
drop1(mod2rm, test = "Chi")

chkres(mod2rm, grouped.dat$Type)


# combined

mod3rm <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = ss_landfill_mot)

summary(mod3rm)
drop1(mod3rm, test = "Chi")

chkres(mod3rm)

plot(ss_landfill_mot$Species ~ ss_landfill_mot$Type)
# appears that target sites have higher SR than either other category...


## repeat with estimated richness

# matched

hist(matched.dat$chao)

mod1em <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = matched.dat)

summary(mod1em)
drop1(mod1em, test = "Chi")

chkres(mod1em)


# grouped 

mod2em <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = grouped.dat)

summary(mod2em)
drop1(mod2em, test = "Chi")

chkres(mod2em)


# combined

mod3em <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "identity"),
                data = ss_landfill_mot)


summary(mod3em)
drop1(mod3em, test = "Chi")

chkres(mod3em)

plot(log(ss_landfill_mot$chao) ~ ss_landfill_mot$Type)


# this model looks overfitted but Poisson is definitely the theoretically most-appropriate model

overdisp.glmer(mod3em)

# this indicates *underdispersion*
# one way to tackle this is with a QP model


model3emQ <- glmmPQL(chao ~ Type,
                   random = list(~1|Group),
                   family = quasipoisson (link = "log"),
                   data = ss_landfill_mot)


summary(model3emQ)
Anova(model3emQ, type = "III")  # drop1 doesn't work properly with this model class so we use a Type III Anova instead

# check residuals
chkres.PQL(model3emQ)


# this seems to confirm the effect, with much more plausible errors

# this also suggests that target sites have HIGHER SR than matched habitat, mot in this case possibly LOWER than surrounding landscape

# do the QP model for the separate pairwise


model1emQ <- glmmPQL(chao ~ Type,
                     random = list(~1|Group),
                     family = quasipoisson (link = "log"),
                     data = matched.dat)


summary(model1emQ)
Anova(model1emQ, type = "III")  # drop1 doesn't work properly with this model class so we use a Type III Anova instead

# check residuals
chkres.PQL(model1emQ)

model2emQ <- glmmPQL(chao ~ Type,
                     random = list(~1|Group),
                     family = quasipoisson (link = "log"),
                     data = grouped.dat)


summary(model2emQ)
Anova(model2emQ, type = "III")  # drop1 doesn't work properly with this model class so we use a Type III Anova instead

# check residuals
chkres.PQL(model2emQ)


# some sign of the bias here that was present in the plant data, so
# let's check sampling completeness out

hist(ss_landfill_mot$Completeness)

# matched

mod1cm <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
               data = matched.dat)

summary(mod1cm)
drop1(mod1cm, test = "Chi")

chkres(mod1cm, matched.dat$Type)



# grouped

mod2cm <- glmer(Completeness ~ Type + (1|Group),
               family = binomial (link = "logit"),
               weights = n,
               data = grouped.dat)

summary(mod2cm)
drop1(mod2cm, test = "Chi")

chkres(mod2cm, grouped.dat$Type)



# combined

mod3cm <- glmer(Completeness ~ Type + (1|Group), 
               family = binomial (link = "logit"),
               weights = n,
               data = ss_landfill_mot)

summary(mod3cm)
drop1(mod3cm, test = "Chi")

chkres(mod3cm)

plot(ss_landfill_mot$Completeness ~ ss_landfill_mot$Type)

# some cause for concern here
# maybe lower in target sites



## species richness is significantly higher in brownfield-type habitats than other habitat types, 
# and significantly higher in brownfield squares than non-brownfields of same type


# repeat for Irr

# matched

hist(matched.dat$Irr.cjm)

mod1im <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data = matched.dat)

summary(mod1im)
drop1(mod1im, test = "Chi")

chkres(mod1im, matched.dat$Type)


# grouped

mod2im <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data=grouped.dat)

summary(mod2im)
drop1(mod2im, test = "Chi")

chkres(mod2im)


# combined

mod3im <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data=ss_landfill_mot)

summary(mod3im)
drop1(mod3im, test = "Chi")

chkres(mod3im)

plot(ss_landfill_mot$Irr.cjm ~ ss_landfill_mot$Type)


# like plants, no sign of a rarity effect


## the next set of questions relates to the effect of *area of brownfield* within target sites


ss_brownfield_mot <- ss_landfill_mot[which(ss_landfill_mot$PERC > 5), ]


# first, estimated SR

mod4srm <- glm(round(chao) ~ LANDFILL_AREA,
               family = poisson (link = "log"),
               data = ss_brownfield_mot)

summary(mod4srm)
drop1(mod4srm, test = "Chi")

chkres(mod4srm)


plot(log(ss_brownfield_mot$chao) ~ log(ss_brownfield_mot$LANDFILL_AREA))


# next, completeness (just to check)

mod4cm <-  glm(Completeness ~ LANDFILL_AREA,
               data = ss_brownfield_mot)

summary(mod4cm)
drop1(mod4cm, test = "Chi")

chkres(mod4cm)

plot(ss_brownfield_mot$Completeness ~ log(ss_brownfield_mot$LANDFILL_AREA))

# no problem here


mod4im <-  glm(Irr.cjm ~ LANDFILL_AREA,
               family = binomial (link = "logit"),
               data = ss_brownfield_mot)


summary(mod4im)
drop1(mod4im, test = "Chi")

chkres(mod4im)


plot(ss_brownfield_mot$Irr.cjm ~ log(ss_brownfield_mot$LANDFILL_AREA))


## for completeness, do these analyses with the percentage landfill by area (instead of raw area)

# first, estimated SR

mod5srm <- glm(round(chao) ~ log(PERC),
               family = poisson (link = "log"),
               data = ss_brownfield_mot)

summary(mod5srm)
drop1(mod5srm, test = "Chi")

chkres(mod5srm)

plot(log(ss_brownfield_mot$chao) ~ log(ss_brownfield_mot$PERC))




# and observed
mod5som <- glm(Species ~ log(PERC),
               family = poisson (link = "log"),
               data = ss_brownfield_mot)

summary(mod5som)
drop1(mod5som, test = "Chi")

chkres(mod5som)

plot(log(ss_brownfield_mot$Species) ~ log(ss_brownfield_mot$PERC))




# next, completeness (just to check)

mod5cm <-  glm(Completeness ~ log(PERC),
               data = ss_brownfield_mot)

summary(mod5cm)
drop1(mod5cm, test = "Chi")

chkres(mod5cm)

plot(ss_brownfield_mot$Completeness ~ log(ss_brownfield_mot$PERC))


# unlike last time, this appears to be fine
# so, moving on to Irr

mod5im <-  glm(Irr.cjm ~ log(PERC),
               family = binomial (link = "logit"),
               data = ss_brownfield_mot)


summary(mod5im)
drop1(mod5im, test = "Chi")

chkres(mod5im)


plot(ss_brownfield_mot$Irr.cjm ~ log(ss_brownfield_mot$PERC))


## and now we have a set of questions relating to the age of the brownfield sites

ss_brownfield_age_mot <- ss_brownfield_mot[which(!is.na(ss_brownfield_mot$TimeSinceClosure)), ]


# first, estimated SR

mod6srm <- glm(round(chao) ~ TSCY,
               family = poisson (link = "log"),
               data = ss_brownfield_age_mot)

summary(mod6srm)
drop1(mod6srm, test = "Chi")

chkres(mod6srm)


plot(log(ss_brownfield_age_mot$chao) ~ ss_brownfield_age_mot$TSCY)

# and observed

mod6som <- glm(Species ~ TSCY,
               family = poisson (link = "log"),
               data = ss_brownfield_age_mot)

summary(mod6som)
drop1(mod6som, test = "Chi")

chkres(mod6som)

plot(log(ss_brownfield_age_mot$Species) ~ ss_brownfield_age_mot$TSCY)


# next, completeness (just to check)

mod6cm <-  glm(Completeness ~ TSCY,
               data = ss_brownfield_age_mot)

summary(mod6cm)
drop1(mod6cm, test = "Chi")

chkres(mod6cm)

plot(ss_brownfield_age_mot$Completeness ~ ss_brownfield_age_mot$TSCY)



# no problem here
# so, moving on to Irr

mod6im <-  glm(Irr.cjm ~ TSCY,
               family = binomial (link = "logit"),
               data = ss_brownfield_age_mot)


summary(mod6im)
drop1(mod6im, test = "Chi")

chkres(mod6im)


plot(ss_brownfield_age_mot$Irr.cjm ~ ss_brownfield_age_mot$TimeSinceClosure)



# dragonflies


data_locs_dra <- data_locs[which(data_locs$GROUP=="Dragonfly"), ]


## now we can start to pick out actual stats on these squares
# let's start very simple:
# figure out species richness (raw and extrapolated)


all.sqs.dra <- levels(droplevels(data_locs_dra$GRID_REF))

richness_dra <- data.frame()


for(x in all.sqs.dra){
  sq.dat <- data_locs_dra[which(data_locs_dra$GRID_REF == x), ]
  print(x)
  
  # first crack raw spp richness
  sq.dat$PRESENT <- 1
  
  sq.spp <- ddply(sq.dat, .(COMMON_NAME, BMSCODE), summarise,
                  RECS = sum(PRESENT),
                  TOT = sum(COUNT),
                  MEAN = mean(COUNT))
  
  # now extrapolate spp richness across samples
  
  sq.pres <- sq.dat[,c(4,6,12)]
  
  # wrangle the data to be in wide form (one row per sample)
  
  sq.mat <- dcast(sq.pres, VISIT_DATE ~ COMMON_NAME,
                  value.var = "COUNT",
                  fun.aggregate = sum)
  
  # get rid of the sample ID column, otherwise it gets treated as an extra species!
  sq.mat <- sq.mat[,-1]
  
  out <- specpool(sq.mat)
  out$Square <- x
  
  richness_dra <- rbind(richness_dra, out)
  
}

summary(richness_dra)



#### figure out species rarity 
# use R package "Rarity"


### need a full site-level matrix (presence/absence, one row per square)

# trim out only the columns we need
data_locs_mat_dra <- data_locs_dra[,c(2,6,12)]
data_locs_mat_dra$presence <- 1


data_mat_dra <- dcast(data_locs_mat_dra, GRID_REF ~ COMMON_NAME,
                      value.var = "presence",
                      fun.aggregate = sum)

# make rownames from col 1

rownames(data_mat_dra) <- data_mat_dra[,1]

data_mat_dra <- data_mat_dra[,-1]

# finally transpose the matrix so that assemblages are in columns and species in rows

data_mat_dra_t <- t(data_mat_dra)


### and a list of species with number of square-level occurrences

# first narrow down to one row per species*square
spp_occ_sq_dra <- ddply(data_locs_mat_dra, .(GRID_REF, COMMON_NAME), summarise,
                        presence = mean(presence))

summary(spp_occ_sq_dra)

# then summarise the data for each species
spp_occ_dra <- ddply(spp_occ_sq_dra, .(COMMON_NAME), summarise,
                     squares = sum(presence))

# need to rejig the structure of this slightly

spp_occ_mat_dra <- spp_occ_dra$squares
names(spp_occ_mat_dra) <- spp_occ_dra$COMMON_NAME
head(spp_occ_mat_dra)

### then need to apply the functions rWeights and Irr

rarity.weights_dra <- rWeights(spp_occ_mat_dra)

square_rarity_dra <- data.frame(Irr(assemblages = data_mat_dra_t, W = rarity.weights_dra))
hist(square_rarity_dra$Irr)

# this actually causes some problems further down the line, because this index generates a lot of species with weight zero
# so I'm going to try and generate my own simpler version for comparison...
# there are 2096 squares in total, so each species can be observed in 1-2096 squares
# for my index I'll take 1-(obs.sqs/2096)

rarity.weights_dra$CJM <- 1-(rarity.weights_dra$Q/2096)

# normalise these between 0 and 1 (slightly boosting the rare species and taking away from the common)
minW <- min(rarity.weights_dra$CJM)
maxW <- max(rarity.weights_dra$CJM)

rarity.weights_dra$CJMz <- (rarity.weights_dra$CJM - minW)/(maxW - minW)

summary(rarity.weights_dra$CJMz)

rarity.weights.cjm_dra <- rarity.weights_dra[,c(1,2,6,4)]
colnames(rarity.weights.cjm_dra)[3] <- "W"

square_rarity.cjm_dra <- data.frame(Irr(assemblages = data_mat_dra_t, W = rarity.weights.cjm_dra))
hist(square_rarity.cjm_dra$Irr)


### now merge all these indices back together so that we have one row per square in a single dataframe
square_rarity_dra$Square <- rownames(square_rarity_dra)
square_rarity.cjm_dra$Square <- rownames(square_rarity.cjm_dra)

colnames(square_rarity.cjm_dra) <- c("Irr.cjm", "Richness.cjm", "Square")

square_rr_dra <- merge(square_rarity_dra, square_rarity.cjm_dra)

square_stats_dra <- merge(richness_dra, square_rr_dra)

# restore squares with zero richness, and fill in blanks
square_stats_complete_dra <- merge(all.sqs, square_stats_dra, all = T)




### now we are ready to start some basic analyses
# first pick out the stats for the squares in pairs

# several ways to tackle this
# I'll want a long-form version (still one row for every square, dra labelling which group they belong to)
# and possibly a wide-form version (one row per brownfield square, columns for )

# first, add a factor label to each square

ss_label_dra <- data.frame()

y <- 0

for (x in squares){
  print(x)
  y <- y+1
  
  sq.pair <- pairs[which(pairs$SQUARE == x), ]
  
  sq.sqs <- append(x, levels(droplevels(sq.pair$NEIGHBOUR)))
  
  sq.dat <- square_stats_complete_dra[which(square_stats_complete_dra$Square %in% sq.sqs), ]
  
  
  sq.pair.base <- sq.pair[,c('NEIGHBOUR','TYPE')]
  sq.pair.base$NEIGHBOUR <- as.character(sq.pair.base$NEIGHBOUR)
  sq.pair.base[5,] <- c(x, "Target")
  colnames(sq.pair.base) <- c("Square", "Type")
  
  
  sq.dat$Group <- y
  
  sq.type <- merge(sq.pair.base, sq.dat)
  
  
  ss_label_dra <- rbind(ss_label_dra, sq.type)
  
  
}

summary(ss_label_dra)

ss_label_dra$Type <- as.factor(ss_label_dra$Type)
ss_label_dra$Group <- as.factor(ss_label_dra$Group)

ss_label_dra$Completeness <- ss_label_dra$Species/ss_label_dra$chao

summary(ss_label_dra)


# tack in the area of brownfield site

ss_landfill_dra <- merge(WCBS_sum, ss_label_dra, all.y = T)


ss_landfill_dra$LANDFILL_AREA[is.na(ss_landfill_dra$LANDFILL_AREA)] <- 0
ss_landfill_dra$PERC[is.na(ss_landfill_dra$PERC)] <- 0


# write this out for easy figure-plotting

write.csv(ss_landfill_dra, "../../PreppedData/all_dragons.csv")


### we are finally ready to begin interrogating the data!

# first a very simple comparison

ss_landfill_dra$Type <- relevel(ss_landfill_dra$Type, ref = "Target")

matched.dat <- ss_landfill_dra[which(ss_landfill_dra$Type != "Neighbouring"), ]
grouped.dat <- ss_landfill_dra[which(ss_landfill_dra$Type != "Matched"), ]


# matched

hist(matched.dat$Species)

mod1rd <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = matched.dat)

summary(mod1rd)
drop1(mod1rd, test = "Chi")

chkres(mod1rd, matched.dat$Type)


# grouped 

mod2rd <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = grouped.dat)

summary(mod2rd)
drop1(mod2rd, test = "Chi")

chkres(mod2rd, grouped.dat$Type)


# combined

mod3rd <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = ss_landfill_dra)

summary(mod3rd)
drop1(mod3rd, test = "Chi")

chkres(mod3rd)

plot(ss_landfill_dra$Species ~ ss_landfill_dra$Type)


## repeat with estimated richness

# matched

hist(matched.dat$chao)

mod1ed <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = matched.dat)

summary(mod1ed)
drop1(mod1ed, test = "Chi")

chkres(mod1ed, matched.dat$Type)


# grouped 

mod2ed <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = grouped.dat)

summary(mod2ed)
drop1(mod2ed, test = "Chi")

chkres(mod2ed, grouped.dat$Type)


# combined

mod3ed <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = ss_landfill_dra)

summary(mod3ed)
drop1(mod3ed, test = "Chi")

chkres(mod3ed)

plot(ss_landfill_dra$chao ~ ss_landfill_dra$Type)

# this also suggests that target sites have HIGHER SR than matched habitat, dra in this case possibly LOWER than surrounding landscape

# some sign of the bias here that was present in the plant data, so
# let's check sampling completeness out

hist(ss_landfill_dra$Completeness)

# matched

mod1cd <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
               data = matched.dat)

summary(mod1cd)
drop1(mod1cd, test = "Chi")

chkres(mod1cd, matched.dat$Type)



# grouped

mod2cd <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
               data = grouped.dat)

summary(mod2cd)
drop1(mod2cd, test = "Chi")

chkres(mod2cd, grouped.dat$Type)



# combined

mod3cd <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
               data = ss_landfill_dra)

summary(mod3cd)
drop1(mod3cd, test = "Chi")

chkres(mod3cd)

plot(ss_landfill_dra$Completeness ~ ss_landfill_dra$Type)

# no cause for concern here
# maybe slightly lower (but not significant) in neighbouring sites



## species richness is significantly lower in brownfield-type habitats than other habitat types, 
# but significantly higher in brownfield squares than non-brownfields of same type

# in other words, the habitats associated with brownfield sites have lower species richness than other local habitat types
# but brownfield sites have higher SR than other sites with similar habitat...?


# repeat for Irr

# matched

hist(matched.dat$Irr.cjm)

mod1id <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data = matched.dat)

summary(mod1id)
drop1(mod1id, test = "Chi")

chkres(mod1id, matched.dat$Type)


# grouped

mod2id <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data=grouped.dat)

summary(mod2id)
drop1(mod2id, test = "Chi")

chkres(mod2id)


# combined

mod3id <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data=ss_landfill_dra)

summary(mod3id)
drop1(mod3id, test = "Chi")

chkres(mod3id)

plot(ss_landfill_dra$Irr.cjm ~ ss_landfill_dra$Type)


# like plants, no sign of a rarity effect


## the next set of questions relates to the effect of *area of brownfield* within target sites
# we can also try including non-target sites with 0% < x < 5% brownfield

ss_brownfield_dra <- ss_landfill_dra[which(ss_landfill_dra$PERC > 5), ]


# first, estimated SR

mod4srd <- glm(round(chao) ~ LANDFILL_AREA,
               family = poisson (link = "log"),
               data = ss_brownfield_dra)

summary(mod4srd)
drop1(mod4srd, test = "Chi")

chkres(mod4srd)


plot(log(ss_brownfield_dra$chao) ~ log(ss_brownfield_dra$LANDFILL_AREA))


# next, completeness (just to check)

mod4cd <-  glm(Completeness ~ LANDFILL_AREA,
               data = ss_brownfield_dra)

summary(mod4cd)
drop1(mod4cd, test = "Chi")

chkres(mod4cd)

plot(ss_brownfield_dra$Completeness ~ log(ss_brownfield_dra$LANDFILL_AREA))

# no problem here


mod4id <-  glm(Irr.cjm ~ LANDFILL_AREA,
               family = binomial (link = "logit"),
               data = ss_brownfield_dra)


summary(mod4id)
drop1(mod4id, test = "Chi")

chkres(mod4id)


plot(ss_brownfield_dra$Irr.cjm ~ log(ss_brownfield_dra$LANDFILL_AREA))


## for completeness, do these analyses with the percentage landfill by area (instead of raw area)

# first, estimated SR

mod5srd <- glm(round(chao) ~ log(PERC),
               family = poisson (link = "log"),
               data = ss_brownfield_dra)

summary(mod5srd)
drop1(mod5srd, test = "Chi")

chkres(mod5srd, variable1 = log(ss_brownfield_dra$PERC))

plot(log(ss_brownfield_dra$chao) ~ log(ss_brownfield_dra$PERC))


# and observed


mod5sod <- glm(Species ~ log(PERC),
               family = poisson (link = "log"),
               data = ss_brownfield_dra)

summary(mod5sod)
drop1(mod5sod, test = "Chi")

chkres(mod5sod, variable1 = log(ss_brownfield_dra$PERC))

plot(log(ss_brownfield_dra$Species) ~ log(ss_brownfield_dra$PERC))



# next, completeness (just to check)

mod5cd <-  glm(Completeness ~ log(PERC),
               data = ss_brownfield_dra)

summary(mod5cd)
drop1(mod5cd, test = "Chi")

chkres(mod5cd, variable1 = log(ss_brownfield_dra$PERC))

plot(ss_brownfield_dra$Completeness ~ log(ss_brownfield_dra$PERC))


# unlike last time, this appears to be fine
# so, moving on to Irr

mod5id <-  glm(Irr.cjm ~ log(PERC),
               family = binomial (link = "logit"),
               data = ss_brownfield_dra)


summary(mod5id)
drop1(mod5id, test = "Chi")

chkres(mod5id, log(ss_brownfield_dra$PERC))


plot(ss_brownfield_dra$Irr.cjm ~ log(ss_brownfield_dra$PERC))



## and now we have a set of questions relating to the age of the brownfield sites

ss_brownfield_age_dra <- ss_brownfield_dra[which(!is.na(ss_brownfield_dra$TimeSinceClosure)), ]


# first, estimated SR

mod6srd <- glm(round(chao) ~ TSCY,
               family = poisson (link = "log"),
               data = ss_brownfield_age_dra)

summary(mod6srd)
drop1(mod6srd, test = "Chi")

chkres(mod6srd)


plot(log(ss_brownfield_age_dra$chao) ~ ss_brownfield_age_dra$TSCY)

# and observed

mod6sod <- glm(Species ~ TSCY,
               family = poisson (link = "log"),
               data = ss_brownfield_age_dra)

summary(mod6sod)
drop1(mod6sod, test = "Chi")

chkres(mod6sod)

plot(log(ss_brownfield_age_dra$Species) ~ ss_brownfield_age_dra$TSCY)


# next, completeness (just to check)

mod6cd <-  glm(Completeness ~ TSCY,
               data = ss_brownfield_age_dra)

summary(mod6cd)
drop1(mod6cd, test = "Chi")

chkres(mod6cd)

plot(ss_brownfield_age_dra$Completeness ~ ss_brownfield_age_dra$TSCY)



# no problem here
# so, moving on to Irr

mod6id <-  glm(Irr.cjm ~ TSCY,
               family = binomial (link = "logit"),
               data = ss_brownfield_age_dra)


summary(mod6id)
drop1(mod6id, test = "Chi")

chkres(mod6id)


plot(ss_brownfield_age_dra$Irr.cjm ~ ss_brownfield_age_dra$TimeSinceClosure)









# damselflies


data_locs_dam <- data_locs[which(data_locs$GROUP=="Damselfly"), ]


## now we can start to pick out actual stats on these squares
# let's start very simple:
# figure out species richness (raw and extrapolated)


all.sqs.dam <- levels(droplevels(data_locs_dam$GRID_REF))

richness_dam <- data.frame()


for(x in all.sqs.dam){
  sq.dat <- data_locs_dam[which(data_locs_dam$GRID_REF == x), ]
  print(x)
  
  # first crack raw spp richness
  sq.dat$PRESENT <- 1
  
  sq.spp <- ddply(sq.dat, .(COMMON_NAME, BMSCODE), summarise,
                  RECS = sum(PRESENT),
                  TOT = sum(COUNT),
                  MEAN = mean(COUNT))
  
  # now extrapolate spp richness across samples
  
  sq.pres <- sq.dat[,c(4,6,12)]
  
  # wrangle the data to be in wide form (one row per sample)
  
  sq.mat <- dcast(sq.pres, VISIT_DATE ~ COMMON_NAME,
                  value.var = "COUNT",
                  fun.aggregate = sum)
  
  # get rid of the sample ID column, otherwise it gets treated as an extra species!
  sq.mat <- sq.mat[,-1]
  
  out <- specpool(sq.mat)
  out$Square <- x
  
  richness_dam <- rbind(richness_dam, out)
  
}

summary(richness_dam)



#### figure out species rarity 
# use R package "Rarity"


### need a full site-level matrix (presence/absence, one row per square)

# trim out only the columns we need
data_locs_mat_dam <- data_locs_dam[,c(2,6,12)]
data_locs_mat_dam$presence <- 1


data_mat_dam <- dcast(data_locs_mat_dam, GRID_REF ~ COMMON_NAME,
                      value.var = "presence",
                      fun.aggregate = sum)

# make rownames from col 1

rownames(data_mat_dam) <- data_mat_dam[,1]

data_mat_dam <- data_mat_dam[,-1]

# finally transpose the matrix so that assemblages are in columns and species in rows

data_mat_dam_t <- t(data_mat_dam)


### and a list of species with number of square-level occurrences

# first narrow down to one row per species*square
spp_occ_sq_dam <- ddply(data_locs_mat_dam, .(GRID_REF, COMMON_NAME), summarise,
                        presence = mean(presence))

summary(spp_occ_sq_dam)

# then summarise the data for each species
spp_occ_dam <- ddply(spp_occ_sq_dam, .(COMMON_NAME), summarise,
                     squares = sum(presence))

# need to rejig the structure of this slightly

spp_occ_mat_dam <- spp_occ_dam$squares
names(spp_occ_mat_dam) <- spp_occ_dam$COMMON_NAME
head(spp_occ_mat_dam)

### then need to apply the functions rWeights and Irr

rarity.weights_dam <- rWeights(spp_occ_mat_dam)

square_rarity_dam <- data.frame(Irr(assemblages = data_mat_dam_t, W = rarity.weights_dam))
hist(square_rarity_dam$Irr)

# this actually causes some problems further down the line, because this index generates a lot of species with weight zero
# so I'm going to try and generate my own simpler version for comparison...
# there are 2096 squares in total, so each species can be observed in 1-2096 squares
# for my index I'll take 1-(obs.sqs/2096)

rarity.weights_dam$CJM <- 1-(rarity.weights_dam$Q/2096)

# normalise these between 0 and 1 (slightly boosting the rare species and taking away from the common)
minW <- min(rarity.weights_dam$CJM)
maxW <- max(rarity.weights_dam$CJM)

rarity.weights_dam$CJMz <- (rarity.weights_dam$CJM - minW)/(maxW - minW)

summary(rarity.weights_dam$CJMz)

rarity.weights.cjm_dam <- rarity.weights_dam[,c(1,2,6,4)]
colnames(rarity.weights.cjm_dam)[3] <- "W"

square_rarity.cjm_dam <- data.frame(Irr(assemblages = data_mat_dam_t, W = rarity.weights.cjm_dam))
hist(square_rarity.cjm_dam$Irr)


### now merge all these indices back together so that we have one row per square in a single dataframe
square_rarity_dam$Square <- rownames(square_rarity_dam)
square_rarity.cjm_dam$Square <- rownames(square_rarity.cjm_dam)

colnames(square_rarity.cjm_dam) <- c("Irr.cjm", "Richness.cjm", "Square")

square_rr_dam <- merge(square_rarity_dam, square_rarity.cjm_dam)

square_stats_dam <- merge(richness_dam, square_rr_dam)

# restore squares with zero richness, and fill in blanks
square_stats_complete_dam <- merge(all.sqs, square_stats_dam, all = T)




### now we are ready to start some basic analyses
# first pick out the stats for the squares in pairs

# several ways to tackle this
# I'll want a long-form version (still one row for every square, dam labelling which group they belong to)
# and possibly a wide-form version (one row per brownfield square, columns for )

# first, add a factor label to each square

ss_label_dam <- data.frame()

y <- 0

for (x in squares){
  print(x)
  y <- y+1
  
  sq.pair <- pairs[which(pairs$SQUARE == x), ]
  
  sq.sqs <- append(x, levels(droplevels(sq.pair$NEIGHBOUR)))
  
  sq.dat <- square_stats_complete_dam[which(square_stats_complete_dam$Square %in% sq.sqs), ]
  
  
  sq.pair.base <- sq.pair[,c('NEIGHBOUR','TYPE')]
  sq.pair.base$NEIGHBOUR <- as.character(sq.pair.base$NEIGHBOUR)
  sq.pair.base[5,] <- c(x, "Target")
  colnames(sq.pair.base) <- c("Square", "Type")
  
  
  sq.dat$Group <- y
  
  sq.type <- merge(sq.pair.base, sq.dat)
  
  
  ss_label_dam <- rbind(ss_label_dam, sq.type)
  
  
}

summary(ss_label_dam)

ss_label_dam$Type <- as.factor(ss_label_dam$Type)
ss_label_dam$Group <- as.factor(ss_label_dam$Group)

ss_label_dam$Completeness <- ss_label_dam$Species/ss_label_dam$chao

summary(ss_label_dam)


# tack in the area of brownfield site

ss_landfill_dam <- merge(WCBS_sum, ss_label_dam, all.y = T)


ss_landfill_dam$LANDFILL_AREA[is.na(ss_landfill_dam$LANDFILL_AREA)] <- 0
ss_landfill_dam$PERC[is.na(ss_landfill_dam$PERC)] <- 0


# write this out for easy figure-plotting

write.csv(ss_landfill_dam, "../../PreppedData/all_damsels.csv")


### we are finally ready to begin interrogating the data!

# first a very simple comparison

ss_landfill_dam$Type <- relevel(ss_landfill_dam$Type, ref = "Target")

matched.dat <- ss_landfill_dam[which(ss_landfill_dam$Type != "Neighbouring"), ]
grouped.dat <- ss_landfill_dam[which(ss_landfill_dam$Type != "Matched"), ]


# matched

hist(matched.dat$Species)

mod1rl <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = matched.dat)

summary(mod1rl)
drop1(mod1rl, test = "Chi")

chkres(mod1rl, matched.dat$Type)


# grouped 

mod2rl <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = grouped.dat)

summary(mod2rl)
drop1(mod2rl, test = "Chi")

chkres(mod2rl, grouped.dat$Type)


# combined

mod3rl <- glmer(Species ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = ss_landfill_dam)

summary(mod3rl)
drop1(mod3rl, test = "Chi")

chkres(mod3rl)
chkconv(mod3rl)

plot(ss_landfill_dam$Species ~ ss_landfill_dam$Type)
# appears that target sites have higher SR than either other category...


## repeat with estimated richness

# matched

hist(matched.dat$chao)

mod1el <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = matched.dat)

summary(mod1el)
drop1(mod1el, test = "Chi")

chkres(mod1el, matched.dat$Type)


# grouped 

mod2el <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = grouped.dat)

summary(mod2el)
drop1(mod2el, test = "Chi")

chkres(mod2el, grouped.dat$Type)


# combined

mod3el <- glmer(round(chao) ~ Type + (1|Group),
                family = poisson (link = "log"),
                data = ss_landfill_dam)

summary(mod3el)
drop1(mod3el, test = "Chi")

chkres(mod3el)

plot(ss_landfill_dam$chao ~ ss_landfill_dam$Type)

# this also suggests that target sites have HIGHER SR than matched habitat, dam in this case possibly LOWER than surrounding landscape

# some sign of the bias here that was present in the plant data, so
# let's check sampling completeness out

hist(ss_landfill_dam$Completeness)

# matched

mod1cl <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
               data = matched.dat)

summary(mod1cl)
drop1(mod1cl, test = "Chi")

chkres(mod1cl, matched.dat$Type)



# grouped

mod2cl <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
               data = grouped.dat)

summary(mod2cl)
drop1(mod2cl, test = "Chi")

chkres(mod2cl, grouped.dat$Type)



# combined

mod3cl <- glmer(Completeness ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
               data = ss_landfill_dam)

summary(mod3cl)
drop1(mod3cl, test = "Chi")

chkres(mod3cl)

plot(ss_landfill_dam$Completeness ~ ss_landfill_dam$Type)

# no major cause for concern here
# maybe slightly lower (but not significant) in neighbouring sites



## species richness is significantly lower in brownfield-type habitats than other habitat types, 
# but significantly higher in brownfield squares than non-brownfields of same type

# in other words, the habitats associated with brownfield sites have lower species richness than other local habitat types
# but brownfield sites have higher SR than other sites with similar habitat...?


# repeat for Irr

# matched

hist(matched.dat$Irr.cjm)

mod1il <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data = matched.dat)

summary(mod1il)
drop1(mod1il, test = "Chi")

chkres(mod1il, matched.dat$Type)


# grouped

mod2il <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data=grouped.dat)

summary(mod2il)
drop1(mod2il, test = "Chi")

chkres(mod2il)


# combined

mod3il <- glmer(Irr.cjm ~ Type + (1|Group),
                family = binomial (link = "logit"),
                weights = n,
                data=ss_landfill_dam)

summary(mod3il)
drop1(mod3il, test = "Chi")

chkres(mod3il)

plot(ss_landfill_dam$Irr.cjm ~ ss_landfill_dam$Type)


# like plants, no sign of a rarity effect


## the next set of questions relates to the effect of *area of brownfield* within target sites

ss_brownfield_dam <- ss_landfill_dam[which(ss_landfill_dam$PERC > 5), ]


# first, estimated SR

mod4srl <- glm(round(chao) ~ LANDFILL_AREA,
               family = poisson (link = "log"),
               data = ss_brownfield_dam)

summary(mod4srl)
drop1(mod4srl, test = "Chi")

chkres(mod4srl)


plot(log(ss_brownfield_dam$chao) ~ log(ss_brownfield_dam$LANDFILL_AREA))


# next, completeness (just to check)

mod4cl <-  glm(Completeness ~ LANDFILL_AREA,
               data = ss_brownfield_dam)

summary(mod4cl)
drop1(mod4cl, test = "Chi")

chkres(mod4cl)

plot(ss_brownfield_dam$Completeness ~ log(ss_brownfield_dam$LANDFILL_AREA))

# no problem here


mod4il <-  glm(Irr.cjm ~ LANDFILL_AREA,
               family = binomial (link = "logit"),
               data = ss_brownfield_dam)


summary(mod4il)
drop1(mod4il, test = "Chi")

chkres(mod4il)


plot(ss_brownfield_dam$Irr.cjm ~ log(ss_brownfield_dam$LANDFILL_AREA))


## for completeness, do these analyses with the percentage landfill by area (instead of raw area)

# first, estimated SR

mod5srl <- glm(round(chao) ~ log(PERC),
               family = poisson (link = "log"),
               data = ss_brownfield_dam)

summary(mod5srl)
drop1(mod5srl, test = "Chi")

chkres(mod5srl)

plot(log(ss_brownfield_dam$chao) ~ log(ss_brownfield_dam$PERC))

# and observed

mod5sol <- glm(Species ~ log(PERC),
               family = poisson (link = "log"),
               data = ss_brownfield_dam)

summary(mod5sol)
drop1(mod5sol, test = "Chi")

chkres(mod5sol)

plot(log(ss_brownfield_dam$chao) ~ log(ss_brownfield_dam$PERC))


# next, completeness (just to check)

mod5cl <-  glm(Completeness ~ log(PERC),
               data = ss_brownfield_dam)

summary(mod5cl)
drop1(mod5cl, test = "Chi")

chkres(mod5cl)

plot(ss_brownfield_dam$Completeness ~ log(ss_brownfield_dam$PERC))


# unlike last time, this appears to be fine
# so, moving on to Irr

mod5il <-  glm(Irr.cjm ~ log(PERC),
               family = binomial (link = "logit"),
               data = ss_brownfield_dam)


summary(mod5il)
drop1(mod5il, test = "Chi")

chkres(mod5il)


plot(ss_brownfield_dam$Irr.cjm ~ log(ss_brownfield_dam$PERC))



## and now we have a set of questions relating to the age of the brownfield sites

ss_brownfield_age_dam <- ss_brownfield_dam[which(!is.na(ss_brownfield_dam$TimeSinceClosure)), ]


# first, estimated SR

mod6srl <- glm(round(chao) ~ TSCY,
               family = poisson (link = "log"),
               data = ss_brownfield_age_dam)

summary(mod6srl)
drop1(mod6srl, test = "Chi")

chkres(mod6srl)


plot(log(ss_brownfield_age_dam$chao) ~ ss_brownfield_age_dam$TSCY)

# and observed

mod6sol <- glm(Species ~ TSCY,
               family = poisson (link = "log"),
               data = ss_brownfield_age_dam)

summary(mod6sol)
drop1(mod6sol, test = "Chi")

chkres(mod6sol)

plot(log(ss_brownfield_age_dam$Species) ~ ss_brownfield_age_dam$TSCY)


# next, completeness (just to check)

mod6cl <-  glm(Completeness ~ TSCY,
               data = ss_brownfield_age_dam)

summary(mod6cl)
drop1(mod6cl, test = "Chi")

chkres(mod6cl)

plot(ss_brownfield_age_dam$Completeness ~ ss_brownfield_age_dam$TSCY)



# no problem here
# so, moving on to Irr

mod6il <-  glm(Irr.cjm ~ TSCY,
               family = binomial (link = "logit"),
               data = ss_brownfield_age_dam)


summary(mod6il)
drop1(mod6il, test = "Chi")

chkres(mod6il)


plot(ss_brownfield_age_dam$Irr.cjm ~ ss_brownfield_age_dam$TimeSinceClosure)






### this script stops here -
# it's not possible to do any within-square analyses with the WCBS data

