########################################################
####   Script for constructing cross-taxon figures  ####
########################################################

##### script setup ####

### Clear the current workspace (DOESN'T CLEAR LOADED LIBRARIES)
rm(list=ls())

### install if necessary and then load the libraries you need

j <- c("rstudioapi","plyr","raster","dplyr","vegan","reshape2","Rarity","lme4","MASS","glmmTMB","effects","ggplot2","ggstance","scales","gridExtra","ggpubr")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

### set working directory to script's saved location

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## read in some functions we might need

source("CheckResidsFunction.R")
source("CheckConvergenceFunction.R")

### read in the datasets

birds <- read.csv("PreppedData/all_birds.csv")
plants <- read.csv("PreppedData/all_plants.csv")
wcb <- read.csv("PreppedData/all_wc_butterflies.csv")
hsb <- read.csv("PreppedData/all_hs_butterflies.csv")
moths <- read.csv("PreppedData/all_moths.csv")
dragons <- read.csv("PreppedData/all_dragons.csv")
damsels <- read.csv("PreppedData/all_damsels.csv")


# create a list of datasets that can be looped over

datasets <- list(birds,plants,wcb,hsb,moths,dragons,damsels)

# and for convenience a list of those taxa

dataset_names <- c("Birds","Plants","Wider countryside butterflies","Habitat specialist butterflies","Moths","Dragonflies","Damselflies")

## loop over the datasets - for each one, run the same analyses and extract predictions, effect sizes etc. for figure plotting

categorical_effects <- data.frame(stringsAsFactors = F)
categorical_estimates <- data.frame(stringsAsFactors = F)

continuous_effects <- data.frame(stringsAsFactors = F)
continuous_estimates <- data.frame(stringsAsFactors = F)

age_effects <- data.frame(stringsAsFactors = F)
age_estimates <- data.frame(stringsAsFactors = F)


for (x in 1:7){
  data <- datasets[[x]]
  name <- dataset_names[[x]]
  
  print(name)
  
  
  # first build the models we need
  # then extract (a) effect sizes and (b) predictions from each
  data$Type <- relevel(as.factor(data$Type), ref = "Target")
  data$Group <- as.factor(data$Group)
  
  ### observed richness
  mod1r <- glmer(Species ~ Type + (1|Group),
                 family = poisson (link = "log"),
                 data = data)
  
  # effect sizes
  
  matched <- cbind(name,"Observed richness", "Matched", as.numeric(summary(mod1r)$coefficients[[2,1]]), as.numeric(summary(mod1r)$coefficients[[2,2]]))
  neighbouring <- cbind(name,"Observed richness", "Neighbouring", as.numeric(summary(mod1r)$coefficients[[3,1]]), as.numeric(summary(mod1r)$coefficients[[3,2]]))
  
  out <- rbind(matched,neighbouring)
  
  categorical_effects <- rbind(categorical_effects, out)

  
  # estimates
  
  errors <- data.frame(effect(c("Type"),mod1r))
  errors$Taxon <- name
  errors$Response <- "Observed richness"
  
  categorical_estimates <- rbind(categorical_estimates, errors)
  
  
  ### estimated richness
  # we need to set up an ifelse because we want to treat moths differently (due to the underdispersion)
  
  if (name == "Moths"){
    # moths 
    mod1e <- glmmPQL(chao ~ Type,
                     random = list(~1|Group),
                     family = quasipoisson (link = "log"),
                     data = data)
    
    # effect sizes
    matched <- cbind(name,"Estimated richness", "Matched", as.numeric(summary(mod1e)$tTable[[2,1]]), as.numeric(summary(mod1e)$tTable[[2,2]]))
    neighbouring <- cbind(name,"Estimated richness", "Neighbouring", as.numeric(summary(mod1e)$tTable[[3,1]]), as.numeric(summary(mod1e)$tTable[[3,2]]))
    
    out <- rbind(matched,neighbouring)
    
    categorical_effects <- rbind(categorical_effects, out)
    
    # estimates
    
    errors <- data.frame(effect(c("Type"),mod1e))
    
    errors$Taxon <- name
    errors$Response <- "Estimated richness"
    
    categorical_estimates <- rbind(categorical_estimates, errors)
    
  } else {
  
    # everything else
    mod1e <- glmer(round(chao) ~ Type + (1|Group),
                 family = poisson (link = "log"),
                 data = data)
  
    # effect sizes
    matched <- cbind(name,"Estimated richness", "Matched", as.numeric(summary(mod1e)$coefficients[[2,1]]), as.numeric(summary(mod1e)$coefficients[[2,2]]))
    neighbouring <- cbind(name,"Estimated richness", "Neighbouring", as.numeric(summary(mod1e)$coefficients[[3,1]]), as.numeric(summary(mod1e)$coefficients[[3,2]]))
  
    out <- rbind(matched,neighbouring)
  
    categorical_effects <- rbind(categorical_effects, out)
  
    # estimates
  
    errors <- data.frame(effect(c("Type"),mod1e))
    errors$Taxon <- name
    errors$Response <- "Estimated richness"
  
    categorical_estimates <- rbind(categorical_estimates, errors)
  }
  
  
  
  ### sampling completeness
  
  mod1c <- glmer(Completeness ~ Type + (1|Group),
                 family = binomial (link = "logit"),
                 weights = n,
                data = data)
  
  summary(mod1c)
  
  # effect sizes
  matched <- cbind(name,"Completeness", "Matched", as.numeric(summary(mod1c)$coefficients[[2,1]]), as.numeric(summary(mod1c)$coefficients[[2,2]]))
  neighbouring <- cbind(name,"Completeness", "Neighbouring", as.numeric(summary(mod1c)$coefficients[[3,1]]), as.numeric(summary(mod1c)$coefficients[[3,2]]))
  
  out <- rbind(matched,neighbouring)
  
  categorical_effects <- rbind(categorical_effects, out)
  
  # estimates
  
  errors <- data.frame(effect(c("Type"),mod1c))
  errors$Taxon <- name
  errors$Response <- "Completeness"
  
  categorical_estimates <- rbind(categorical_estimates, errors)
  
  
  
  
  ### Irr
  
  mod1i <- glmer(Irr.cjm ~ Type + (1|Group),
                 family = binomial (link = "logit"),
                 weights = n,
                 data=data)
  
  summary(mod1i)
  
  # effect sizes
  matched <- cbind(name,"Rarity index", "Matched", as.numeric(summary(mod1i)$coefficients[[2,1]]), as.numeric(summary(mod1i)$coefficients[[2,2]]))
  neighbouring <- cbind(name,"Rarity index", "Neighbouring", as.numeric(summary(mod1i)$coefficients[[3,1]]), as.numeric(summary(mod1i)$coefficients[[3,2]]))
  
  out <- rbind(matched,neighbouring)
  
  categorical_effects <- rbind(categorical_effects, out)
  
  # estimates
  
  errors <- data.frame(effect(c("Type"),mod1i))
  errors$Taxon <- name
  errors$Response <- "Rarity index"
  
  categorical_estimates <- rbind(categorical_estimates, errors)
  
  
  
  
  ##### the next set of questions relates to the effect of *area of brownfield* within target sites

  data_bf <- data[which(data$PERC > 5), ]
  
  
  ### first, estimated SR
  
  mod2sr <- glm(round(chao) ~ log(PERC),
                family = poisson (link = "log"),
                data = data_bf)
  
  
  # effect sizes
  out <- cbind(name,"Estimated richness", 
               as.numeric(summary(mod2sr)$coefficients[[2,1]]), as.numeric(summary(mod2sr)$coefficients[[2,2]]))
  
  continuous_effects <- rbind(continuous_effects, out)
  
  # estimates

  
  newdata <- expand.grid(PERC = c(seq(max(round(min(data_bf$PERC),2),0.002),round(max(data_bf$PERC),2),0.1)), chao = 0)
  preds <- predict(mod2sr, newdata = newdata, type="response", se.fit = T)
  
  newdata$chao <- preds[[1]]
  newdata$se <- preds[[2]]
  
  newdata$Taxon <- name
  newdata$Variable <- "Estimated richness"
  
  continuous_estimates <- rbind(continuous_estimates, newdata)
  
  
  
  
  
  ### and observed
  
  mod2so <- glm(Species ~ log(PERC),
                family = poisson (link = "log"),
                data = data_bf)

  
  # effect sizes
  out <- cbind(name,"Observed richness", 
               as.numeric(summary(mod2so)$coefficients[[2,1]]), as.numeric(summary(mod2so)$coefficients[[2,2]]))
  
  continuous_effects <- rbind(continuous_effects, out)
  
  # estimates
  
  newdata <- expand.grid(PERC = seq(max(round(min(data_bf$PERC),2),0.002),round(max(data_bf$PERC),2),0.1), chao = 0)
  preds <- predict(mod2so, newdata = newdata, type="response", se.fit = T)
  
  newdata$chao <- preds[[1]]
  newdata$se <- preds[[2]]
  
  newdata$Taxon <- name
  newdata$Variable <- "Observed richness"
  
  continuous_estimates <- rbind(continuous_estimates, newdata)
  
  
  
  ### Irr
  
  mod2i <-  glm(Irr.cjm ~ log(PERC),
                family = binomial (link = "logit"),
                data = data_bf)
  
  
  # effect sizes
  out <- cbind(name,"Rarity index", 
               as.numeric(summary(mod2i)$coefficients[[2,1]]), as.numeric(summary(mod2i)$coefficients[[2,2]]))
  
  continuous_effects <- rbind(continuous_effects, out)
  
  # estimates
  
  newdata <- expand.grid(PERC = seq(max(round(min(data_bf$PERC),2),0.002),round(max(data_bf$PERC),2),0.1), chao = 0)
  preds <- predict(mod2i, newdata = newdata, type="response", se.fit = T)
  
  newdata$chao <- preds[[1]]
  newdata$se <- preds[[2]]
  
  newdata$Taxon <- name
  newdata$Variable <- "Rarity index"
  
  continuous_estimates <- rbind(continuous_estimates, newdata)
  
  ####### the last set of questions relates to the effect of *age of brownfield* within target sites
  
  data_bf_age <- data_bf[which(!is.na(data_bf$TSCY)), ]
  
  
  ### first, estimated SR
  
  mod3sr <- glm(round(chao) ~ TSCY,
                family = poisson (link = "log"),
                data = data_bf_age)
  
  
  # effect sizes
  out <- cbind(name,"Estimated richness", 
               as.numeric(summary(mod3sr)$coefficients[[2,1]]), as.numeric(summary(mod3sr)$coefficients[[2,2]]))
  
  age_effects <- rbind(age_effects, out)
  
  # estimates
  
  
  newdata <- expand.grid(TSCY = c(seq(max(round(min(data_bf_age$TSCY),2),0.002),round(max(data_bf_age$TSCY),2),0.1)), chao = 0)
  preds <- predict(mod3sr, newdata = newdata, type="response", se.fit = T)
  
  newdata$chao <- preds[[1]]
  newdata$se <- preds[[2]]
  
  newdata$Taxon <- name
  newdata$Variable <- "Estimated richness"
  
  age_estimates <- rbind(age_estimates, newdata)
  
  
  
  
  
  ### and observed
  
  mod3so <- glm(Species ~ TSCY,
                family = poisson (link = "log"),
                data = data_bf_age)
  
  
  # effect sizes
  out <- cbind(name,"Observed richness", 
               as.numeric(summary(mod3so)$coefficients[[2,1]]), as.numeric(summary(mod3so)$coefficients[[2,2]]))
  
  age_effects <- rbind(age_effects, out)
  
  # estimates
  
  newdata <- expand.grid(TSCY = c(seq(max(round(min(data_bf_age$TSCY),2),0.002),round(max(data_bf_age$TSCY),2),0.1)), chao = 0)
  preds <- predict(mod3so, newdata = newdata, type="response", se.fit = T)
  
  newdata$chao <- preds[[1]]
  newdata$se <- preds[[2]]
  
  newdata$Taxon <- name
  newdata$Variable <- "Observed richness"
  
  age_estimates <- rbind(age_estimates, newdata)
  
  
  
  ### Irr
  
  mod3i <-  glm(Irr.cjm ~ TSCY,
                family = binomial (link = "logit"),
                data = data_bf_age)
  
  
  # effect sizes
  out <- cbind(name,"Rarity index", 
               as.numeric(summary(mod3i)$coefficients[[2,1]]), as.numeric(summary(mod3i)$coefficients[[2,2]]))
  
  age_effects <- rbind(age_effects, out)
  
  # estimates
  
  newdata <- expand.grid(TSCY = c(seq(max(round(min(data_bf_age$TSCY),2),0.002),round(max(data_bf_age$TSCY),2),0.1)), chao = 0)
  preds <- predict(mod3i, newdata = newdata, type="response", se.fit = T)
  
  newdata$chao <- preds[[1]]
  newdata$se <- preds[[2]]
  
  newdata$Taxon <- name
  newdata$Variable <- "Rarity index"
  
  age_estimates <- rbind(age_estimates, newdata)
  
  ##
  
}


colnames(categorical_effects) <- c("Taxon","Variable","Type","Effect","SE")

# function to convert multiple columns from factor to numeric, and from character to factor

Funcall <- function(f, ...) f(...)
translate.cols <- function(f){
  Reduce(Funcall, list(as.numeric, as.character), f, right = TRUE)
}



categorical_effects[,4:5] <- sapply( categorical_effects[,4:5], translate.cols)
categorical_effects$Taxon <- as.factor(categorical_effects$Taxon)
categorical_effects$Variable <- as.factor(categorical_effects$Variable)
categorical_effects$Type <- as.factor(categorical_effects$Type)

summary(categorical_effects)

categorical_estimates$Taxon <- as.factor(categorical_estimates$Taxon)
categorical_estimates$Response <- as.factor(categorical_estimates$Response)

summary(categorical_estimates)

colnames(continuous_effects) <- c("Taxon","Variable","Effect","SE")
continuous_effects[,3:4] <- sapply( continuous_effects[,3:4], translate.cols)

continuous_effects$Taxon <- as.factor(continuous_effects$Taxon)
continuous_effects$Variable <- as.factor(continuous_effects$Variable)
summary(continuous_effects)

continuous_estimates$lower <- continuous_estimates$chao - 2*continuous_estimates$se
continuous_estimates$upper <- continuous_estimates$chao + 2*continuous_estimates$se

continuous_estimates$Taxon <- as.factor(continuous_estimates$Taxon)
continuous_estimates$Variable <- as.factor(continuous_estimates$Variable)

summary(continuous_estimates)


colnames(age_effects) <- c("Taxon", "Variable","Effect","SE")
age_effects[,3:4] <- sapply( age_effects[,3:4], translate.cols)
age_effects$Taxon <- as.factor(age_effects$Taxon)
age_effects$Variable <- as.factor(age_effects$Variable)
summary(age_effects)

age_estimates$lower <- age_estimates$chao - 2*age_estimates$se
age_estimates$upper <- age_estimates$chao + 2*age_estimates$se

age_estimates$Taxon <- as.factor(age_estimates$Taxon)
age_estimates$Variable <- as.factor(age_estimates$Variable)

summary(age_estimates)


## now we can start to construct some figures!

# first, let's do a sample forest plot for cross-taxon effects

categorical_effects$Min <- categorical_effects$Effect - 2*categorical_effects$SE
categorical_effects$Max <- categorical_effects$Effect + 2*categorical_effects$SE

categorical_effects$Taxon <- factor(categorical_effects$Taxon, levels = rev(c("Birds","Plants","Wider countryside butterflies",
                                                                             "Habitat specialist butterflies","Moths",
                                                                             "Dragonflies","Damselflies")))

categorical_effects$Variable <- factor(categorical_effects$Variable, levels = c("Observed richness","Estimated richness","Completeness","Rarity index"))


# rename comparisons

categorical_effects$Type[which(categorical_effects$Type == "Neighbouring")] <- "Different land-use"
categorical_effects$Type[which(categorical_effects$Type == "Matched")] <- "Matched land-use"

categorical_effects$Type <- factor(categorical_effects$Type, levels = c("Matched land-use","Different land-use"))


g1 <- ggplot()+
  geom_point(data = categorical_effects[-c(29:30),],
             aes(x = Effect, y = Taxon, group = Type, colour = Type),
             position = position_dodgev(height = -0.5))+
  geom_errorbarh(data = categorical_effects[-c(29:30),],
                 aes(xmin = Min, xmax = Max, y = Taxon, group = Type, colour = Type),
                 position = position_dodgev(height = -0.5), height = 0)+
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(strip.background = element_rect(fill="white", linetype = "blank"))+
  labs(x = "Effect size", colour = "Comparison with")+
  scale_color_manual(values = c("royalblue","darkgoldenrod"))+
  facet_wrap(~Variable, ncol = 4, scales = "free_x")

g1

ggsave("Figures/BrownfieldYN.svg", plot = g1, device = svg, width = 250, height = 100, units = "mm", limitsize = F)

# and for continuous/age effects?


continuous_effects$Min <- continuous_effects$Effect - 2*continuous_effects$SE
continuous_effects$Max <- continuous_effects$Effect + 2*continuous_effects$SE

continuous_effects$Taxon <- factor(continuous_effects$Taxon, levels = rev(c("Birds","Plants","Wider countryside butterflies",
                                                                              "Habitat specialist butterflies","Moths",
                                                                              "Dragonflies","Damselflies")))

continuous_effects$Variable <- factor(continuous_effects$Variable, levels = c("Observed richness","Estimated richness","Rarity index"))



age_effects$Min <- age_effects$Effect - 2*age_effects$SE
age_effects$Max <- age_effects$Effect + 2*age_effects$SE

age_effects$Taxon <- factor(age_effects$Taxon, levels = rev(c("Birds","Plants","Wider countryside butterflies",
                                                                            "Habitat specialist butterflies","Moths",
                                                                            "Dragonflies","Damselflies")))

age_effects$Variable <- factor(age_effects$Variable, levels = c("Observed richness","Estimated richness","Rarity index"))


# bind these two together, first labelling them

continuous_effects$Type <- "% cover of brownfield"
age_effects$Type <- "Time since closure"

# and removing HS butterflies from the continuous effects since they obscure the figure
continuous_effects_trim <- continuous_effects[which(continuous_effects$Taxon != "Habitat specialist butterflies"), ]


ca_effects <- rbind(continuous_effects_trim, age_effects)


# now plot the figure
g2 <- ggplot()+
  geom_point(data = ca_effects,
             aes(x = Effect, y = Taxon, group = Type, colour = Type),
             position = position_dodgev(height = -0.5))+
  geom_errorbarh(data = ca_effects,
                 aes(xmin = Min, xmax = Max, y = Taxon, group = Type, colour = Type),
                 position = position_dodgev(height = -0.5), height = 0)+
  geom_vline(xintercept = 0)+
  theme_bw()+
  theme(strip.background = element_rect(fill="white", linetype = "blank"))+
  scale_color_manual(values = c("royalblue","darkgoldenrod"))+
  labs(x = "Effect size", colour = "Explanatory variable")+
  facet_wrap(~Variable, ncol = 4, scales = "free_x")

g2

ggsave("Figures/BrownfieldPerc.svg", plot = g2, device = svg, width = 250, height = 100, units = "mm", limitsize = F)




##### taxon-specific plots
# probably for the SI, we'll want one multiplot per taxon summarising the estimates rather than the effects, 
# ideally with underlying data also shown 

percs <- data.frame()
ages <- data.frame()

categorical_estimates$Type <- as.character(categorical_estimates$Type)

categorical_estimates$Type[which(categorical_estimates$Type == "Neighbouring")] <- "Different land-use"
categorical_estimates$Type[which(categorical_estimates$Type == "Matched")] <- "Matched land-use"

categorical_estimates$Type <- factor(categorical_estimates$Type, levels = c("Matched land-use","Different land-use","Target"))


# build the loop

for (x in 1:7){

taxon <- dataset_names[[x]]

data <- categorical_estimates[which(categorical_estimates$Taxon == taxon), ]
full_dat <- datasets[[x]]

full_dat$Type[which(full_dat$Type == "Matched")] <- "Matched land-use"
full_dat$Type[which(full_dat$Type == "Neighbouring")] <- "Different land-use"

full_dat$Type <- factor(full_dat$Type, levels = c("Target","Matched land-use","Different land-use"))

# set up axis limits for SR plots
max <- max(full_dat$chao[!is.na(full_dat$chao)])

max_age <- max(full_dat[which('PERC' > 5),'TSCY'])

# start plotting

gp1 <- ggplot(data[which(data$Response == "Observed richness"), ])+
  geom_jitter(data = full_dat, aes(x = Type, y = Richness), colour = "grey50", alpha = 0.5, width = 0.2, height = 0)+
  geom_point(aes(x = Type, y = fit), size = 3)+
  geom_errorbar(aes(x = Type, ymin = lower, ymax = upper), width = 0.2)+
  theme_bw()+
  ylim(0,max)+
  labs(tag = "a", y = "Observed richness")


gp1



gp2 <- ggplot(data[which(data$Response == "Estimated richness"), ])+
  geom_jitter(data = full_dat, aes(x = Type, y = chao), colour = "grey50", alpha = 0.5, width = 0.2, height = 0)+
  geom_point(aes(x = Type, y = fit), size = 3)+
  geom_errorbar(aes(x = Type, ymin = lower, ymax = upper), width = 0.2)+
  theme_bw()+
  ylim(0,max)+
  labs(tag = "b", y = "Estimated richness")


gp2


# for completeness, we don't want error bars on the habitat specialist butterflies since we couldn't fit a model

if (taxon == "Habitat specialist butterflies"){
  
  gp3 <- ggplot(data[which(data$Response == "Completeness"), ])+
    geom_jitter(data = full_dat, aes(x = Type, y = Completeness), colour = "grey50", alpha = 0.5, width = 0.2, height = 0)+
    scale_y_continuous(limits = c(0,1), oob=squish)+
    theme_bw()+
    labs(tag = "c", y = "Completeness")
  
  
  gp3  
} else {
  
  gp3 <- ggplot(data[which(data$Response == "Completeness"), ])+
    geom_jitter(data = full_dat, aes(x = Type, y = Completeness), colour = "grey50", alpha = 0.5, width = 0.2, height = 0)+
    geom_point(aes(x = Type, y = fit), size = 3)+
    geom_errorbar(aes(x = Type, ymin = lower, ymax = upper), width = 0.2)+
    scale_y_continuous(limits = c(0,1), oob=squish)+
    theme_bw()+
    labs(tag = "c", y = "Completeness")
  
  
  gp3
}



gp4 <- ggplot(data[which(data$Response == "Rarity index"), ])+
  geom_jitter(data = full_dat, aes(x = Type, y = Irr.cjm), colour = "grey50", alpha = 0.5, width = 0.2, height = 0)+
  geom_point(aes(x = Type, y = fit), size = 3)+
  geom_errorbar(aes(x = Type, ymin = lower, ymax = upper), width = 0.2)+
  scale_y_continuous(limits = c(0,1), oob=squish)+
  theme_bw()+
  labs(tag = "d", y = "Rarity index")


gp4


mp1 <- grid.arrange(gp1,gp2,gp3,gp4)


## and now try to plot the continuous effects as line graphs in the same way
# first for area

continuous_estimates <- continuous_estimates[which(continuous_estimates$PERC > 5), ]

data_cont <- continuous_estimates[which(continuous_estimates$Taxon == taxon), ]
full_dat_cont <- full_dat[which(full_dat$PERC > 5), ]

## at this point, rip some data out to use down the way

dat_percs <- full_dat_cont[,c('Square','PERC')]
dat_percs$TAXON <- as.factor(taxon)
dat_percs$DATASET <- as.factor(ifelse(x == 1, "BBS",
                            ifelse(x == 2, "NPMS",
                                   "WCBS")))

percs <- rbind(percs,dat_percs)

# back to plotting

gp5 <- ggplot()+
  geom_point(data = full_dat_cont,
             aes(x = PERC, y = Species), colour = "grey50", alpha = 0.5)+
  geom_line(data = data_cont[which(data_cont$Variable == "Observed richness"), ],
            aes(x = PERC, y = chao))+
  geom_line(data = data_cont[which(data_cont$Variable == "Observed richness"), ],
            aes(x = PERC, y = lower),
            linetype = "dashed")+
  geom_line(data = data_cont[which(data_cont$Variable == "Observed richness"), ],
            aes(x = PERC, y = upper),
            linetype = "dashed")+
  scale_x_log10(breaks = c(5,10,20,50))+
  theme_bw()+
  ylim(0,max)+
  labs(tag = "e", y = "Observed richness", x = "% cover of brownfield")

gp5

gp6 <- ggplot()+
  geom_point(data = full_dat_cont,
             aes(x = PERC, y = chao), colour = "grey50", alpha = 0.5)+
  geom_line(data = data_cont[which(data_cont$Variable == "Estimated richness"), ],
            aes(x = PERC, y = chao))+
  geom_line(data = data_cont[which(data_cont$Variable == "Estimated richness"), ],
            aes(x = PERC, y = lower),
            linetype = "dashed")+
  geom_line(data = data_cont[which(data_cont$Variable == "Estimated richness"), ],
            aes(x = PERC, y = upper),
            linetype = "dashed")+
  scale_x_log10(breaks = c(5,10,20,50))+
  theme_bw()+
  ylim(0,max)+
  labs(tag = "g", y = "Estimated richness", x = "% cover of brownfield")

gp6



gp7 <- ggplot()+
  geom_point(data = full_dat_cont,
             aes(x = PERC, y = Irr.cjm), colour = "grey50", alpha = 0.5)+
  geom_line(data = data_cont[which(data_cont$Variable == "Rarity index"), ],
            aes(x = PERC, y = chao))+
  geom_line(data = data_cont[which(data_cont$Variable == "Rarity index"), ],
            aes(x = PERC, y = lower),
            linetype = "dashed")+
  geom_line(data = data_cont[which(data_cont$Variable == "Rarity index"), ],
            aes(x = PERC, y = upper),
            linetype = "dashed")+
  labs(tag = "i", y = "Rarity index", x = "% cover of brownfield")+
  scale_x_log10(breaks = c(5,10,20,50))+
  ylim(0,1)+
  theme_bw()

gp7


# then for age

data_cont <- age_estimates[which(age_estimates$Taxon == taxon), ]
full_dat_cont <- full_dat[which(full_dat$PERC > 5), ]

## at this point, rip some data out to use down the way

dat_ages <- full_dat_cont[,c('Square','TSCY')]
dat_ages$TAXON <- as.factor(taxon)
dat_ages$DATASET <- as.factor(ifelse(x == 1, "BBS",
                                      ifelse(x == 2, "NPMS",
                                             "WCBS")))

ages <- rbind(ages,dat_ages)

# back to plotting

gp8 <- ggplot()+
  geom_point(data = full_dat_cont,
             aes(x = TSCY, y = Species), colour = "grey50", alpha = 0.5)+
  geom_line(data = data_cont[which(data_cont$Variable == "Observed richness"), ],
            aes(x = TSCY, y = chao))+
  geom_line(data = data_cont[which(data_cont$Variable == "Observed richness"), ],
            aes(x = TSCY, y = lower),
            linetype = "dashed")+
  geom_line(data = data_cont[which(data_cont$Variable == "Observed richness"), ],
            aes(x = TSCY, y = upper),
            linetype = "dashed")+
  scale_x_continuous(limits = c(0,max_age),breaks = c(0,20,40,60,80))+
  theme_bw()+
  ylim(0,max)+
  labs(tag = "f", y = "Observed richness", x = "Time since closure (years)")

gp8


gp9 <- ggplot()+
  geom_point(data = full_dat_cont,
             aes(x = TSCY, y = chao), colour = "grey50", alpha = 0.5)+
  geom_line(data = data_cont[which(data_cont$Variable == "Estimated richness"), ],
            aes(x = TSCY, y = chao))+
  geom_line(data = data_cont[which(data_cont$Variable == "Estimated richness"), ],
            aes(x = TSCY, y = lower),
            linetype = "dashed")+
  geom_line(data = data_cont[which(data_cont$Variable == "Estimated richness"), ],
            aes(x = TSCY, y = upper),
            linetype = "dashed")+
  scale_x_continuous(limits = c(0,max_age),breaks = c(0,20,40,60,80))+
  theme_bw()+
  ylim(0,max)+
  labs(tag = "h", y = "Estimated richness", x = "Time since closure (years)")

gp9



gp10 <- ggplot()+
  geom_point(data = full_dat_cont,
             aes(x = TSCY, y = Irr.cjm), colour = "grey50", alpha = 0.5)+
  geom_line(data = data_cont[which(data_cont$Variable == "Rarity index"), ],
            aes(x = TSCY, y = chao))+
  geom_line(data = data_cont[which(data_cont$Variable == "Rarity index"), ],
            aes(x = TSCY, y = lower),
            linetype = "dashed")+
  geom_line(data = data_cont[which(data_cont$Variable == "Rarity index"), ],
            aes(x = TSCY, y = upper),
            linetype = "dashed")+
  labs(tag = "j", y = "Rarity index", x = "Time since closure (years)")+
  scale_x_continuous(limits = c(0,max_age),breaks = c(0,20,40,60,80))+
  ylim(0,1)+
  theme_bw()

gp10


mp2 <- grid.arrange(gp5, gp8, gp6, gp9, gp7, gp10, ncol = 2)

mp3 <- grid.arrange(mp1, mp2)

ggsave(filename = paste0("Figures/Taxon plots/",taxon,".svg"), plot = mp3, device = svg, width = 200, height = 240, units = "mm", limitsize = F)
ggsave(filename = paste0("Figures/Taxon plots/",taxon,".png"), plot = mp3, device = png, width = 20000, height = 24000, units = "mm", limitsize = F)



}


### now to do something with those percs values - 
# we want a histogram or something showing the spread of different levels of brownfield coverage across different datasets

percs$Square <- as.factor(percs$Square)
summary(percs)


# lots of squares are duplicated, sometimes multiple times, within the WCBS dataset
# and we want to whittle this away

percs_trim <- ddply(percs, .(Square,DATASET), numcolwise(mean))


# build up a plot

theme_set(theme_classic())


g3 <- ggplot(percs_trim, aes(x = PERC))+
  geom_histogram(aes(fill = DATASET),
    binwidth = 1, color = "black")+
  scale_fill_brewer(palette = "Paired")+
  scale_x_continuous(breaks = c(5,10,20,30,40,50,60,70))+
  labs(x = "% cover of brownfield", y = "No. target squares", fill = "Dataset")


g3

ggsave(filename = "Figures/Percs.svg", plot = g3, device = svg, width = 200, height = 120, units = "mm", limitsize = F)


### and the same with the ages values - 

ages$Square <- as.factor(ages$Square)
summary(ages)


# lots of squares are duplicated, sometimes multiple times, within the WCBS dataset
# and we want to whittle this away

ages_trim <- ddply(ages, .(Square,DATASET), numcolwise(mean))


# build up a plot

theme_set(theme_classic())


g4 <- ggplot(ages_trim, aes(x = TSCY))+
  geom_histogram(aes(fill = DATASET),
                 binwidth = 1, color = "black")+
  scale_fill_brewer(palette = "Paired")+
  scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80))+
  labs(x = "Time since closure (years)", y = "No. target squares", fill = "Dataset")


g4

ggsave(filename = "Figures/Ages.svg", plot = g4, device = svg, width = 200, height = 120, units = "mm", limitsize = F)
