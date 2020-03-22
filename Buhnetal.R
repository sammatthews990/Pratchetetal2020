
##### LOAD PACKAGES #####

library(dplyr)
# install.packages("chron")
# install.packages("nnet")
# install.packages("rstan")
# install.packages("brms")
# install.packages("vegan")
# install.packages("MASS")
# install.packages("purrr")
# install.packages("lme4")
# install.packages("export")
library(nnet)
library(MASS)
library(vegan)
library(chron)
library(lme4)
library(brms)
library(export)
library(RColorBrewer)
library(tidyverse)
library(emmeans)
library(lsmeans)
library(brms)
library(tidybayes)
library(brmstools)
library(gridExtra)


##### LOAD DATA #####

move = read.csv("COTSmovement.csv", header = T)
levels(move$Behaviour)
move = move %>% 
  filter(!Behaviour %in% c("", "(s)spawnng")) %>%
  mutate(Date=as.Date(Date, format="%d/%m/%Y"),
         Time=chron::chron(times=as.character(Time)),
         Depth=as.integer(gsub("m", "",gsub(" ", "",Depth))),
         Diam_cm=as.integer(as.character(Diam_cm)),
         DateTime = as.POSIXct(paste(Date,Time), format="%Y-%m-%d %H:%M:%S"),
         Behaviour = factor(Behaviour),
         Time_Cat = factor(Time_Cat, levels = c("Morning", "Midday", "Afternoon", "Night")),
         Behaviour = relevel(Behaviour, ref="Resting"))

##### HELPER FUNCTIONS #####

# Creates Effects Plot of BRMS Model 
ploteffects.brms = function(model, title=NULL, plotlabels=NULL, exp=FALSE) {
  summy = as.data.frame(summary(model)$fixed)
  summy = summy %>% 
    mutate(Variable = factor(rownames(summy))) %>%
    filter(!Variable %in% "Intercept") %>%
    mutate(Sig = factor(ifelse(sign(`l-95% CI`)==sign(`u-95% CI`), "Y", "N"), levels = c("Y", "N")))
  if(exp==TRUE){
    summy[1:4] = exp(summy[1:4])
  }
  ggplot(summy,aes(x=Variable,y=Estimate))  + 
    geom_hline(yintercept = 0) +
    geom_linerange(aes(ymin=`l-95% CI`,ymax=`u-95% CI`)) +
    geom_point(size=3, aes(colour=Sig)) + 
    scale_color_manual(breaks = c("Y", "N"), values = c("black", "grey")) +
    guides(colour=FALSE) +
    ggtitle(title) +
    labs(x=plotlabels[2], y=plotlabels[1]) +
    coord_flip() 
}

# Saves all marginal effects plots of BRMS model
marginaleffects.brms = function(model, Effects) {
  plots = list()
  for (i in 1: length(Effects)){
    plots[i] = plot(marginal_effects(model, effects=Effects[i]))
  }
  return(plots)
}

# Generates pvalue from BRMS model
mcmcpvalue <- function(samp) {
  ## elementary version that creates an empirical p-value for the
  ## hypothesis that the columns of samp have mean zero versus a general
  ## multivariate distribution with elliptical contours.
  
  ## differences from the mean standardized by the observed
  ## variance-covariance factor
  
  ## Note, I put in the bit for single terms
  if (length(dim(samp)) == 0) {
    std <- backsolve(chol(var(samp)), cbind(0, t(samp)) - mean(samp),
                     transpose = TRUE)
    sqdist <- colSums(std * std)
    sum(sqdist[-1] > sqdist[1])/length(samp)
  } else {
    std <- backsolve(chol(var(samp)), cbind(0, t(samp)) - colMeans(samp),
                     transpose = TRUE)
    sqdist <- colSums(std * std)
    sum(sqdist[-1] > sqdist[1])/nrow(samp)
  }
  
}





##### SUMMARY TABLES #####

# Convert Depth to Habitat
move$Habitat = ifelse(move$Depth == 2, "Flat", 
                      ifelse(move$Depth %in% c(3,8), "Crest", "Slope"))

# COTS Summary - To plot and model Density "COTSPer100"
data.COTS = move %>% 
  dplyr::group_by(Reef, Location, Date, Time_Cat, Transect, Depth, Coral.Cover) %>%
  summarise(TotalCOTS = n(),
            MeanSize = mean(Diam_cm, na.rm=T)) %>%
  mutate(Habitat = ifelse(Depth == 2, "Flat", 
                          ifelse(Depth %in% c(3,8), "Crest", "Slope")),
         # Create COTS density column
         COTSPer100 = ifelse(Reef %in% "Gili Lankanfushi", TotalCOTS/2,
                             ifelse(Depth == 6, TotalCOTS/5, TotalCOTS/2.5)))

# COTS SUmmary Size by Habitat
data.COTS.HAB = move %>% 
  dplyr::group_by(Reef,Habitat) %>%
  summarise(TotalCOTS = n(),
            MeanSize = mean(Diam_cm, na.rm=T),
            sesize = sd(Diam_cm, na.rm=T)/sqrt(n()))

# COTS SUmmary Size by Time
data.COTS.TIME = move %>% 
  dplyr::group_by(Reef,Time_Cat) %>%
  summarise(TotalCOTS = n(),
            MeanSize = mean(Diam_cm, na.rm=T),
            sesize = sd(Diam_cm, na.rm=T)/sqrt(n()))

# COTS SUmmary Behaviour for Model
data.COTS.beh = move %>% filter(!Habitat %in% "Flat") %>%
  dplyr::group_by(Reef,Behaviour, Transect) %>%
  summarise(n = n()) %>%
  left_join(move %>% 
              dplyr::group_by(Reef,Transect) %>%
              summarise(nTot = n())) %>%
  mutate(Freq = n/nTot)


# COTS SUmmary Behaviour for summary stats (i.e Mean proportion)
data.COTS.beh2 = move %>% filter(!Habitat %in% "Flat") %>%
  dplyr::group_by(Reef,Behaviour, Transect) %>%
  summarise(n = n()) %>%
  left_join(move %>% 
              dplyr::group_by(Reef,Transect) %>%
              summarise(nTot = n())) %>%
  mutate(Freq = n/nTot) %>%
  ungroup() %>% group_by(Reef, Behaviour)%>%
  summarise(MeanFreq = mean(Freq),
            SEFreq = (sd(Freq)/sqrt(n())))
# Summary behaviour by Reef and time category
data.COTS.beh3 = move %>% filter(!Habitat %in% "Flat") %>%
  dplyr::group_by(Reef,Behaviour, Time_Cat, Transect) %>%
  summarise(n = n()) %>%
  left_join(move %>% 
              dplyr::group_by(Reef,Time_Cat,Transect) %>%
              summarise(nTot = n())) %>%
  mutate(Freq = n/nTot)



# COTS SUmmary Density by Habitat
data.COTS.sum = data.COTS %>% group_by(Reef, Habitat) %>%
  summarise(nTransects = n(),
            MeanCOTS = mean(COTSPer100),
            SECOTS = sd(COTSPer100)/sqrt(n()))

# # COTS SUmmary Density by Reef
data.COTS.ReefSum = data.COTS %>% filter(!Habitat %in% "Flat") %>% 
  group_by(Reef) %>%
  summarise(nTransects = n(),
            MeanCOTS = mean(COTSPer100),
            SECOTS = sd(COTSPer100)/sqrt(n()))

#### COTS PLOTS ####

# Density by Habitat
p1 = ggplot(data.COTS %>% filter(!Habitat %in% "Flat"), aes(x=Habitat,y=COTSPer100, colour=Reef)) + geom_boxplot() + theme_bw(base_size = 14) + facet_wrap(~Reef, ncol = 2) +
  theme(legend.position = "none") +
  xlab("Habitat Type") + ylab(bquote('COTS Density (per 100'~m^2~')')) +
  scale_color_manual(values = c("orange","seagreen"))

# temporary dataframe for mean size
df = move %>% filter(!Habitat %in% "Flat") %>% group_by(Reef) %>% 
  summarise(mean = mean(Diam_cm, na.rm=T),
            median = median(Diam_cm, na.rm=T),
            se = sd(Diam_cm, na.rm = T)/sqrt(n()))

# Size Histogram
p2 = ggplot(data = move %>% left_join(df), aes(x=Diam_cm, fill=Reef)) + geom_histogram(bins=7, color="white") + facet_wrap(~Reef, ncol = 2) + 
  theme_bw(base_size = 14)+ 
  xlab("Diameter (cm)") + ylab('Number of starfish') +
  geom_vline(data=filter(move, Reef=="Rib Reef"), 
             aes(xintercept=31), colour="black", linetype = "dashed") +
  geom_vline(data=filter(move, Reef=="Gili Lankanfushi"), 
             aes(xintercept=40.89), colour="black", linetype = "dashed") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("orange","seagreen"))

# Diamter by Time
p3  = ggplot(move %>% filter(!Habitat %in% "Flat"), aes(x=Time_Cat,y=Diam_cm, colour=Reef)) + geom_boxplot() + theme_bw(base_size = 14) + facet_wrap(~Reef, ncol = 2) +
  theme(legend.position = "none") +
  xlab("Observation Time") + ylab("Diameter (cm)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("orange","seagreen"))

# Diameter by Habitat
p4 = ggplot(move %>% filter(!Habitat %in% "Flat"), aes(x=Habitat,y=Diam_cm, colour=Reef)) + geom_boxplot() + theme_bw(base_size = 14) + facet_wrap(~Reef, ncol = 2) +
  theme(legend.position = "none") +
  xlab("Habitat Type") + ylab("Diameter (cm)") +
  scale_color_manual(values = c("orange","seagreen"))


# Behaviour Proportionss by Reef
p5 = ggplot(data.COTS.beh, aes(x=Behaviour,y=Freq, colour=Reef)) + geom_boxplot() + theme_bw(base_size = 14) + facet_wrap(~Reef, ncol = 2) +
  theme(legend.position = "none") +
  xlab("Habitat Type") + ylab("Diameter (cm)") +
  scale_color_manual(values = c("orange","seagreen"))

# Arrange Plots
gridExtra::grid.arrange(p1 + ggtitle("(a)"),p2 + ggtitle("(b)"), p4 + ggtitle("(c)"), p3 + ggtitle("(d)"))
# Save Plots
graph2eps(file="Figure1.eps", aspectr=3, font = "Arial", height = 8, width=8, bg = "transparent")
graph2tif(file="Figure1.tif",dpi=300, font = "Arial", height=8, width=8)

#### COTS Density and Size Models ####

# DENSITY MODEL ----
lm1 = lm(COTSPer100~Reef + Habitat + Time_Cat, data = data.COTS %>% filter(!Habitat %in% "Flat"))
# Density Full Facotial
lm2 = lm(COTSPer100~Reef*Habitat  + Reef*Time_Cat, data = data.COTS %>% filter(!Habitat %in% "Flat"))

# Compare models
MuMIn::AICc(lm1, lm2)
# lm1 wins
# ANOVA
anova(lm1)
# PostHOC Pairwise comparisons
lsmeans(lm1, list(pairwise ~ Time_Cat|Reef))
lsmeans(lm1, list(pairwise ~ Habitat+Reef))

# DIAMETER MODEL ----
lm3 = lm(Diam_cm~Reef + Time_Cat + Habitat, data = move %>% filter(!Habitat %in% "Flat"))
# Diameter FUll Factorial
lm4 = lm(Diam_cm~Reef*Time_Cat + Reef*Habitat, data = move %>% filter(!Habitat %in% "Flat"))
# Compare models
MuMIn::AICc(lm3, lm4)
# lm4 wins
# ANOVA
anova(lm4)
# PostHOC Pairwise comparisons
lsmeans(lm4, list(pairwise ~ Habitat|Reef))
lsmeans(lm4, list(pairwise ~ Time_Cat|Reef))


#### CORAL COMMUNITY COMPOSITION ANALYSIS ####


#### Prepare Data ####

data.CORAL = read.csv("Maldives_PIT.csv", na.strings = "") %>% 
  mutate_at(.funs = as.character,.vars = c("subtrate", "genus", "growth.form"))
data.CORAL$genus[which(is.na(data.CORAL$genus))] =  data.CORAL$subtrate[which(is.na(data.CORAL$genus))]

data.CORAL.GBR = read.csv("COTS_PITRIB.csv") %>% dplyr::select(-c(6:15), -Date, -Location) %>%
  dplyr::rename(Montipora = "Montipora.ecnructing",
                Pocillopora = "Pocillopora.spp.")

# Captitalize all first words
data.CORAL$genus = tools::toTitleCase(data.CORAL$genus)
data.CORAL$growth.form = tools::toTitleCase(data.CORAL$growth.form)

unique(data.CORAL$genus)
unique(data.CORAL$growth.form)

# Create Proportions
data.PROP = data.CORAL %>%
  group_by(site, transect, depth) %>%
  mutate(Points = n()) %>%
  group_by(site, transect, depth, Points, genus) %>%
  summarise(n = n()) %>%
  mutate(Prop = n/Points)

data.PROP.wide = data.PROP %>% dplyr::select(-Prop) %>% spread(genus, n)

# MAtch maldives data to GBR data
data.CORAL.MAL = data.PROP.wide %>% dplyr::select(-c(Algae:Anemone), -CCA, -`Dead Coral`, -Halimeda,
                                                  -`Macro Algae`, -c(Rock:Sand), -Sponge) %>%
  mutate(Other.HC = Agaricia + Astreopora + Cyphastrea + Diploastrea + Favia + Galaxea + Goniastrea + Merulina + Other + Pavona + Symphyllia + Tubastrea) %>%
  ungroup() %>%
  dplyr::select(`site`, `transect`, `depth`,Acropora, Montipora, Pocillopora, Porites, Other.HC, Sarcophyton, -Points) %>%
  rename(Soft.coral = "Sarcophyton",
         Reef = "site",
         Transect = "transect", Depth = "depth") %>% mutate(Depth = as.factor(Depth))

# Create coral data set for composition analysis
data.CORAL = rbind(data.CORAL.GBR, data.CORAL.MAL)
data.CORAL[is.na(data.CORAL)]=0

# Summary for Pie Chart
data.CORAL.sum = data.CORAL %>% 
  mutate(NonCoral = (100-rowSums(.[4:9]))) %>%
  gather(key = "Coral", value = "PercCover", -c(Reef:Depth)) %>%
  group_by(Reef, Coral) %>%
  summarise(Mean = mean(PercCover),
            SE = sd(PercCover)/sqrt(n()))

#### Coral Composition Pie Chart ####
data.CORAL.sum$Reef = factor(data.CORAL.sum$Reef, levels = c("Rib Reef","Lankanfushi island"), labels = c("Rib Reef","Gili Lankanfushi"))
data.CORAL.sum$Coral = factor(data.CORAL.sum$Coral, 
                              levels = c("NonCoral", "Acropora", "Porites", "Pocillopora", "Montipora", "Soft.coral", "Other.HC"), labels = c("Non Coral", "Acropora", "Porites", "Pocillopora", "Montipora", "Soft Coral", "Other Hard Coral"))
p6 = ggplot(data=data.CORAL.sum, aes(x="", y=Mean, fill = Coral)) + 
  geom_bar(width = 1, stat = "identity", color = "white") +
  facet_wrap(~Reef) + 
  
  # geom_text(aes(label = paste0(round(Mean), "%")), 
  #             position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("#a9a9a9", brewer.pal(6, "Spectral"))) +
  coord_polar(theta="y") +
  # ylab("") + 
  xlab("") +
  theme_bw(base_size = 14) + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(), 
        axis.title.y=element_text(size = 30), 
        axis.title.x=element_blank(),
        axis.ticks = element_blank()) +
  labs(fill="Benthic Category") + ggtitle("(a)")
p6


#### NMDS COmpositoin Plot ####

library(vegan)
library(MASS)

# Dissimilarity Matrix
CORAL.dis <- vegdist(data.CORAL[,-c(1:3)], "bray")
# NMDS
CORAL.nmds <- metaMDS(data.CORAL[,-c(1:3)])

# Analysis of Similarity - ANOSIM
CORAL.ano <- anosim(CORAL.dis, data.CORAL$Reef)
# Print pVal
CORAL.ano$signif
CORAL.ano$statistic


# Get Scores to build in GGplot
data.scores <- as.data.frame(scores(CORAL.nmds))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$site <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$grp <- data.CORAL$Reef  #  add the grp variable created earlier
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(CORAL.nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)

# Assign groups to Reef to plot convex hull
grp.a <- data.scores[data.scores$grp == "Rib Reef", ][chull(data.scores[data.scores$grp == 
                                                                          "Rib Reef", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data.scores[data.scores$grp == "Lankanfushi island", ][chull(data.scores[data.scores$grp == 
                                                                                    "Lankanfushi island", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b

# CHange name of other "Species"
species.scores$species[5:6] = c("Other Hard Coral", "Soft Coral")

p7 =ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.30) + # add the convex hulls
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=4) + # add the point markers
  geom_label(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  scale_colour_manual(labels = c("Rib Reef", "Gili Lankanfushi"), values=c( "seagreen", "orange")) +
  scale_fill_manual(labels = c("Rib Reef", "Gili Lankanfushi"),values=c( "seagreen", "orange")) +
  scale_shape_manual(labels = c("Rib Reef", "Gili Lankanfushi"), values = 16:17) +
  coord_equal() +
  guides(fill=guide_legend(title="Reef", reverse =T)) +
  guides(color=guide_legend(title="Reef", reverse =T)) +
  guides(shape=guide_legend(title="Reef", reverse =T)) +
  # geom_rect(aes(ymin = 0.75, ymax = 0.95, xmin = 0.55, xmax = 0.99), fill= "white", color = "black") +
  annotate(geom="text", x=0.99, y=0.9, label="Stress = 0.099", hjust=1) +
  annotate(geom="text", x=0.99, y=0.83, label="R = 0.863", hjust=1) +
  annotate(geom="text", x=0.99, y=0.76, label="P = 0.001", hjust=1) +
  theme_bw(base_size=14) +
  ggtitle("(b)")+
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        # axis.title.x = element_text(size=18), # remove x-axis labels
        # axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())

# Arrange Pie and NMDS Plot
gridExtra::grid.arrange(p6, p7)
# Export Plot
graph2eps(file="Figure2.eps", aspectr=3, font = "Arial", height = 8, width=8, bg = "transparent")
graph2tif(file="Figure2.tiff",dpi=300, font = "Arial", height=8, width=8)


#### EXPOSURE BRMS MODEL ####

# Add density and coral to dataset
data.MOVE = move %>% left_join(dplyr::select(data.COTS, Transect, Habitat, COTSPer100)) %>%
  left_join(dplyr::select(data.CORAL, Transect, Acropora, Pocillopora)) %>%
  filter(!Habitat %in% "Flat") %>%
  mutate(Feed = ifelse(Behaviour %in% "Feeding", 1, 0))

# Fit Models
model.brms = brm(formula = Exposure_binary ~ Habitat + Diam_cm + Time_Cat*Reef + Behaviour + Acropora,
                 data = data.MOVE,
                 family = bernoulli(),
                 warmup = 200, iter = 1000, chains = 3, cores = 3,
                 control = list(adapt_delta = 0.95))
model.brms2 = brm(formula = Exposure_binary ~ Habitat + Diam_cm*Reef + Time_Cat*Reef + Acropora + COTSPer100 + Behaviour,
                  data = data.MOVE,
                  family = bernoulli(),
                  warmup = 200, iter = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.95))

model.brms3 = brm(formula = Exposure_binary ~ Habitat + Diam_cm*Reef + Time_Cat*Reef + Behaviour,
                  data = data.MOVE,
                  family = bernoulli(),
                  warmup = 200, iter = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.95))
model.brms4 = brm(formula = Exposure_binary ~ Habitat + Diam_cm + Time_Cat*Acropora + Behaviour,
                  data = data.MOVE,
                  family = bernoulli(),
                  warmup = 200, iter = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.95))
model.brms5 = brm(formula = Exposure_binary ~ Habitat + Diam_cm + Time_Cat*Reef + Behaviour,
                  data = data.MOVE,
                  family = bernoulli(),
                  warmup = 200, iter = 1000, chains = 3, cores = 3,
                  control = list(adapt_delta = 0.95))
model.brms6 = brm(formula = Exposure_binary ~ Habitat + Diam_cm*Time_Cat + Time_Cat*Reef + Behaviour + Acropora,
                 data = data.MOVE,
                 family = bernoulli(),
                 warmup = 200, iter = 1000, chains = 3, cores = 3,
                 control = list(adapt_delta = 0.95))


summary(model.brms)
summary(model.brms2)
summary(model.brms3)
summary(model.brms4)
summary(model.brms5)
summary(model.brms6)
# Compare
WAIC(model.brms, model.brms2, model.brms3, model.brms4, model.brms5, model.brms6)

# Model 1 wins
summy = summary(model.brms)$fixed
# Create Effects Plot
pa = ploteffects.brms(model.brms) + theme_bw(base_size = 14) + 
  scale_x_discrete(limits = rev(rownames(summy)[-1]),
                   labels = c("Night - Morning : Reef(Rib)", "Afternoon - Morning : Reef(Rib)",
                              "Midday - Morning : Reef(Rib)","% Acropora", "Moving - Resting",
                              "Feeding - Resting",
                              "Rib Reef - Gili Lankanfushi", "Night - Morning", "Afternoon - Morning",
                              "Midday - Morning", "Diameter of Starfish (cm)", "Slope - Crest")) +
  ylab("Effect Size (Log Odds-Ratio)") + ggtitle("(a)") + xlab("Variable")

# Diameter Marginal Effects plot
pb = plot(marginal_effects(model.brms, "Diam_cm"), points=F)[[1]] 
pb = pb + theme_bw(base_size = 14) + xlab("Diameter of Starfish (cm)") + 
  ylab("P(Exposure)") + ggtitle("(b)")
# Time:Reef Marginal effects plot
pc = plot(marginal_effects(model.brms, "Time_Cat:Reef"), points=F)[[1]] 
pc = pc + theme_bw(base_size = 14) + xlab("Time of Observation") + 
  ylab("P(Exposure)") + scale_color_manual(values = c( "seagreen", "orange")) +
  theme(legend.position = "left") + ggtitle("(c)")
# Behaviou Marginal effects plot
pd = plot(marginal_effects(model.brms, "Behaviour"), points=F)[[1]] 
pd = pd + theme_bw(base_size = 14) + xlab("Behaviour") + 
  ylab("P(Exposure)") + ggtitle("(d)")

 

# Arrange PLots
gridExtra::grid.arrange(pa + theme(plot.background = element_rect(size=0.5,linetype="solid",color="grey")),
                        pb + theme(plot.background = element_rect(size=0.5,linetype="solid",color="grey")),
                        pc + theme(plot.background = element_rect(size=0.5,linetype="solid",color="grey")),
                        pd + theme(plot.background = element_rect(size=0.5,linetype="solid",color="grey")),
                        ncol=2, widths = c(3,2))
# Save plots
graph2eps(file="Figure3.eps", aspectr=3, font = "Arial", height = 8, width=10, bg = "transparent")
graph2tif(file="Figure3.tiff",dpi=300, font = "Arial", height=8, width=10)

#### Post Hos Comparisons ####

#get the adjusted means
model_em <- emmeans(model.brms,  ~ Time_Cat | Reef)
model_em2 <- emmeans(model.brms,  ~  Reef | Time_Cat)
model_em3 <- emmeans(model.brms,  ~  Behaviour)


#get all possible contrasts
cont <- contrast(model_em, "tukey")
cont2 <- contrast(model_em2, "tukey")
cont3 <- contrast(model_em3, "tukey")
cont
cont2
cont3

# Pvalues
mcmcpvalue(as.matrix(model.brms)[, "b_Diam_cm"])
mcmcpvalue(as.matrix(model.brms)[, "b_Time_CatNight"])
mcmcpvalue(as.matrix(model.brms)[, "b_Acropora"])
mcmcpvalue(as.matrix(model.brms)[, "b_HabitatSlope"])

#### YAY WE DID IT! ####


#### Behaviour Model ####

# BEHAVIOURAL PROPORTIONS MODEL ----
# lm5 = lm(Freq~Reef*Behaviour, data = data.COTS.beh)

# z proportional test
table(data.MOVE$Behaviour, data.MOVE$Reef)
table(data.MOVE$Reef)
#Feeding
prop.test(c(70,288), c(159,823), alternative = "greater")
#Resting
prop.test(c(74,469), c(159,823), alternative = "less")
# Moving
prop.test(c(15,66), c(159,823), alternative = "less")



data.MOVE.MAL = data.MOVE %>% dplyr::filter(Reef %in% "Gili Lankanfushi")
data.MOVE.RIB = data.MOVE %>% dplyr::filter(Reef %in% "Rib Reef") %>%
  mutate(Size = ifelse(Diam_cm < 25, "u25", "o25"))
data.MOVE.RIB.u25 = data.MOVE %>% dplyr::filter(Reef %in% "Rib Reef" & Diam_cm <25)
data.MOVE.RIB.o25 = data.MOVE %>% dplyr::filter(Reef %in% "Rib Reef" & Diam_cm >25)

# p8 = ggplot(data.COTS.beh3 %>% filter(Behaviour %in% "Feeding"), 
#        aes(x=Time_Cat,y=Freq, colour=Reef)) + geom_boxplot() + 
#   theme_bw(base_size = 11) + facet_wrap(~Reef, ncol = 2) +
#   theme(legend.position = "none", axis.title.x = element_blank()) +
#   ylab("Proportion of COTS Feeding") +
#   scale_color_manual(values = c("seagreen","orange"))
# p9 = ggplot(data.COTS.beh3 %>% filter(Behaviour %in% "Moving"), 
#        aes(x=Time_Cat,y=Freq, colour=Reef)) + geom_boxplot() + 
#   theme_bw(base_size = 11) + facet_wrap(~Reef, ncol = 2) +
#   theme(legend.position = "none",axis.title.x = element_blank()) +
#   xlab("Observation Time") + ylab("Proportion of COTS Moving") +
#   scale_color_manual(values = c("seagreen","orange"))
# p10= ggplot(data.COTS.beh3 %>% filter(Behaviour %in% "Resting"), 
#             aes(x=Time_Cat,y=Freq, colour=Reef)) + geom_boxplot() + 
#   theme_bw(base_size = 11) + facet_wrap(~Reef, ncol = 2) +
#   theme(legend.position = "none") +
#   xlab("Observation Time") + ylab("Proportion of COTS Moving") +
#   scale_color_manual(values = c("seagreen","orange"))
# 
# grid.arrange(p8,p9,p10, ncol=1)
# # Save plots
# graph2eps(file="Figure4.eps", aspectr=3, font = "Arial", height = 8, width=6, bg = "transparent")
# graph2tif(file="Figure4.tiff",dpi=300, font = "Arial", height=8, width=6)

# Rib Feeding
prop.table(table(data.MOVE.RIB$Behaviour, data.MOVE.RIB$Time_Cat),margin = 2)
table(data.MOVE.RIB$Behaviour, data.MOVE.RIB$Time_Cat) 
table(data.MOVE.RIB$Time_Cat)
#Feeding
prop.test(c(172,43), c(263,161), alternative = "greater")

# Rib Feeding by Size - Under 25
prop.table(table(data.MOVE.RIB.u25$Behaviour, data.MOVE.RIB.u25$Time_Cat),margin = 2)
table(data.MOVE.RIB.u25$Behaviour, data.MOVE.RIB.u25$Time_Cat) 
table(data.MOVE.RIB.u25$Time_Cat)

# Rib Feeding by Size - Over 25
prop.table(table(data.MOVE.RIB.o25$Behaviour, data.MOVE.RIB.o25$Time_Cat),margin = 2)
table(data.MOVE.RIB.u25$Behaviour, data.MOVE.RIB.u25$Time_Cat) 
table(data.MOVE.RIB.u25$Time_Cat)
#Feeding
prop.test(c(172,43), c(263,161), alternative = "greater")
#Feeding
prop.test(c(172,43), c(263,161), alternative = "greater")


# Gili Feeding
prop.table(table(data.MOVE.MAL$Behaviour, data.MOVE.MAL$Time_Cat),margin = 2)
table(data.MOVE.MAL$Behaviour, data.MOVE.MAL$Time_Cat)
table(data.MOVE.MAL$Time_Cat)
#Feeding
prop.test(c(20,13), c(31,34), alternative = "greater")
# model2.brms1 = brm(Feed~Habitat + Diam_cm + Time_Cat*Reef + Acropora + 
#              COTSPer100 + Coral.Cover, data = data.MOVE, family = bernoulli(),
#            warmup = 200, iter = 1000, chains = 3, cores = 3,
#            control = list(adapt_delta = 0.95))
# model2.brms2 = brm(Feed~Habitat + Diam_cm*Time_Cat + Time_Cat*Reef + Acropora + 
#              COTSPer100 + Coral.Cover, data = data.MOVE, family = bernoulli(),
#            warmup = 200, iter = 1000, chains = 3, cores = 3,
#            control = list(adapt_delta = 0.95))
# summary(model2.brms1)
# summary(model2.brms2)
# plot(marginal_effects(model2.brms2, "Diam_cm:Time_Cat"), points=F)[[1]]
# WAIC(model2.brms1, model2.brms2)
# model2_em <- emmeans(model2.brms2,  ~ Time_Cat | Diam_cm)
# model2_em2 <- emmeans(model2.brms2,  ~  Reef | Time_Cat)
# model2_em3 <- emmeans(model2.brms2,  ~  Habitat)
# 
# 
# #get all possible contrasts
# cont <- contrast(model2_em, "tukey")
# cont2 <- contrast(model2_em2, "tukey")
# cont3 <- contrast(model2_em3, "tukey")
# cont
# cont2
# cont3
# 
# # plot Results
# summy = summary(model2.brms2)$fixed
# p2a = ploteffects.brms(model2.brms2) + theme_bw(base_size = 14) + 
#   # scale_x_discrete(limits = rev(rownames(summy)[-1]),
#   #                  labels = c("Night - Morning : Reef(Rib)", "Afternoon - Morning : Reef(Rib)",
#   #                             "Midday - Morning : Reef(Rib)","% Acropora", "Moving - Resting",
#   #                             "Feeding - Resting",
#   #                             "Rib Reef - Gili Lankanfushi", "Night - Morning", "Afternoon - Morning",
#   #                             "Midday - Morning", "Diameter of Starfish (cm)", "Slope - Crest")) +
#   ylab("Effect Size (Log Odds-Ratio)") + ggtitle("(a)") + xlab("Variable")
# p2a
# # Diameter Marginal Effects plot
# p2b = plot(marginal_effects(model.brms, "Diam_cm:Time_Cat"), points=F)[[1]] 
# p2b = p2b + theme_bw(base_size = 14) + xlab("Diameter of Starfish (cm)") + 
#   ylab("P(Exposure)") + ggtitle("(b)")
# # Time:Reef Marginal effects plot
# p2c = plot(marginal_effects(model.brms, "Time_Cat:Reef"), points=F)[[1]] 
# pc = pc + theme_bw(base_size = 14) + xlab("Time of Observation") + 
#   ylab("P(Exposure)") + scale_color_manual(values = c( "seagreen", "orange")) +
#   theme(legend.position = "left") + ggtitle("(c)")
# # Behaviou Marginal effects plot
# pd = plot(marginal_effects(model.brms, "Behaviour"), points=F)[[1]] 
# pd = pd + theme_bw(base_size = 14) + xlab("Behaviour") + 
#   ylab("P(Exposure)") + ggtitle("(d)")
# 
# # Arrange PLots
# gridExtra::grid.arrange(pa + theme(plot.background = element_rect(size=0.5,linetype="solid",color="grey")),
#                         pb + theme(plot.background = element_rect(size=0.5,linetype="solid",color="grey")),
#                         pc + theme(plot.background = element_rect(size=0.5,linetype="solid",color="grey")),
#                         pd + theme(plot.background = element_rect(size=0.5,linetype="solid",color="grey")),
#                         ncol=2, widths = c(3,2))
# # Save plots
# graph2eps(file="Figure3.eps", aspectr=3, font = "Arial", height = 8, width=10, bg = "transparent")
# graph2tif(file="Figure3.tiff",dpi=300, font = "Arial", height=8, width=10)
