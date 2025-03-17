library(googlesheets4)
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(ggeffects)
library(dplyr)
library(vegan)
library(sf)
library(tmap)
library(leaflet)
library(cowplot)

#This is the code for 2024 seed count fitness data.
#We collected morphology and seed count data from 100 infected trillium (only 81 could have their seeds counted).
#Paired with each infected trillium, we measured the morphology of the closest asymptomatic trillium. 
#To predict seed counts of asymptomatic trillium based on size, we measured and counted the seeds of 100 trillium chosen to maximize natural variation.
#Here, we predict the number of seeds the symptomatic and asymptomatic plants made based on their size. 

### ALL DATA UPLOADS & MANIPULATION ###
disease_data <- read_sheet("https://docs.google.com/spreadsheets/d/1gXR8N1B6rqh5kZQsVuEUjdG6C9Wib4BAoBgX3hq3ah0/edit?gid=1676613822#gid=1676613822", sheet = "disease")

healthy_data <- read_sheet("https://docs.google.com/spreadsheets/d/1gXR8N1B6rqh5kZQsVuEUjdG6C9Wib4BAoBgX3hq3ah0/edit?gid=1676613822#gid=1676613822", sheet = "healthy")

all_data <- rbind(disease_data, healthy_data) %>%
  filter(!is.na(Date))

summary(all_data$pair_number)
#adding a column to pair the data
all_data$pair_number <- all_data$ID
all_data$pair_number <- substr(all_data$pair_number, start = 2, stop = nchar(all_data$pair_number))
all_data$pair_number <- as.numeric(all_data$pair_number)
# MB: making this a character string instead of numeric values just in case that was causing the problem with the models

#adding a column to make D=1, H=0
all_data$DH_number <- all_data$DH_type
all_data$DH_number[all_data$DH_number == "D"] <- 1
all_data$DH_number[all_data$DH_number == "H"] <- 0
all_data$DH_number = as.numeric(all_data$DH_number)

### WORKING ON TRANSECTS ###
# transects - make tall dataframe w counts 
transects = all_data %>% 
  select(ID, transect_number, R_6m:l_4m) %>% 
  # mutate_all(as.character()) %>% 
  pivot_longer(!c(ID, transect_number)) %>% 
  separate(value,into = c("red", "infected", "adult", "juvenile"), fill = "left", sep = ",") %>% 
  mutate(adult = as.numeric(str_trim(adult))) %>% 
  mutate(juvenile = as.numeric(str_trim(juvenile)))

summary(transects)

transects$red_count = if_else(!is.na(transects$red)&str_detect(transects$red, "red"), as.numeric(str_trim(str_remove(transects$red, "red"))), 0)
transects$red_count = if_else(!is.na(transects$infected)&str_detect(transects$infected, "red"), as.numeric(str_trim(str_remove(transects$infected, "red"))), transects$red_count)
#1st - if inf has red - str.detect, then infected count = 0, then transects_infected_count = str.remove(transects_infected"i") - this removes the i

transects$infected_count = if_else(!is.na(transects$infected)&str_detect(transects$infected, "i"), as.numeric(str_trim(str_remove(transects$infected, "i"))), 0)


### !! SEED COUNT DATA !! ###

seed_count <- read_sheet("https://docs.google.com/spreadsheets/d/1GVjfkGL2Zw95SNeJ0TAr3I1OgY0e0J5ucuaSv9F2QoU/edit?gid=0#gid=0", sheet = "Sheet1")

seed_sample <- read_sheet("https://docs.google.com/spreadsheets/d/1LAnaoIIN2YkGvapiwCqjUMaH5jZ6WNj4UQi-IdbSG8Q/edit?gid=0#gid=0", sheet = "Sheet1") %>% 
  left_join(seed_count)


###seed sample data analysis (cor & reg)

ggplot(data=seed_sample, aes(x=fruit_width, y=seed_count))+
  geom_point()+
  geom_smooth(method="lm")
#width is slightly stronger based on cor.test

ggplot(data=seed_sample, aes(x=fruit_length, y=seed_count))+
  geom_point()+
  geom_smooth(method="lm")

cor.test(seed_sample$stem_diameter_mm, seed_sample$leaf_length_cm)
#length & seed count = 0.8593
#width & seed count = 0.8763
#leaf length & seed count = 0.5890
#stem diam & seed count = 0.6409

reg1 <- glm(seed_count~fruit_width, data = seed_sample, family="poisson")
summary(reg1)
#significant <0.001

reg2 <- lm(seed_count~fruit_length, data=seed_sample)
summary(reg2)
#significant <0.001

reg3 <- glm(seed_count~lf_length, data=seed_sample, family="poisson")
summary(reg3)
#significant <0.001 (not as sig as length & width)

reg4 <- lm(seed_count~stem_diameter_mm, data=seed_sample)
summary(reg4)
#significant

reg5 <- glm(seed_count~stem_diameter_mm + lf_length, data=seed_sample, family="poisson")
summary(reg5)

#trying to get predicted fruit size values from reg2
predicted_length <- predict.glm(reg1, type="response")
predicted_length2 <- predict(reg2)
predicted_length3 <- predict.glm(reg3, type="response")
predicted_length4 <- predict(reg4)
predicted_length5 <- predict.glm(reg5, type="response")

r_squared <- sapply(list(reg1, reg2, reg3, reg4, reg5), function(model) summary(model)$r.squared)
print(r_squared)
#highest r2 is reg1 - so WIDTH is best predictor, reg3/leaf length is worst
# are the r2 values for leaf & stem diam too low 
  #to use to get a measurement for healthy indiv?

seed_sample$predicted_count=predicted_length5

ggplot(data=seed_sample, aes(predicted_count, seed_count))+
  geom_point()+
  geom_abline()

seed_sample$predicted_w=predicted_length

ggplot(data=seed_sample, aes(predicted_w, seed_count))+
  geom_point()+
  geom_abline()


#predicted seed count for H&D from leaf
all_data$predicted_seed_count_l <- predict.glm(reg3, type="response", newdata=all_data)

ggplot(data=all_data, aes(lf_length, predicted_seed_count_l))+
  geom_point()+
  geom_abline()

#replacing w actual seed counts in ONLY the HEALTHY datasheet (manually did disease)
all_data$seed_count[is.na(all_data$seed_count)] <- all_data$predicted_seed_count_l[is.na(all_data$seed_count)]
print(all_data %>% select(ID,predicted_seed_count_l))
write.csv(all_data %>% select(ID,predicted_seed_count_l),"trillium_seed_test.csv")
write.csv(all_data %>% select(ID,pair_number),"trillium_seed_pairs.csv")

all_data$seed_count[all_data$seed_count_NOTES == "not relocated"] <- NA

anyNA(all_data$seed_count)
as.numeric(all_data$seed_count)
summary(all_data$seed_count)
#some are character bc I typed NA into google sheets for the plants that were
  #not relocated, so they aren't 0, but blank instead

# PPP: seed_count ----
p3 <- lmer(data=all_data, seed_count ~ DH_type + (1|pair_number)) #
summary(p3)
plot(ggpredict(p3, terms = c("DH_type")))

pre3 <- ggpredict(p3, terms = c("DH_type"))

p4 <- lmer(data=all_data, predicted_seed_count_l ~ DH_type + (1|pair_number))
summary(p4)
plot(ggpredict(p4, terms = c("DH_type")))
pre4 <- ggpredict(p4, terms = "DH_type")

ggplot(data = all_data, aes(x = DH_type, y = seed_count)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.4, size = 0.8) +
  geom_pointrange(data = filter(pre3, x=="D"), aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low, 
                                   color = "Actual values"), linewidth = 0.7, size = 0.5, 
                  position = position_nudge(x = 0.05)) +
  geom_pointrange(data = pre4, aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low, 
                                   color = "Predicted values"),linewidth = 0.7, size = 0.5, 
                  position = position_nudge(x = -0.05)) +
  scale_color_manual(values = c("Predicted values" = "lightblue", "Actual values" = "darkgreen"),
                     name = "Predictions",
                     labels = c("Actual values", "Predicted values")) +theme_cowplot() +
  xlab("Infection status") +
  ylab("Seed count") +
  scale_x_discrete(labels = c("D" = "Symptomatic", "H" = "Asymptomatic"))+
  #ggtitle("Disease vs Healthy Seed Counts") +
  theme(legend.position="top")+
  guides(color = guide_legend(title.position = "top"))



##################################################################################################

all_data %>% ggplot(aes(x=, y=))
