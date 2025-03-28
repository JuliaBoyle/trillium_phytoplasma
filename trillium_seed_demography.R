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

cor.test(seed_sample$lf_length, seed_sample$seed_count) #leaf length & seed count = 0.575

reg3 <- glm(seed_count~lf_length, data=seed_sample, family="poisson")
summary(reg3)

#trying to get predicted fruit size values from model
predicted_length3 <- predict.glm(reg3, type="response")

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

seed_tril_plot<-ggplot(data = all_data, aes(x = DH_type, y = seed_count)) +
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
# Demography data from 2024 with more age classes
#juvenile count= no. plants with 3 leaves but no flower
#adult count = no. plants with 3 leaves and flower but no infection
#infected adult count = no. plants with 3 leaves and flower and infection
#one leaf count = no. plants at 1 leaf stage

all_data<-all_data %>% mutate(Total_trilliums=(infected_adult+adult_count+juvenile_count+one_leaf_count+cot_count)) %>% mutate(Freq_infected_adult=infected_adult/(infected_adult+adult_count)) %>% mutate(count_3leaved=(juvenile_count+infected_adult+adult_count))
all_data<-all_data %>% mutate(Ratio_one_to_repro=(one_leaf_count/(infected_adult+adult_count)))
all_data<-all_data %>% mutate(Ratio_one_to_three=(one_leaf_count/(infected_adult+adult_count+juvenile_count)))
all_data<-all_data %>% mutate(Ratio_veg_to_repro=(juvenile_count/(infected_adult+adult_count)))
#all_data$Freq_infected_adult[is.na(all_data$Freq_infected_adult)] <- 0 
#all_data$Freq_infected_adult[is.nan(all_data$Freq_infected_adult)] <- 0 
all_data<-all_data %>% mutate(across(where(is.numeric), \(x) ifelse(is.nan(x), NA, x)))
all_data$Ratio_one_to_repro[is.infinite(all_data$Ratio_one_to_repro)] <- NA
all_data$Ratio_veg_to_repro[is.infinite(all_data$Ratio_veg_to_repro)] <- NA

cor.test(all_data$Freq_infected_adult, all_data$Ratio_one_to_repro)
library(bestNormalize)
hist(predict(bestNormalize(all_data$Ratio_one_to_repro)))
demography2024_one<-ggplot(data=all_data, aes(Freq_infected_adult, Ratio_one_to_repro))+geom_point()+geom_smooth(method="lm", color="darkgreen")+theme_cowplot()+ylab("1-leaved:flowering ratio")+xlab("Proportion of flowering trillium infected")+labs(subtitle = "r = -0.082, p = 0.336")

cor.test(all_data$Freq_infected_adult, all_data$Ratio_veg_to_repro)
demography_2024_veg<-all_data %>% ggplot(aes(Freq_infected_adult, Ratio_veg_to_repro))+geom_point()+geom_smooth(method="lm", color="darkgreen")+theme_cowplot()+ylab("3-leaved non-flowering:flowering ratio")+xlab("Proportion of flowering trillium infected")+ggtitle("2024")+labs(subtitle = "r = -0.079, p = 0.352")

cor.test(all_data$Freq_infected_adult, all_data$Ratio_one_to_three)
ggplot(data=all_data, aes(Freq_infected_adult, Ratio_one_to_three))+geom_point()+geom_smooth(method="lm", color="darkgreen")+theme_cowplot()+ylab("1-leaved:3-leaved ratio")+xlab("Proportion of flowering trillium infected")+labs(subtitle = "r = -0.023, p = 0.784")
