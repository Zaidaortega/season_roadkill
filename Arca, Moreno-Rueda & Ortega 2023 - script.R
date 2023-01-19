#################################################################################################################################################################
######################### The distribution of vertebrate roadkill varies by season, surrounding environment and animal class ####################################
############################################# Jesús Arca, Gregorio Moreno-Rueda & Zaida Ortega ##################################################################
############################################# Script author: Zaida Ortega (zaidaortega@usal.es) #################################################################
########################################################## Date: 19-Jan-2023 ####################################################################################
#################################################################################################################################################################

## The study is in review process in European Journal of Wildlife Research, you can read the preprint version here: https://www.researchsquare.com/article/rs-2246792/v1

## Load packages
library(ggcorrplot) # To plot correlations
library(tidyverse)
library(readxl)
library(ggplot2)
library(dplyr)
library(knitr)
library(lme4)
library(forcats) # To manipulate factor variables
library(ggmosaic) # To make mosaic plots
library(ggpubr) # To use ggarrange function to join plots
library(rstatix) # For basic statistics as a chi-squared test
library(sjPlot) # To visualize model results
library(unmarked) # To fit hierarchical models N-mixture
library(nnet) # To fit multinomial regression models
library(broom) # To arrange model results
library(MNLpred) # To estimate predictions from multinomial models


## Load data
roadkill <- read_excel("Data roadkill.xlsx")
str(roadkill)

## 1. Data curation (translating from Spanish to English and giving format to the variables) ###################################################################
################################################################################################################################################################
roadkill <- roadkill %>%
  mutate(section = as.factor(tramo),
         month = as.factor(mes),
         tortuosity = tortuosidad,
         water = agua,
         temp = temperatura,
         humidity = humedad,
         precip = precipitacion,
         radiation = radiacion,
         class = as.factor(atropellos),
         fence = as.factor(if_else(vallado == "si", "Fenced", "Unfenced")),
         month = recode_factor(month, MARZO = "Mar", ABRIL = "Apr", MAYO = "May", JUNIO = "Jun", 
                             JULIO = "Jul", AGOSTO = "Aug", SEPTIEMBRE = "Sep", OCTUBRE = "Oct", 
                             NOVIEMBRE = "Nov", DICIEMBRE = "Dec",
                             ENERO = "Jan", FEBRERO = "Feb"),
         class = recode_factor(class, anfibio = "Amphibians", reptil = "Reptiles", ave = "Birds", mamifero = "Mammals")) %>%
  mutate(season = as.factor(if_else(month == "Mar"| month == "Apr" | month == "May", "Spring", 
                                      if_else(month == "Jun" | month == "Jul" | month == "Aug", "Summer",
                                              if_else(month == "Sep" | month == "Oct" | month == "Nov", "Autumn",
                                                      "Winter"))))) %>%
  mutate(season = fct_relevel(season, c("Spring", "Summer", "Autumn", "Winter")),
         fence = fct_relevel(fence, c("Unfenced", "Fenced"))) %>%
  select(id, section, tortuosity, fence, month, season, class, water, temp, tmax, tmin, precip, humidity, radiation) %>%
  drop_na(class)

str(roadkill)
summary(roadkill)

## 2. Analyses #################################################################################################################################################
################################################################################################################################################################

## 2.1. Roadkills by class #####################################################################################################################################
################################################################################################################################################################
roadkill_by_class <- roadkill %>%
  count(month, section, class) 

roadkill_by_class_summary <- roadkill_by_class %>%
  group_by(class) %>%
  summarise(mean = mean(n),
            sd = sd(n))

Figure2A <- ggplot(roadkill_by_class, aes(x = class, y = n, fill = class)) +
  geom_bar(stat = "identity") +
  labs(y = "Total roadkill events", x = element_blank(), tag = "A") +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = (element_text(size = 12)),
        axis.title = (element_text(size = 14)),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.2, 0.85),
        legend.text = element_text(size = 14),
        legend.background = element_blank(),
        strip.background =  element_blank(),
        strip.text = element_text(size = 12))

roadkill_by_class_summary  <- roadkill %>%
  count(class) %>%
  mutate(total = sum(n)) %>%
  add_row(class = "Total", n = 413) %>%
  select(class, n)

kable(roadkill_by_class_summary , digits = 2, format = "simple", col.names = c("Class", "N"),
      caption = "Total counts of roadkilled vertebrates by class")

roadkill_by_class_summary_month <- roadkill_by_class %>%
  group_by(month, class) %>%
  summarise(mean = mean(n),
            sd = sd(n))

Figure2B <- ggplot(roadkill_by_class_summary_month, aes(x = month, y = mean, group = class, color = class)) +
  geom_line(size = 2, alpha = .5) +
  geom_point(size = 3) +
  labs(y = "Average roadkill events per section", tag = "B") +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = (element_text(size = 14)),
        axis.title = (element_text(size = 14)),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.88, 0.86),
        legend.text = element_text(size = 14),
        legend.background = element_blank(),
        strip.background =  element_blank(),
        strip.text = element_text(size = 12))

Figure2 <- ggarrange(Figure2A, Figure2B, common.legend = TRUE, widths = c(0.75, 1), legend="bottom")
ggsave("Figure 2.jpeg", Figure2, device = "jpeg", width = 10, height = 4)


kable(roadkill_by_class_summary_month, digits = 2, format = "simple", col.namonth = c("Month", "Class", "Mean", "SD"),
      caption = "Annual evolution of roadkill by class")

## 2.2. Factors related to the road ############################################################################################################################
################################################################################################################################################################

## We have two factors that differ between roads: tortuosity (real length / Euclidean distance start - end points) and fence (yes/no), we fit a hierarchical model (N-mixture or Royle-Nichols model) 
str(roadkill)
roadkill_mean_water <- roadkill %>%
  group_by(month, section, .drop = FALSE) %>%
  summarise(water = mean(water, na.remove = FALSE)) 

roadkill_count_full <- roadkill %>%
  count(month, section, .drop = FALSE) 

roadkill_count_water <- left_join(x = roadkill_count_full, y = roadkill_mean_water, by = c("month", "section"))

roadkill_variables <- roadkill %>%
  select(section, month, tortuosity, fence, season, temp, precip)

roadkill_hier_model <- merge(x = roadkill_count_water, y = roadkill_variables, by.x = c("month", "section"), all.x = TRUE)
roadkill_hier_model <- unique(roadkill_hier_model)

roadkill_hier_model_scaled <- roadkill_hier_model %>%
  mutate(water = scale(water),
         tortuosity = scale(tortuosity),
         temp = scale(temp),
         precip = scale(precip)) %>%
  select(month, section, n, water, tortuosity, fence, season, temp, precip)

roadkill_hier_model_wide <- roadkill_hier_model_scaled %>%
  pivot_wider(id_cols = "section", names_from = month, values_from = n, values_fill = 0)
str(roadkill_hier_model_wide)

roadkill_hier_model_scaled_covariates <- roadkill_hier_model_scaled %>%
  select(section, tortuosity, fence) %>%
  drop_na(tortuosity, fence) %>%
  distinct()

roadkill_hier_model_wide_covariates <- merge(x = roadkill_hier_model_wide, y = roadkill_hier_model_scaled_covariates, by = "section", all.x = TRUE)
str(roadkill_hier_model_wide_covariates)

data_roadkill <- unmarkedFramePCount(y = roadkill_hier_model_wide[,2:12], siteCovs=data.frame(tortuosity = roadkill_hier_model_wide_covariates$tortuosity,
                                                                                              fence = roadkill_hier_model_wide_covariates$fence))
str(data_roadkill)
head(data_roadkill)

RNmodel1 <- pcount(formula = ~tortuosity+fence ~tortuosity+fence, data = data_roadkill, K = 150)
plot(RNmodel)
print(RNmodel)

## 2.3. Factor related to both the environment and the animal: Distance to the water ###########################################################################
################################################################################################################################################################
str(roadkill)
shapiro.test(roadkill$water)
hist(roadkill$water)
ggplot(roadkill, aes(x = water))+
  geom_histogram()+
  facet_wrap(~ class) 

kruskal.test(water ~ class, data = roadkill)
pairwise.wilcox.test(roadkill$water, roadkill$class,
                     p.adjust.method = "bonferroni")
Figure3B <- ggplot(roadkill, aes(x = class, y = water/1000)) +
  geom_boxplot(aes(fill = class)) +
  labs(y = "Distance to a water body (km)", x = element_blank()) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = (element_text(size = 16)),
        axis.title = (element_text(size = 16)),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        strip.background =  element_blank(),
        strip.text = element_text(size = 16))

## 2.4. Factors related to environment #########################################################################################################################
################################################################################################################################################################
roadkill_season <- roadkill %>%
  count(season, class, section)

roadkill_count_season <- roadkill %>%
  count(season, class) 

roadkill_summary_season <- roadkill %>%
  group_by(season, class) %>%
  summarise(n = n()) %>%
  mutate(percentage = (n/(sum(n))*100))

table_abs<-table(roadkill$season, roadkill$class, useNA = "no") # total frequency by class and season
table_abs
round(prop.table(table_abs), 2) # contingency table with total proportions
round(prop.table(table_abs, 1), 2) # contingency table with proportions per row
round(prop.table(table_abs, 2), 2) # contingency table with proportions per column
addmargins(table_abs, c(1, 2)) # contingency table with absolute sums at the marginals
round(addmargins(prop.table(table_abs), c(1, 2)), 2) # contingency table with sums of proportions at the marginals

chi_test <- chisq_test(table_abs)
chisq_descriptives(chi_test)
pairwise_comparisons <- pairwise_chisq_gof_test(table_abs)

Figure3A <- ggplot(data = roadkill) +
  geom_mosaic(aes(x = product(class, season), fill = class)) +
  theme_mosaic() +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        legend.position = "none",
        axis.title = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = (element_text(size = 16)),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        strip.background =  element_blank(),
        strip.text = element_text(size = 16))

## Model interaction class and season
model_season <- glmer(n ~ class*season + (1|section), data = roadkill_season, family = "poisson")
summary(model_season)
#plot(model_season)
qqnorm(resid(model_season))
tab_model(model_season)


## Effect of the different environmental variables on roadkill of the different vertebrate classes
environmental_variables <- roadkill %>%
  mutate(Temperature = temp,
         Precipitation = precip,
         Humidity = humidity,
         Radiation = radiation,
         Tmax = tmax,
         Tmin = tmin) 

environmental_variables <- environmental_variables %>%
  select(Temperature, Tmax, Tmin, Precipitation, Humidity, Radiation)

correlations <- round(cor(environmental_variables), 3)
head(correlations[, 1:6])

plot_environmental_variables <- ggcorrplot(correlations, hc.order = TRUE, type = "lower", lab = TRUE) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = (element_text(size = 12)),
        legend.title = element_text(size = 11),
        legend.key = element_blank(),
        legend.background = element_blank(),
        strip.background =  element_blank(),
        strip.text = element_text(size = 12))
ggsave("ESM-1.jpeg", plot_environmental_variables, device = "jpeg", width = 8, height = 8)

roadkill_by_class_month <- roadkill %>%
  count(month, class)

roadkill_by_class_month_percentage <- roadkill_by_class_month %>%
  mutate(percentage = n/sum(n)*1000)

roadkill_full <- left_join(roadkill, roadkill_by_class_month_percentage) %>%
  select(id, section, tortuosity, fence, month, season, class, water, temp,
         tmax, tmin, precip, humidity, radiation, percentage)
str(roadkill_full)

roadkill_full_long <- roadkill_full %>%
  pivot_longer(cols = c(temp, humidity, radiation, percentage), names_to = "variable", values_to = "valor") %>%
  mutate(variable = as.factor(variable))
str(roadkill_full_long)

Figure4 <- ggplot(roadkill_full_long, aes(x = month, y = valor, group = variable, color = variable)) +
  geom_line(size = 1.5, alpha = .5) +
  geom_point(size = 2) +
  geom_vline(aes(xintercept = which(levels(month) == "Mar") + 0.75), linetype = "dashed") +
  geom_vline(aes(xintercept = which(levels(month) == "Jun") + 0.75), linetype = "dashed") +
  geom_vline(aes(xintercept = which(levels(month) == "Sep") + 0.75), linetype = "dashed") +
  geom_vline(aes(xintercept = which(levels(month) == "Dec") + 0.75), linetype = "dashed") +
  scale_color_manual(values=c("lightsteelblue", "red", "khaki", "wheat3")) +
  labs(y = "Value") +
  facet_wrap(~ class) +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.4, hjust= 0.3),
        axis.text = (element_text(size = 12)),
        axis.title = (element_text(size = 14)),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background =  element_blank(),
        strip.text = element_text(size = 12))
ggsave("Figure 4.jpeg", Figure4, device = "jpeg", width = 8, height = 6)

roadkill_environmental_models <- left_join(roadkill_full, roadkill_by_class_month)

## Multinomial logistic regression

model_multinomial <- multinom(class ~ temp + precip, 
                              data = roadkill_environmental_models, Hess = TRUE)
summary(model_multinomial)
z_ENG <- summary(model_multinomial)$coefficients/summary(model_multinomial)$standard.errors
z_ENG
p_ENG <- (1 - pnorm(abs(z_ENG), 0, 1)) * 2
p_ENG

exp(coef(model_multinomial))

pred_temp <- mnl_pred_ova(model = model_multinomial,
                      data = roadkill_environmental_models,
                      x = "temp",
                      by = 1,
                      seed = "random", # default
                      nsim = 100, # faster
                      probs = c(0.025, 0.975))
pred_temp$plotdata %>% head()

Figure3C <- ggplot(data = pred_temp$plotdata, aes(x = temp, y = mean, ymin = lower, ymax = upper)) +
  geom_ribbon(aes(fill = class), alpha = 0.3) + # Confidence intervals
  geom_line(aes(color = class), size = 1) + # Mean
  facet_wrap(.~ class, ncol = 4) +
  labs(y = "Roadkill probability", x = "Mean temperature (º C)") +
  theme(panel.border = element_rect(fill = NA, colour = "black"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, vjust = 1, size = 20),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        strip.background =  element_blank(),
        strip.text = element_text(size = 18))

Figure3 <- ggarrange(ggarrange(Figure3A, Figure3B, ncol = 2, labels = c("A", "B"), font.label = list(size = 20)), 
                              Figure3C, nrow = 2, labels = c("","C"), font.label = list(size = 20)) 
ggsave("Figure 3.jpeg", Figure3, device = "jpeg", width = 12, height = 12)