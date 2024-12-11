rm(list=ls())
library(splines)
library(ggplot2)
library(mgcv)
library(gratia)
library(stats)
library(dplyr)

setwd("C:/University of Michigan Dropbox/Deng Feiyang/Feiyang Deng的文件/私人/Biostat Study/BIOSTAT 885/Proj1")


######## Read data ########
heart_data <- read.table("Heart_data.txt", header = TRUE, sep = "", 
                         na.strings = c("NA", "<NA>"), skip = 18)

heart_data$survive <- ifelse(heart_data$fustat == 1, 0, 1)  # Create binary response variable
head(heart_data)



######## EDA ########
# Explore distribution of age
summary(heart_data)

# summary table of survival and age by transplant
heart_data %>%
  summarise(mean_age = mean(age, na.rm = TRUE),
            sd_age = sd(age, na.rm = TRUE),
            min_age = min(age, na.rm = TRUE),
            max_age = max(age, na.rm = TRUE),
            mean_survive = mean(survive, na.rm = TRUE),
            sd_survive = sd(survive, na.rm = TRUE),
            min_survive = min(survive, na.rm = TRUE),
            max_survive = max(survive, na.rm = TRUE))

heart_data %>%
  group_by(transplant) %>%
  summarise(mean_age = mean(age, na.rm = TRUE),
            sd_age = sd(age, na.rm = TRUE),
            min_age = min(age, na.rm = TRUE),
            max_age = max(age, na.rm = TRUE),
            mean_survive = mean(survive, na.rm = TRUE),
            sd_survive = sd(survive, na.rm = TRUE),
            min_survive = min(survive, na.rm = TRUE),
            max_survive = max(survive, na.rm = TRUE))



ggplot(heart_data, aes(x = age)) +
  geom_histogram(binwidth = 10, alpha = 0.7) +
  labs(title = "Age Distribution", x = "Age (Years)", y = "Count")




# Transplant indicator vs survival
ggplot(heart_data, aes(x = factor(transplant), fill = factor(fustat))) +
  geom_bar(position = "fill") +
  labs(title = "Transplant vs Survival", x = "Transplant (1 = Yes, 0 = No)", y = "Proportion") +
  scale_fill_manual(values = c(4, 2), labels = c("Alive", "Dead"))

# Age vs survival
# Create age groups (e.g., bins of 10 years)
heart_data <- heart_data %>%
  mutate(age_group = cut(age, breaks = seq(0, 90, by = 10), right = FALSE))

# Compute survival probability (1 - death rate) for each age group and transplant status
survival_prob <- heart_data %>%
  group_by(age_group, transplant) %>%
  summarise(survival_probability = mean(survive), .groups = "drop")

# Plot survival probabilities
ggplot(survival_prob, aes(x = age_group, y = survival_probability, group = transplant, color = as.factor(transplant), linetype = as.factor(transplant))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    # title = "Survival Probability by Age Group and Transplant Status",
    x = "Age Group",
    y = "Survival Probability",
    color = "Transplant",
    linetype = "Transplant"
  ) +
  scale_color_manual(values = c("red", "blue"), labels = c("No Transplant", "Transplant")) +
  scale_linetype_manual(values = c("solid", "dashed"), labels = c("No Transplant", "Transplant")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


######## GAM ########
# Convert variables to factors if needed
heart_data_f <- heart_data
heart_data_f$transplant <- as.factor(heart_data_f$transplant)  # Convert to factor
heart_data_f$fustat <- as.factor(heart_data_f$fustat)  # Ensure binary response is a factor
heart_data_f$survive <- as.factor(heart_data_f$survive)  # Ensure binary response is a factor

# Fit a GAM with a smooth term for age and interaction with transplant
gam_model <- gam(survive ~ s(age, by = transplant), 
                 data = heart_data_f, 
                 family = binomial(link = "logit"))

# View summary of the GAM
summary(gam_model)
draw(gam_model)
appraise(gam_model, method = "simulate")
k.check(gam_model)


# Predicted probabilities
age_seq <- seq(min(heart_data_f$age), max(heart_data_f$age), length = 100)

# plot_data transplant=0
data_trans_0 <- data.frame(age = age_seq, transplant = 0)
data_trans_0$transplant <- as.factor(data_trans_0$transplant)
pred_0_logodds <- predict(gam_model, newdata = data_trans_0, type = "link", se.fit = TRUE)
pred_0_prob <- plogis(pred_0_logodds$fit)
pred_0_prob_ci_upper <- plogis(pred_0_logodds$fit + 1.96 * pred_0_logodds$se.fit)
pred_0_prob_ci_lower <- plogis(pred_0_logodds$fit - 1.96 * pred_0_logodds$se.fit)
data_p1 <- data.frame(age = age_seq,
                      transplant = 0,
                      pred_prob = pred_0_prob, 
                      pred_prob_ci_upper = pred_0_prob_ci_upper, 
                      pred_prob_ci_lower = pred_0_prob_ci_lower)

# plot_data transplant=1
data_trans_1 <- data.frame(age = age_seq, transplant = 1)
data_trans_1$transplant <- as.factor(data_trans_1$transplant)
pred_1_logodds <- predict(gam_model, newdata = data_trans_1, type = "link", se.fit = TRUE)
pred_1_prob <- plogis(pred_1_logodds$fit)
pred_1_prob_ci_upper <- plogis(pred_1_logodds$fit + 1.96 * pred_1_logodds$se.fit)
pred_1_prob_ci_lower <- plogis(pred_1_logodds$fit - 1.96 * pred_1_logodds$se.fit)
data_p2 <- data.frame(age = age_seq, 
                      transplant = 1,
                      pred_prob = pred_1_prob, 
                      pred_prob_ci_upper = pred_1_prob_ci_upper, 
                      pred_prob_ci_lower = pred_1_prob_ci_lower)

data_plot <- rbind(data_p1, data_p2)

# plot
ggplot(data_plot, aes(x = age, y = pred_prob, color = as.factor(transplant))) +
  geom_line(alpha = 0.7) +
  geom_ribbon(aes(ymin = pred_prob_ci_lower, ymax = pred_prob_ci_upper, fill = as.factor(transplant)), alpha = 0.1) +
  labs(title = "Predicted Probability of Survive", x = "Age (Years)", y = "Probability") +
  scale_color_manual(values = c("red", "blue"), labels = c("Transplant = 0", "Transplant = 1")) +
  scale_fill_manual(values = c("red", "blue"), labels = c("Transplant = 0", "Transplant = 1"))

summary(heart_data_f)
summary(heart_data_f[heart_data_f$transplant == 1, ])


sum(heart_data_f[heart_data_f$transplant == 0, ]$survive==1)










