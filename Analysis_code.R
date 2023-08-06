

# 1 - INITIALIZE ----------------------------------------------------------

# Set working directory
setwd("/Users/iremduman/Desktop/KU Leuven/Thesis /data_23.06.2023")

### Upload necessary packages
library(readbulk)
library(dplyr)
library(reshape2)
library(lme4)
library(ggplot2)
library(robustbase)
library(ggResidpanel)
library(lmerTest)
library(car)
library(estimability)
install.packages("emmeans_1.8.7.tar",
                 repos = NULL, type = "source")
library(emmeans)
library(arm)




### Read in all CSV files in the current working directory
file_list <- list.files(pattern = "\\.csv$")
# Create an empty list to store the modified data frames
csv_list <- list()
# Loop through the CSV files
for (i in seq_along(file_list)) {
  # Read in the CSV file
  csv <- read.csv(file_list[i], header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  # Check the number of columns in the CSV file
  if (ncol(csv) == 41) {
    # Add a missing column with default values
    csv$missing_column <- NA
  }
  # Add the modified CSV file to the list
  csv_list[[i]] <- csv
}
# Combine all the CSV files into a single data frame
df <- do.call(rbind, csv_list)

# Change NAs in 'coherence' column
table(is.na(df$coherence))
df <- df[!is.na(df$coherence),]

# Change columns 'coherence' and 'condition' into factors
# coherence = [0.15, 0.4, 1]
# condition = [1, 2] --> 1: stable, 2: changing
df$coherence <- as.factor(df$coherence)
df$condition <- as.factor(df$condition)

# Change the column name 'loop' to 'block'
colnames(df)[colnames(df) == 'loop'] <- 'block'

# Change 'direction1' and 'direction2' columns to numeric
df$direction1 <- as.numeric(df$direction1)
df$direction2 <- as.numeric(df$direction2)

# Create a variable for the final direction 
# i.e., correct response
df$direction <- ifelse(df$condition == 1, df$direction1, df$direction2)

# Create a variable for accurate responses for indicating the direction
# i.e., direction accuracy
df$acc <- ifelse(df$direction == df$direction_response, 1, 0)

# Create a variable for accurate responses for indicating whether the
# direction of the dots changed or not during the trial
# i.e., change accuracy
df$acc_change <- ifelse((!is.na(df$direction2) & grepl("y", df$change_response)) | 
                          (is.na(df$direction2) & grepl("n", df$change_response)), 1, 0)

### Create a final data frame, ready to use!
# Delete practice trials by deleting rows where 'p_block' is not NA
df_cleaned <- df[is.na(df$p_block), ]
# Remove coherence = 1
df_cleaned <- df_cleaned[df_cleaned$coherence != 1, ]
# Get the index of the 'sub' column
sub_col_index <- which(colnames(df_cleaned) == 'sub')
# Delete the columns before 'sub'
df_cleaned <- df_cleaned[, (sub_col_index):ncol(df_cleaned)]
# Remove the participants that start with 1
df_cleaned <- df_cleaned[!grepl("^1", df_cleaned$sub), ]
# Change the wrongly named participant
df_cleaned$sub[df_cleaned$sub == '54064'] <- '515'
# Get all the unique participants
unique_subs <- unique(df_cleaned$sub)
# Number of participants
N <- length(unique_subs)



# 2 - SNEAK PEAK ----------------------------------------------------------

# Visualize RTs to direction and change questions
hist(df$direction_RT, breaks = 50)
hist(df$change_RT, breaks = 50)

# Calculate and plot mean accuracy of direction question, for each participant
meanerr <- with(df_cleaned,aggregate(acc,by=list(sub=sub),mean))
plot(meanerr$x)

### Check whether anyone is at chance level in direction accuracy
at_chance <- matrix(NA, nrow = N, ncol = 2)
for (i in 1:N) {
  tempDat <- subset(df_cleaned, sub == unique_subs[i])
  chisq_result <- chisq.test(table(tempDat$acc))
  at_chance[i, ] <- c(unique_subs[i], chisq_result$p.value)
}
print(paste('There are', sum(as.numeric(at_chance[,2])>.05)," people at chance level"))
# Find out which participant to exclude
exclude <- at_chance[as.numeric(at_chance[,2]) > 0.05, 1]
# Remove rows where 'sub' column matches values in 'at_chance' 
df_cleaned <- df_cleaned[!(df_cleaned$sub %in% exclude), ]
# Check unique participants
length(unique(df_cleaned$sub))

### Create plots per participant
par(mfrow=c(2,2))
for(i in 1:N){
  tempDat <- subset(df_cleaned, sub == unique_subs[i])
  #main exp
  acc_block <- aggregate(tempDat$acc,by=list(block=tempDat$block),mean)
  change_acc_block <- with(tempDat,aggregate(acc_change,by=list(block=block),mean))
  plot(acc_block,ylab="acc (.) acc_change (x)",frame=F,ylim=c(0,1),pch=19);abline(h=.5,lty=2,col="grey")
  points(change_acc_block,pch=4,col="grey")
  plot(tempDat$direction_RT,frame=F,main=paste('subject',i),ylab="RT",ylim=c(0,5))
  plot(tempDat$change_confidence,frame=F,ylab="change_conf")
}

### Demographics
# Subset the data for unique participant values:
df_unique <- df_cleaned %>% 
  distinct(sub, .keep_all = T)
# Handedness
as.data.frame(table(df_unique$handedness))
# Age
as.data.frame(table(df_unique$age))
mean(df_unique$age)
sd(df_unique$age)
# Gender
as.data.frame(table(df_unique$sex))



# 3 - ANALYSES ------------------------------------------------------------

##### 3.1 - Direction Accuracy --------------------------------------------

### Mixed models for direction accuracy -----------------------------------
# Step 1- Selection of the fixed model: start with direction accuracy by condition, coherence level, and participant
m1_da <- glmer(acc ~ condition * coherence + (1 | sub), data = df_cleaned, family = binomial)
m2_da <- glmer(acc ~ condition * coherence + acc_change + (1 | sub), data = df_cleaned, family = binomial)
m3_da <- glmer(acc ~ condition * coherence * acc_change + (1 | sub), data = df_cleaned, family = binomial)
m4_da <- glmer(acc ~ condition * coherence + change_confidence + (1 | sub), data = df_cleaned, family = binomial)
m5_da <- glmer(acc ~ condition * coherence * change_confidence + (1 | sub), data = df_cleaned, family = binomial) # convergence issue
# Compare
anova(m1_da, m2_da) # Significant; BIC 31536
anova(m1_da, m3_da) # Significant; BIC 31254 
anova(m1_da, m4_da) # Significant; BIC 33147 --> m3_da to the win!
# Step 2- Selection of the random effects:
m6_da <- glmer(acc ~ condition * coherence * acc_change + (condition | sub), data = df_cleaned, family = binomial)
m7_da <- glmer(acc ~ condition * coherence * acc_change + (coherence | sub), data = df_cleaned, family = binomial)
m8_da <- glmer(acc ~ condition * coherence * acc_change + (acc_change | sub), data = df_cleaned, family = binomial)
# Compare
anova(m3_da, m6_da) # Significant; BIC 30879
anova(m3_da, m7_da) # Significant; BIC 31064
anova(m3_da, m8_da) # Significant; BIC 30474 --> model: m8_da
# Step 2 continues...
m9_da <- glmer(acc ~ condition * coherence * acc_change + (acc_change + condition | sub), data = df_cleaned, family = binomial)
m10_da <- glmer(acc ~ condition * coherence * acc_change + (acc_change + coherence | sub), data = df_cleaned, family = binomial)
# Compare
anova(m8_da, m9_da) # Significant; BIC 30275
anova(m8_da, m10_da) # Significant; BIC 30361 --> model: m9_da
# Step 2 still continues...
m11_da <- glmer(acc ~ condition * coherence * acc_change + (acc_change + condition + coherence | sub), 
                data = df_cleaned, family = binomial) # convergence issue

### Model assumptions
# (1) Linearity
# Categorical/dummy coded predictors -> assumption met by definition
# (2) Homogeneity of variance
resid_panel(m9_da)
# (3) No multicolinearity (cutoff=5)
vif(m9_da)

### Independence of errors
### Gray lines indicate plus and minus 2 standard-error bounds (around 95% of residuals)
binnedplot(fitted(m9_da), 
           residuals(m9_da, type = "response"), 
           nclass = NULL, 
           xlab = "Expected Values", 
           ylab = "Average residual", 
           main = "Binned residual plot", 
           cex.pts = 0.8, 
           col.pts = 1, 
           col.int = "gray")

### Model interpretation
Anova(m9_da)
summary(m9_da)
# Odds ratio
oddrat_da <- coef(summary(m9_da))
oddrat_da <- cbind(oddrat_da, exp(oddrat_da[,1]))
colnames(oddrat_da)[ncol(oddrat_da)] <- "odds_ratio"
lower_b <- exp(data.frame(oddrat_da)$Estimate + (-1.96)*data.frame(oddrat_da)$Std..Error)
upper_b <- exp(data.frame(oddrat_da)$Estimate - (-1.96)*data.frame(oddrat_da)$Std..Error)
oddrat_da <- cbind(oddrat_da, lower_b, upper_b)

### Contrasts
# Coherence
cont_coh_m9_da <- data.frame(pairs(emmeans(m9_da, "coherence"), reverse = T))
lower_b <- exp(data.frame(cont_coh_m9_da)$estimate + (-1.96)*data.frame(cont_coh_m9_da)$SE)
upper_b <- exp(data.frame(cont_coh_m9_da)$estimate - (-1.96)*data.frame(cont_coh_m9_da)$SE)
cont_coh_m9_da <- cbind(cont_coh_m9_da, exp(cont_coh_m9_da$estimate))
colnames(cont_coh_m9_da)[ncol(cont_coh_m9_da)] <- "odds_ratio"
cont_coh_m9_da <- cbind(cont_coh_m9_da, lower_b, upper_b)
# Condition
cont_con_m9_da <- data.frame(pairs(emmeans(m9_da, "condition"), reverse = T))
lower_b <- exp(data.frame(cont_con_m9_da)$estimate + (-1.96)*data.frame(cont_con_m9_da)$SE)
upper_b <- exp(data.frame(cont_con_m9_da)$estimate - (-1.96)*data.frame(cont_con_m9_da)$SE)
cont_con_m9_da <- cbind(cont_con_m9_da, exp(cont_con_m9_da$estimate))
colnames(cont_con_m9_da)[ncol(cont_con_m9_da)] <- "odds_ratio"
cont_con_m9_da <- cbind(cont_con_m9_da, lower_b, upper_b)


### SDT --> bias and sensitivity measure for direction accuracy
# Sensitivity --> dprime
# Bias --> c value
# Iterate over unique participant numbers
dprime_da <- matrix(NA, ncol = 7, nrow = 4 * length(unique(df_cleaned$sub)))
for (participant in unique(df_cleaned$sub)) {
  for (coherence in unique(df_cleaned$coherence)) {
    for (condition in unique(df_cleaned$condition)) {
      # Subset the data for the current participant
      df_current <- df_cleaned[df_cleaned$sub == participant & df_cleaned$coherence == coherence & df_cleaned$condition == condition, ]
      # Compute Hit Rate (H) and False Alarm Rate (FA)
      H <- length(df_current$direction_response[df_current$direction == 0 & df_current$direction_response == 0]) / sum(df_current$direction == 0)
      FA <- length(df_current$direction_response[df_current$direction == 180 & df_current$direction_response == 0]) / sum(df_current$direction == 0)
      # Trick
      H[H==1] = .99
      FA[FA==0] = .01
      # Compute the desired values
      p_c <- mean(df_current$acc)
      print(p_c)
      dprime <- qnorm(H) - qnorm(FA)
      c_value <- -.5 * (qnorm(H) + qnorm(FA))
      # Name the condition (or aggregate)
      agg <- ''
      if (coherence == 0.15 & condition == 1) {agg <- 'low_coh_s'}
      else if (coherence == 0.15 & condition == 2) {agg <- 'low_coh_c'}
      else if (coherence == 0.4 & condition == 1) {agg <- 'high_coh_s'}
      else if (coherence == 0.4 & condition == 2) {agg <- 'high_coh_c'}
      coh <- ''
      if (coherence == 0.15) {coh <- 'low'}
      else if (coherence == 0.4) {coh <- 'high'}
      cond <- ''
      if (condition == 1) {cond <- 'stable'}
      else if (condition == 2) {cond <- 'changing'}
      # Print the results for the current participant
      print(paste('Participant:', participant))
      print(paste('p(c) =', p_c, '; dprime =', dprime, '; c =', c_value))
      # Store the results in the matrix
      result_row <- c(participant, agg, cond, coh, p_c, dprime, c_value)
      dprime_da <- rbind(dprime_da, result_row)
    }}}
# Remove empty rows from dprime_da
dprime_da <- dprime_da[complete.cases(dprime_da), ]
# Change the column names
colnames(dprime_da) <- c('sub', 'agg', 'condition', 'coherence',
                             'p_c', 'dprime', 'c_value')
# Change into dataframe
dprime_da <- as.data.frame(dprime_da)
dprime_da$dprime <- as.numeric(dprime_da$dprime)
dprime_da$c_value <- as.numeric(dprime_da$c_value)

### Ttest for dprime
# low coherence changing
t.test(dprime_da$dprime[dprime_da$agg == "low_coh_c"])
# high coherence changing
t.test(dprime_da$dprime[dprime_da$agg == "high_coh_c"])
# low coherence stable
t.test(dprime_da$dprime[dprime_da$agg == "low_coh_s"])
# high coherence stable
t.test(dprime_da$dprime[dprime_da$agg == "high_coh_s"])
# changing vs. stable
t.test(dprime_da$dprime[dprime_da$condition == "changing"], dprime_da$dprime[dprime_da$condition == "stable"], paired = T)
# high coherence vs. low coherence
t.test(dprime_da$dprime[dprime_da$coherence == "high"], dprime_da$dprime[dprime_da$coherence == "low"], paired = T)
# high coherence changing vs. high coherence stable
t.test(dprime_da$dprime[dprime_da$agg == "high_coh_c"], dprime_da$dprime[dprime_da$agg == "high_coh_s"], paired = T)
# low coherence changing vs. low coherence stable
t.test(dprime_da$dprime[dprime_da$agg == "low_coh_c"], dprime_da$dprime[dprime_da$agg == "low_coh_s"], paired = T)
# high coherence changing vs. low coherence changing
t.test(dprime_da$dprime[dprime_da$agg == "high_coh_c"], dprime_da$dprime[dprime_da$agg == "low_coh_c"], paired = T)
# high coherence stable vs. low coherence stable
t.test(dprime_da$dprime[dprime_da$agg == "high_coh_s"], dprime_da$dprime[dprime_da$agg == "low_coh_s"], paired = T)
# high coherence stable vs. low coherence changing
t.test(dprime_da$dprime[dprime_da$agg == "high_coh_s"], dprime_da$dprime[dprime_da$agg == "low_coh_c"], paired = T)
# high coherence changing vs. low coherence stable
t.test(dprime_da$dprime[dprime_da$agg == "high_coh_c"], dprime_da$dprime[dprime_da$agg == "low_coh_s"], paired = T)



# Ttest for c_value
# low coherence changing
t.test(dprime_da$c_value[dprime_da$agg == "low_coh_c"])
# high coherence changing
t.test(dprime_da$c_value[dprime_da$agg == "high_coh_c"])
# low coherence stable
t.test(dprime_da$c_value[dprime_da$agg == "low_coh_s"])
# high coherence stable
t.test(dprime_da$c_value[dprime_da$agg == "high_coh_s"])
# changing vs. stable
t.test(dprime_da$c_value[dprime_da$condition == "changing"], dprime_da$c_value[dprime_da$condition == "stable"], paired = T)
# high coherence vs. low coherence
t.test(dprime_da$c_value[dprime_da$coherence == "high"], dprime_da$c_value[dprime_da$coherence == "low"], paired = T)
# high coherence changing vs. high coherence stable
t.test(dprime_da$c_value[dprime_da$agg == "high_coh_c"], dprime_da$c_value[dprime_da$agg == "high_coh_s"], paired = T)
# low coherence changing vs. low coherence stable
t.test(dprime_da$c_value[dprime_da$agg == "low_coh_c"], dprime_da$c_value[dprime_da$agg == "low_coh_s"], paired = T)
# high coherence changing vs. low coherence changing
t.test(dprime_da$c_value[dprime_da$agg == "high_coh_c"], dprime_da$c_value[dprime_da$agg == "low_coh_c"], paired = T)
# high coherence stable vs. low coherence stable
t.test(dprime_da$c_value[dprime_da$agg == "high_coh_s"], dprime_da$c_value[dprime_da$agg == "low_coh_s"], paired = T)
# high coherence stable vs. low coherence changing
t.test(dprime_da$c_value[dprime_da$agg == "high_coh_s"], dprime_da$c_value[dprime_da$agg == "low_coh_c"], paired = T)



### Mixed models for dprime 
# Choose a model for dprime
dprime1_da <- lmer(dprime~coherence*condition+(1|sub),dprime_da)
dprime2_da <- lmer(dprime~coherence*condition+(coherence|sub),dprime_da) # convergence issue
### Model assumptions
# (1) Linearity
resid_panel(dprime1_da)
# (2) Homogeneity of variance
resid_panel(dprime1_da)
# (3) No multicolinearity (cutoff=5)
vif(dprime1_da)
# (4) Normality
resid_panel(dprime1_da)
### Model interpretation
Anova(dprime1_da)
summary(dprime1_da)
# Odds ratio
oddrat_dprime_da <- coef(summary(dprime1_da))
oddrat_dprime_da <- cbind(oddrat_dprime_da, exp(oddrat_dprime_da[,1]))
colnames(oddrat_dprime_da)[ncol(oddrat_dprime_da)] <- "odds_ratio"
lower_b <- exp(data.frame(oddrat_dprime_da)$Estimate + (-1.96)*data.frame(oddrat_dprime_da)$Std..Error)
upper_b <- exp(data.frame(oddrat_dprime_da)$Estimate - (-1.96)*data.frame(oddrat_dprime_da)$Std..Error)
oddrat_dprime_da <- cbind(oddrat_dprime_da, lower_b, upper_b)
### Contrasts
# Coherence
pairs(emmeans(dprime1_da, "coherence"))
# Condition
pairs(emmeans(dprime1_da, "condition"))

# Mixed models for c value
# Choose a model for c value
cvalue1_da <- lmer(c_value ~ coherence*condition+ (1|sub), dprime_da)
cvalue2_da <- lmer(c_value~coherence*condition+(coherence|sub),dprime_da) # convergence issue
cvalue3_da <- lmer(c_value~coherence*condition+(condition|sub),dprime_da)
# Compare
anova(cvalue1_da, cvalue3_da) # log-likelihood ratio ~ -3 --> model: cvalue3_da 
### Model assumptions
# (1) Linearity
resid_panel(cvalue3_da)
# (2) Homogeneity of variance
resid_panel(cvalue3_da)
# (3) No multicolinearity (cutoff=5)
vif(cvalue3_da)
# (4) Normality
resid_panel(cvalue3_da)
### Model interpretation
Anova(cvalue3_da)
summary(cvalue3_da)
### Odds ratio
oddrat_cvalue_da <- coef(summary(cvalue3_da))
oddrat_cvalue_da <- cbind(oddrat_cvalue_da, exp(oddrat_cvalue_da[,1]))
colnames(oddrat_cvalue_da)[ncol(oddrat_cvalue_da)] <- "odds_ratio"
lower_b <- exp(data.frame(oddrat_cvalue_da)$Estimate + (-1.96)*data.frame(oddrat_cvalue_da)$Std..Error)
upper_b <- exp(data.frame(oddrat_cvalue_da)$Estimate - (-1.96)*data.frame(oddrat_cvalue_da)$Std..Error)
oddrat_cvalue_da <- cbind(oddrat_cvalue_da, lower_b, upper_b)
### Contrasts
# Coherence
pairs(emmeans(cvalue3_da, "coherence"))
# Condition
pairs(emmeans(cvalue3_da, "condition"))


##### 3.2 - Change Accuracy --------------------------------------------

### Mixed models for change accuracy -----------------------------------

# Step 1- Selection of the fixed model: start with change accuracy by condition, coherence level, and participant
m1_ca <- glmer(acc_change ~ condition * coherence + (1 | sub), data = df_cleaned, family = binomial)
m2_ca <- glmer(acc_change ~ condition * coherence + acc + (1 | sub), data = df_cleaned, family = binomial)
m3_ca <- glmer(acc_change ~ condition * coherence * acc + (1 | sub), data = df_cleaned, family = binomial)
m4_ca <- glmer(acc_change ~ condition * coherence + change_confidence + (1 | sub), data = df_cleaned, family = binomial)
m5_ca <- glmer(acc_change ~ condition * coherence * change_confidence + (1 | sub), 
               data = df_cleaned, family = binomial) # convergence issue
# Compare
anova(m1_ca, m2_ca) # Significant; BIC 31774
anova(m1_ca, m3_ca) # Significant; BIC 31486 --> model: m3_ca
anova(m1_ca, m4_ca) # Significant; BIC 33353 
# Step 2- Selection of the random effects:
m6_ca <- glmer(acc_change ~ condition * coherence * acc + (condition | sub), data = df_cleaned, family = binomial)
m7_ca <- glmer(acc_change ~ condition * coherence * acc + (coherence | sub), data = df_cleaned, family = binomial)
m8_ca <- glmer(acc_change ~ condition * coherence * acc + (acc | sub), data = df_cleaned, family = binomial)
# Compare
anova(m3_ca, m6_ca) # Significant; BIC 28682
anova(m3_ca, m7_ca) # Significant; BIC 31306
anova(m3_ca, m8_ca) # Significant; BIC 30696 --> model: m6_ca
# Step 2 continues...
m9_ca <- glmer(acc_change ~ condition * coherence * acc + (condition + acc | sub), 
               data = df_cleaned, family = binomial) # convergence issue
m10_ca <- glmer(acc_change ~ condition * coherence * acc + (condition + coherence | sub), 
                data = df_cleaned, family = binomial) # convergence issue
### Model assumptions
# (1) Linearity
# Categorical/dummy coded predictors -> assumption met by definition
# (2) Homogeneity of variance
resid_panel(m6_ca)
# (3) No multicolinearity (cutoff=5)
vif(m6_ca) # above cutoff
# Control other models to select:
vif(m8_ca) # above cutoff
vif(m7_ca) # above cutoff
vif(m3_ca) # above cutoff
vif(m2_ca) # model: m2_ca
# (2) Homogeneity of variance 
resid_panel(m2_ca)

### Model interpretation
Anova(m2_ca)
summary(m2_ca)
# Odds ratio
oddrat_ca <- coef(summary(m2_ca))
oddrat_ca <- cbind(oddrat_ca, exp(oddrat_ca[,1]))
colnames(oddrat_ca)[ncol(oddrat_ca)] <- "odds_ratio"
lower_b <- exp(data.frame(oddrat_ca)$Estimate + (-1.96)*data.frame(oddrat_ca)$Std..Error)
upper_b <- exp(data.frame(oddrat_ca)$Estimate - (-1.96)*data.frame(oddrat_ca)$Std..Error)
oddrat_ca <- cbind(oddrat_ca, lower_b, upper_b)

### Contrasts
# Coherence
cont_coh_m2_ca <- data.frame(pairs(emmeans(m2_ca, "coherence")))
lower_b <- exp(data.frame(cont_coh_m2_ca)$estimate + (-1.96)*data.frame(cont_coh_m2_ca)$SE)
upper_b <- exp(data.frame(cont_coh_m2_ca)$estimate - (-1.96)*data.frame(cont_coh_m2_ca)$SE)
cont_coh_m2_ca <- cbind(cont_coh_m2_ca, exp(cont_coh_m2_ca$estimate))
colnames(cont_coh_m2_ca)[ncol(cont_coh_m2_ca)] <- "odds_ratio"
cont_coh_m2_ca <- cbind(cont_coh_m2_ca, lower_b, upper_b)

# Condition
cont_cond_m2_ca <- data.frame(pairs(emmeans(m2_ca, "condition")))
lower_b <- exp(data.frame(cont_cond_m2_ca)$estimate + (-1.96)*data.frame(cont_cond_m2_ca)$SE)
upper_b <- exp(data.frame(cont_cond_m2_ca)$estimate - (-1.96)*data.frame(cont_cond_m2_ca)$SE)
cont_cond_m2_ca <- cbind(cont_cond_m2_ca, exp(cont_cond_m2_ca$estimate))
colnames(cont_cond_m2_ca)[ncol(cont_cond_m2_ca)] <- "odds_ratio"
cont_cond_m2_ca <- cbind(cont_cond_m2_ca, lower_b, upper_b)



### SDT --> bias and sensitivity measure for change accuracy

# Sensitivity --> dprime
# Bias --> c value
dprime_ca <- matrix(NA, ncol = 7, nrow = 4 * length(unique(df_cleaned$sub)))
for (participant in unique(df_cleaned$sub)) {
  for (coherence in unique(df_cleaned$coherence)) {
    # Subset the data for the current participant
    df_current <- df_cleaned[df_cleaned$sub == participant & df_cleaned$coherence == coherence, ]
    # Compute Hit Rate (H) and False Alarm Rate (FA)
    H <- length(df_current$direction_response[df_current$condition == 2 & df_current$change_response == "['y']"]) / sum(df_current$condition == 2)
    FA <- length(df_current$direction_response[df_current$condition == 1 & df_current$change_response == "['y']"]) / sum(df_current$condition == 1)
    # Trick
    H[H==1] = .99
    FA[FA==0] = .01
    # Compute the desired values
    p_c <- mean(df_current$acc_change)
    dprime <- qnorm(H) - qnorm(FA)
    c_value <- -.5 * (qnorm(H) + qnorm(FA))
    # Name the condition (or aggregate)
    agg <- ''
    if (coherence == 0.15 & condition == 1) {agg <- 'low_coh_s'}
    else if (coherence == 0.15 & condition == 2) {agg <- 'low_coh_c'}
    else if (coherence == 0.4 & condition == 1) {agg <- 'high_coh_s'}
    else if (coherence == 0.4 & condition == 2) {agg <- 'high_coh_c'}
    coh <- ''
    if (coherence == 0.15) {coh <- 'low'}
    else if (coherence == 0.4) {coh <- 'high'}
    cond <- ''
    if (condition == 1) {cond <- 'stable'}
    else if (condition == 2) {cond <- 'changing'}
    # Store the results in the matrix
    result_row <- c(participant, agg, cond, coh, p_c, dprime, c_value)
    dprime_ca <- rbind(dprime_ca, result_row)
    }}
# Remove empty rows from dprime_ca
dprime_ca <- dprime_ca[complete.cases(dprime_ca), ]
# Change the column names
colnames(dprime_ca) <- c('sub', 'agg', 'condition', 'coherence',
                         'p_c', 'dprime', 'c_value')
# Change into dataframe
dprime_ca <- as.data.frame(dprime_ca)
dprime_ca$dprime <- as.numeric(dprime_ca$dprime)
dprime_ca$c_value <- as.numeric(dprime_ca$c_value)
dprime_ca$coherence <- as.factor(dprime_ca$coherence)

### Ttest for dprime
# low coherence
t.test(dprime_ca$dprime[dprime_ca$coherence == "low"])
# high coherence 
t.test(dprime_ca$dprime[dprime_ca$coherence == "high"])
# high coherence vs. low coherence
t.test(dprime_ca$dprime[dprime_ca$coherence == "high"], dprime_ca$dprime[dprime_ca$coherence == "low"], paired = T)

# Ttest for c_value
# low coherence
t.test(dprime_ca$c_value[dprime_ca$coherence == "low"])
# high coherence changing
t.test(dprime_ca$c_value[dprime_ca$coherence == "high"])
# high coherence vs. low coherence
t.test(dprime_ca$c_value[dprime_ca$coherence == "high"], dprime_ca$c_value[dprime_ca$coherence == "low"], paired = T)


### Mixed models for dprime 
# Choose a model for dprime
dprime1_ca <- lmer(dprime~coherence+(1|sub), dprime_ca)
### Model assumptions
# (1) Linearity
resid_panel(dprime1_ca) # not linear --> not used


### Mixed models for c value 
# Choose a model for c value
cvalue1_ca <- lmer(c_value~coherence+(1|sub), dprime_ca)
### Model assumptions
# (1) Linearity
resid_panel(cvalue1_ca) 
# (2) Homogeneity of variance
resid_panel(cvalue1_ca)
# (3) Normality
resid_panel(cvalue1_ca)
### Model interpretation
Anova(cvalue1_ca)
summary(cvalue1_ca)
# Odds ratio
oddrat_cvalue_ca <- coef(summary(cvalue1_ca))
oddrat_cvalue_ca <- cbind(oddrat_cvalue_ca, exp(oddrat_cvalue_ca[,1]))
colnames(oddrat_cvalue_ca)[ncol(oddrat_cvalue_ca)] <- "odds_ratio"
lower_b <- exp(data.frame(oddrat_cvalue_ca)$Estimate + (-1.96)*data.frame(oddrat_cvalue_ca)$Std..Error)
upper_b <- exp(data.frame(oddrat_cvalue_ca)$Estimate - (-1.96)*data.frame(oddrat_cvalue_ca)$Std..Error)
oddrat_cvalue_ca <- cbind(oddrat_cvalue_ca, lower_b, upper_b)
### Contrasts
# Coherence
pairs(emmeans(cvalue1_ca, "coherence"))
cont_coh_d_ca <- data.frame(pairs(emmeans(cvalue1_ca, "coherence")))
lower_b <- exp(data.frame(cont_coh_d_ca)$estimate + (-1.96)*data.frame(cont_coh_d_ca)$SE)
upper_b <- exp(data.frame(cont_coh_d_ca)$estimate - (-1.96)*data.frame(cont_coh_d_ca)$SE)
cont_coh_d_ca <- cbind(cont_coh_d_ca, exp(cont_coh_d_ca$estimate))
colnames(cont_coh_d_ca)[ncol(cont_coh_d_ca)] <- "odds_ratio"
cont_coh_d_ca <- cbind(cont_coh_d_ca, lower_b, upper_b)






##### 3.3 - Change Confidence --------------------------------------------

### Mixed models for change confidence -----------------------------------

### Model 1
# Fixed effects: direction accuracy and change accuracy
m_conf0 <- lmer(change_confidence ~ acc * acc_change + (1 | sub), data = df_cleaned)
# Random effects:
m_conf1 <- lmer(change_confidence ~ acc * acc_change + (acc | sub), data = df_cleaned)
m_conf2 <- lmer(change_confidence ~ acc * acc_change + (acc_change | sub), data = df_cleaned)
m_conf3 <- lmer(change_confidence ~ acc * acc_change + (acc_change + acc | sub), data = df_cleaned) # convergence issue
# Compare
anova(m_conf0, m_conf1) # Significant; BIC -5951 --> selected
anova(m_conf0, m_conf2) # Significant; BIC -5899 

### Model assumptions
# (1) Linearity
resid_panel(m_conf1)
# (2) Homogeneity of variance
resid_panel(m_conf1)
# (3) No multicolinearity (cutoff=5)
vif(m_conf1)
# (4) Normality
resid_panel(m_conf1)

### Model interpretation
Anova(m_conf1)
summary(m_conf1)

### Odds ratio
oddrat_conf <- coef(summary(m_conf1))
oddrat_conf <- cbind(oddrat_conf, exp(oddrat_conf[,1]))
colnames(oddrat_conf)[ncol(oddrat_conf)] <- "odds_ratio"
lower_b <- exp(data.frame(oddrat_conf)$Estimate + (-1.96)*data.frame(oddrat_conf)$Std..Error)
upper_b <- exp(data.frame(oddrat_conf)$Estimate - (-1.96)*data.frame(oddrat_conf)$Std..Error)
oddrat_conf <- cbind(oddrat_conf, lower_b, upper_b)

### Contrasts
# Coherence
pairs(emmeans(m_conf1, "acc"), pbkrtest.limit = 30240)
# Condition
pairs(emmeans(m_conf1, "acc_change"), pbkrtest.limit = 30240)



### Model 2
# Fixed effects: direction accuracy, change accuracy, coherence, and condition
conf1 <- lmer(change_confidence ~ coherence * condition * acc * acc_change + (1 | sub), data = df_cleaned)
conf2 <- lmer(change_confidence ~ coherence * condition * acc * acc_change + (coherence | sub), data = df_cleaned)
conf3 <- lmer(change_confidence ~ coherence * condition * acc * acc_change + (condition | sub), data = df_cleaned)
conf4 <- lmer(change_confidence ~ coherence * condition * acc * acc_change + (acc | sub), data = df_cleaned)
conf5 <- lmer(change_confidence ~ coherence * condition * acc * acc_change + (acc_change | sub), data = df_cleaned)
# Compare
anova(conf1, conf2) # Significant; BIC -8432 --> selected
anova(conf1, conf3) # Significant; BIC -8052
anova(conf1, conf4) # Significant; BIC -7951
anova(conf1, conf5) # Significant; BIC -7915
# Random slopes
conf6 <- lmer(change_confidence ~ coherence * condition * acc * acc_change + (coherence + condition | sub), data = df_cleaned)
conf7 <- lmer(change_confidence ~ coherence * condition * acc * acc_change + (coherence + acc | sub), data = df_cleaned)
conf8 <- lmer(change_confidence ~ coherence * condition * acc * acc_change + (coherence + acc_change | sub), data = df_cleaned)
# Compare
anova(conf2, conf6) # Significant; BIC -8992 --> selected
anova(conf2, conf7) # Significant; BIC -8634
anova(conf2, conf8) # Significant; BIC -8606
# Random slopes
conf9 <- lmer(change_confidence ~ coherence * condition * acc * acc_change + (coherence + condition + acc | sub), data = df_cleaned)
conf10 <- lmer(change_confidence ~ coherence * condition * acc * acc_change + (coherence + condition + acc_change | sub), data = df_cleaned)
# Compare
anova(conf6, conf9) # Significant; BIC -9134
anova(conf6, conf10) # Significant; BIC -9165 --> selected
# Random slopes
conf11 <- lmer(change_confidence ~ coherence * condition * acc * acc_change + (coherence + condition + acc_change + acc | sub), data = df_cleaned) # convergence issue

### Model assumptions
# (1) Linearity
resid_panel(conf10)
# (2) Homogeneity of variance
resid_panel(conf10)
# (3) No multicolinearity (cutoff=5)
vif(conf10) # not met
# (4) Normality
resid_panel(conf10)


### Model 2 == when acc == 1 & acc_change == 1
# Subset the data for direction accuracy == 1 & change accuracy == 1
df_accurate <- df_cleaned[df_cleaned$acc_change == 1 & df_cleaned$acc == 1, ]
# Fixed effects:
m_conf4 <- lmer(change_confidence ~ condition * coherence + (1 | sub), data = df_accurate)
# Random effects:
m_conf5 <- lmer(change_confidence ~ condition * coherence + (condition | sub), data = df_accurate)
m_conf6 <- lmer(change_confidence ~ condition * coherence + (coherence | sub), data = df_accurate)
# Compare
anova(m_conf4, m_conf5) # Significant; -6986 --> selected
anova(m_conf4, m_conf6) # Significant; -6447 
# Random effects:
m_conf7 <- lmer(change_confidence ~ condition * coherence + (condition + coherence | sub), data = df_accurate)
# Compare
anova(m_conf4, m_conf7) # Significant; -7444 --> selected

### Model assumptions
# (1) Linearity
resid_panel(m_conf7)
# (2) Homogeneity of variance
resid_panel(m_conf7)
# (3) No multicolinearity (cutoff=5)
vif(m_conf7)
# (4) Normality
resid_panel(m_conf7)

### Model interpretation
Anova(m_conf7)
summary(m_conf7)

### Odds ratio
oddrat_conf2 <- coef(summary(m_conf7))
oddrat_conf2 <- cbind(oddrat_conf2, exp(oddrat_conf2[,1]))
colnames(oddrat_conf2)[ncol(oddrat_conf2)] <- "odds_ratio"
lower_b <- exp(data.frame(oddrat_conf2)$Estimate + (-1.96)*data.frame(oddrat_conf2)$Std..Error)
upper_b <- exp(data.frame(oddrat_conf2)$Estimate - (-1.96)*data.frame(oddrat_conf2)$Std..Error)
oddrat_conf2 <- cbind(oddrat_conf2, lower_b, upper_b)









