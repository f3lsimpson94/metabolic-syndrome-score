library(tidyverse)
library(readr)
library(psych)
library(moments)
library("reshape2")
library("moments")
library("readr")
library("janitor")
library("tibble")
#library("MASS")


#-----------Metabolic syndrome score----------------------#


setwd("C:/Users/c3371138/Dropbox/activate-food")
act_food = read.csv("act-food-pccv.csv", stringsAsFactors = TRUE, encoding = "UTF-8")


# ------------- NOTE ----------------------------- #

#1) Censor: set all biomarker values for each individual case that are below the clinical threshold to the clinical threshold
#2) Center: subtract relevant clinical threshold from each biomarker value (optional: winsorize to remove any outliers)

# I have winsorized the data, should be relatively easy to do posthoc though. 

#input the relevent clinical thresholds 

waist_f = 80
waist_m = 94
trig = 1.7
hdl_m = 1.3
hdl_f = 1
sBP = 130
dBP = 85
BlGl = 5.6

#work out what each of them are in the data 
colnames(act_food)

new_cols = c("visit_age", 
             "sex",
             "waist_average" ,
             "TRIGLY",
             "sbp_average.x",
             "dbp_average.x",
             "blood_glucose",
             "HDL_mmoll", 
             "record_id")


met_data <- act_food[, new_cols] %>% 
  as_tibble() %>% 
  na.omit() # NAs bugger everything

#we need to work out whether 1 or 2 is female

met_data %>% na.omit(waist_average) %>% group_by(sex) %>% summarise(mean = mean(waist_average))

#looks like it's male = 1, female = 2

summary(met_data)

#Censor  the data 

met_data = met_data %>% 
  mutate(waist_average = if_else(sex == 1, waist_average - waist_m,waist_average),
         waist_average = if_else(sex == 2, waist_average - waist_f, waist_average),
         TRIGLY = TRIGLY-trig, 
         HDL_mmoll = if_else(sex == 1,HDL_mmoll - hdl_m,HDL_mmoll),
         HDL_mmoll = if_else(sex == 2,HDL_mmoll - hdl_f, HDL_mmoll),
         sbp_average.x= sbp_average.x- sBP,
         dbp_average.x = dbp_average.x - dBP,
         blood_glucose = blood_glucose - BlGl)
summary(met_data)

#set all less than threshold data to zero and inverse hdl

met_data = met_data %>% mutate_at(vars(-HDL_mmoll), funs(if_else(. < 0, 0, .))) %>% 
  mutate(HDL_mmoll = abs(if_else(HDL_mmoll > 0, 0, HDL_mmoll)))

summary(met_data)


#3) Standardize: divide each censored and centered biomarker value by
#the matching standard deviation shown in Table A.1

met_data[,3:8] = scale(met_data[,3:8], center = FALSE, scale = apply(met_data[,3:8], 2, sd, na.rm = TRUE))

#check sd = 1
apply(met_data[,3:8], 2, sd)

# 4) Principal Component Analysis was used to extract six 
# uncorrelated component scores which were standardized to a standard deviation of one;

#standardise the scores
scores[,1:3] = scale(scores[,1:3], center = FALSE, scale = apply(scores[,1:3], 2, sd, na.rm = TRUE))

#check that they have an SD of 0
apply(scores,2,sd)

#6) Square: square each individual value

square = function(x){x^2}

square_scores = apply(scores, c(1,2), square)

head(square_scores)

#7) Sum: take the sum of the squared  values for each individual row

row_sum = rowSums(square_scores)

#8) Square root: take the square root of the sum of the squared values

MetSSS = sqrt(row_sum)

met_data = cbind(met_data, MetSSS)

head(met_data)


# Now run it again for the different groups #######

# ----------- women -----------# 

met_f = met_data %>% filter(sex == 2)

#use a 3 factor

pca_met <- principal(met_f[,3:8], rotate = "varimax", nfactors = 3)
pca_met$loadings

rm(met_f)


#loadings = scale(pca_met$loadings[1:6,]) 
#loadings
loadings = pca_met$loadings[1:6,]

#check the loadings
apply(loadings,2,sd)

#remove the attributes
attr(loadings, "scaled:scale") <- NULL
loadings

#calculate scores using the scaled loadings
scores = scale(met_data[,3:8]) %*% loadings 
head(scores)


#standardise the scores
scores[,1:3] = scale(scores[,1:3], center = FALSE, scale = apply(scores[,1:3], 2, sd, na.rm = TRUE))

#check that they have an SD of 0
apply(scores,2,sd)

#6) Square: square each individual value

square = function(x){x^2}

square_scores = apply(scores, c(1,2), square)

head(square_scores)

#7) Sum: take the sum of the squared  values for each individual row

row_sum = rowSums(square_scores)

#8) Square root: take the square root of the sum of the squared values

MetSSS_f = sqrt(row_sum)

met_data = cbind(met_data, MetSSS_f)

# ------- and now for men -------------#

met_m = met_data %>% filter(sex == 1)

#use a 3 factor

pca_met <- principal(met_m[,3:8], rotate = "varimax", nfactors = 3)
pca_met$loadings

rm(met_m)


#loadings = scale(pca_met$loadings[1:6,]) 
#loadings
loadings = pca_met$loadings[1:6,]

#check the loadings
apply(loadings,2,sd)

#remove the attributes
attr(loadings, "scaled:scale") <- NULL
loadings

#calculate scores using the scaled loadings
scores = scale(met_data[,3:8]) %*% loadings 
head(scores)


#standardise the scores
scores[,1:3] = scale(scores[,1:3], center = FALSE, scale = apply(scores[,1:3], 2, sd, na.rm = TRUE))

#check that they have an SD of 0
apply(scores,2,sd)

#6) Square: square each individual value

square = function(x){x^2}

square_scores = apply(scores, c(1,2), square)

head(square_scores)

#7) Sum: take the sum of the squared  values for each individual row

row_sum = rowSums(square_scores)

#8) Square root: take the square root of the sum of the squared values

MetSSS_m = sqrt(row_sum)

met_data = cbind(met_data, MetSSS_m)
head(met_data)
head(MetSSS_m)



# -------- test the correlation ----------- 
# we want to look at the top row, the only one that's below their cutoff is old

cor = cor(met_data[10:14])[,1] %>% print()

summary(met_data$MetSSS)

# sex differences 
t.test(MetSSS ~ sex, data = met_data)

ggplot(data = met_data, aes(y = MetSSS, x = age_y, color = as.factor(sex), group = sex)) + 
  geom_point() + 
  geom_smooth() + 
  scale_y_continuous(limits = c(0,3), expand = c(0,0)) + 
  theme_classic()

#------------- CODE GRAVEYARD -------------- #

# Sex-specific clinical thresholds
# female score is listed first. 

#Waist circumference, cm ≥80 ≥94
#Triglycerides, mmol/L ≥1.7
#HDL cholesterol, mmol/L ≤1.3 ≤1.0
#Systolic BP, mm Hg ≥130
#Diastolic BP, mm Hg ≥85
#Blood glucose, mmol/L ≥5.6

head(met_data)


#SBP_sd = 13.28
#DBP_sd = 3.14
#Trigly_sd = 0.30
#HDL_sd = 0.17
#Waist_sd = 8.60
#BG_sd = 0.13

#met_data = met_data %>% 
#  mutate(waist_average = waist_average/Waist_sd,
#         trigly = trigly/Trigly_sd,
#         sbp_average.x= sbp_average.x/SBP_sd,
#        dbp_average.x = dbp_average.x/DBP_sd, 
blood_glucose = blood_glucose/BG_sd,
HDL_mmoll = HDL_mmoll/HDL_sd)

summary(met_data)

apply(met_data, 2, sd)