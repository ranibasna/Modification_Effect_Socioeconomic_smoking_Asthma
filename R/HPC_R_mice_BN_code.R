# Libraries and data ----
library(dplyr)
library(haven)
library(mice)
library(bnlearn)
library(parallel)

# load the data 
data_all <- read_dta("~/Soci_smok_proj/OriginalData/Analysis_dataset(60)25.2.21.dta.dta")
data_all <- as.data.frame(data_all)

# Functions ----


Prepare_data <- function(Smok_Soci_data){
  # make it a data frame object
  my_f_data <- as.data.frame(Smok_Soci_data)
  # dropping uneccesary variables
  drops <- c("ID","cohort")
  my_f_data <- my_f_data[ , !(names(my_f_data) %in% drops)]
  # dropping categorized variables 
  my_f_data <- my_f_data %>% select(-c(cbmi, agecateg, cduration))
  # drop variables with high correlations
  drops_corr <- c("weight", "quitage")
  my_f_data <- my_f_data[ , !(names(my_f_data) %in% drops_corr)]
  # removing non-related variables
  drops_rel <- c("birthyear","height","age_group")
  my_f_data <- my_f_data[ , !(names(my_f_data) %in% drops_rel)]
  # drop more correlated vars
  #drops_corr_cat_2 <- c("asthma_treatmnt","varq10b","cbc","varq10c")
  drops_corr_cat_2 <- c("asthma_treatmnt","cbc")
  my_f_data <- my_f_data[ , !(names(my_f_data) %in% drops_corr_cat_2)]
  
  # if mice is running then also s_amount
  my_f_data <- my_f_data %>% select(-c(asthma_diagnosed, asthma, s_amount))
  
  # from mice models out variables
  my_f_data <- my_f_data %>% select(-c(e_smoking))
  
  # corr_3 correlated with the outcome c_asthma
  drops_corr_cat_3 <- c("alle_asthma","noalle_asthma","w_asthma")
  my_f_data <- my_f_data[ , !(names(my_f_data) %in% drops_corr_cat_3)]
  # remove any_smp and only_symptoms
  my_f_data <- my_f_data %>% select(-c(any_smp, only_symptoms)) # did help a lot
  # trt_copd
  my_f_data <- my_f_data %>% select(- trt_copd) # did help
  # edu_credits
  my_f_data <- my_f_data %>% select(- edu_credits) # did help
  
  # converting variables to factor and numerical
  num_cols <- colnames(my_f_data %>% select( c(BMI, age,duration,startage)))
  fact_cols <-  setdiff(colnames(my_f_data), num_cols)
  # mutare to factors and num
  my_f_data <- my_f_data %>% mutate_at(num_cols, as.numeric)
  my_f_data <- my_f_data %>% mutate_at(fact_cols, factor)
  
  return(my_f_data)
}

# to conclude  we need to remove the edu_credit, s_amount, her_dis, herditery_puldis, any_symp, only_symp, trt_copd

# data preprocssing

Preprocess_Data <- function(Ready_Raw_data){
  # check the data typee
  if (!(is.data.frame(Ready_Raw_data))){
    stop("The input  is not a dataframe format")
  }
  # ordinal as orderd variables 
  Ready_Raw_data$jabstatus <- ordered(Ready_Raw_data$jabstatus, levels=1:8)
  Ready_Raw_data$sei_class <- ordered(Ready_Raw_data$sei_class, levels=0:7)
  Ready_Raw_data$smoking_status <- ordered(Ready_Raw_data$smoking_status, levels=0:2)
  
  #data_all$education <- ordered(data_all$education, levels=0:2)
  #Ready_Raw_data$edu_credits <- ordered(Ready_Raw_data$edu_credits, levels=0:6)
  #Ready_Raw_data$s_amount <- ordered(Ready_Raw_data$s_amount,levels=0:3)
  Ready_Raw_data$e_amount <- ordered(Ready_Raw_data$e_amount,levels=1:3)
  
  # more ordinal
  Ready_Raw_data$syk_class <- ordered(Ready_Raw_data$syk_class,levels=0:10)
  Ready_Raw_data$SSY_class <- ordered(Ready_Raw_data$SSY_class,levels=0:10) #(2012 classification)
  
  # remove the 2012 classification
  Ready_Raw_data <- Ready_Raw_data %>% select(- SSY_class)
  
  return(Ready_Raw_data)
}

# in this function we engineer the variables cbc, any_smp, and only_symptoms and remoove the other sym variables
Prepare_data_bn <- function(raw_correlated_data, imputed_data){
  # some preprocessing
  data_all <- raw_correlated_data %>% select(w_asthma, noalle_asthma, alle_asthma)
  data_all <- data_all %>% mutate_all(factor)
  # feature engineer the variables cbc, any_smp, and only_symptoms
  #imputed_data <- imputed_data %>% dplyr::mutate(cbc = case_when(varq8a ==1 & varq8b==1 & varq8c==1 ~ 1))
  imputed_data <- imputed_data %>% dplyr::mutate(cbc = dplyr::if_else(varq8a ==1 & varq8b==1 & varq8c==1,1,0)) %>%
    dplyr::mutate( any_smp = dplyr::if_else(varq7 == 1 | varq8a ==1 | varq8b ==1 | varq8c==1 | varq9 ==1 | varq10a==1 |varq10b==1, 1,0)) %>% 
    dplyr::mutate(only_symptoms = dplyr::if_else(c_asthma ==1 & any_smp == 1, 1,0))  
  imputed_data <- imputed_data %>% dplyr::mutate_at(vars(any_smp, only_symptoms, cbc), factor)   
  # bind the two data set
  data_modeling <- cbind(imputed_data, data_all)
  # remove some variables symptoms
  data_modeling <- data_modeling %>% select(- dplyr::matches('varq'))
  return(data_modeling)
}

# Prepare for imputation with the pred matrix ----

my_data <- Prepare_data(Smok_Soci_data = data_all)
my_data <- Preprocess_Data(Ready_Raw_data = my_data)
# creat pred matrix
pred <- make.predictorMatrix(my_data)
# due to high depndencies between some variables we imopse which variables to include in each of the imp$
pred[c("hereditery_asthma","herditery_pulldis"),c("herditery_pulldis","hereditery_asthma")] <- 0
pred[c("hereditery_allergy","herditery_pulldis"),c("herditery_pulldis","hereditery_allergy")] <- 0
pred[c("hereditery_allergy","hereditery_asthma"),c("hereditery_asthma","hereditery_allergy")] <- 0
pred[c("hereditery_asthma"),c("c_asthma")] <- 0
pred[c("her_dis"),c("smoking_status")] <- 0

set.seed(111)
imp_parl_final <- parlmice(data = my_data, n.core = 8, n.imp.core = 15, maxit = 40, predictorMatrix = pred, seed = 11, cluster.seed = 22)
# saving the mice object as rds
saveRDS(object =  imp_parl_final, file = "~/Soci_smok_proj/Results/final_mice_model.rds")

# Start the Bayesian Network structure learning ----

# select one sample from the multiple imputations
imputed_data_sample <- mice::complete(imp_parl_final, 20)
data_modeling <- Prepare_data_bn(raw_correlated_data = data_all, imputed_data = imputed_data_sample) 

# discretize the tow non-normal variables duration and startage
dis_data <- bnlearn::discretize(data = data_modeling %>% select(duration, startage), method = "hartemink", breaks = 10, ordered = FALSE, ibreaks=60, idisc="quantile")
print(colnames(dis_data))
data_modeling <- data_modeling %>% dplyr::select(-c( duration, startage)) %>% cbind(dis_data) 

print(detectCores())

cl = makeCluster(8)

#check it works.
clusterEvalQ(cl, runif(10))

# Bootstrap HC 50000
str.bootstrap_mice_hc = boot.strength(data_modeling, R = 50000, algorithm = "hc", cluster=cl)

# saving the bnlearn object
saveRDS (object = str.bootstrap_mice_hc, file = "~/Soci_smok_proj/Results/boot_mice_final_hc.rds")

# tabu
# Bootstrap 50000
str.bootstrap_mice_tabu = boot.strength(data_modeling, R = 50000, algorithm = "tabu", cluster=cl)

# saving the bnlearn object
saveRDS(object = str.bootstrap_mice_tabu, file = "~/Soci_smok_proj/Results/boot_mice_final_tabu.rds")

# Cross_validation to learn the structure, 100 10-fold CV. with Hill climbing
cv.bn_str_hc <- bn.cv(data = data_modeling, 'hc', runs = 80, method = "k-fold", loss = "pred-lw-cg", loss.args = list(target ="c_asthma"), cluster = cl)
saveRDS(object = cv.bn_str_hc, file = "~/Soci_smok_proj/Results/cv.bn_str_hc.rds")

# Cross_validation to leaarn the structure, 100 10-fold CV. with tabu
cv.bn_str_tabu <- bn.cv(data = data_modeling, 'tabu', runs = 80, method = "k-fold", loss = "pred-lw-cg", loss.args = list(target ="c_asthma"), cluster = cl)
saveRDS(object = cv.bn_str_tabu, file = "~/Soci_smok_proj/Results/cv.bn_str_tabu.rds")

# stop the cluster.
stopCluster(cl)



