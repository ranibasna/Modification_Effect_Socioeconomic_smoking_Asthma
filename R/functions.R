# imputation functions ----

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


# this function will run a minimal case for a computation reasons for a full case we used the HCP cluster
impute_missing_data_mice <- function(processed_missing_data, miss_seed){
  # check the seed
  if(missing(miss_seed)){
    set.seed(11)
  }else{
    set.seed(miss_seed)
  }
  # prepare the pred matrix
  pred <- make.predictorMatrix(processed_missing_data)
  # due to high depndencies between some variables we imopse which variables to include in each of the imputation models.
  pred[c("hereditery_asthma","herditery_pulldis"),c("herditery_pulldis","hereditery_asthma")] <- 0
  pred[c("hereditery_allergy","herditery_pulldis"),c("herditery_pulldis","hereditery_allergy")] <- 0
  pred[c("hereditery_allergy","hereditery_asthma"),c("hereditery_asthma","hereditery_allergy")] <- 0
  pred[c("hereditery_asthma"),c("c_asthma")] <- 0
  pred[c("her_dis"),c("smoking_status")] <- 0
  # n.core is number of cores and n.imp.cores number of imputations per core
  Imputed_data <- mice(data = processed_missing_data, m = 2, maxit = 2, predictorMatrix = pred, seed = miss_seed)

  return(Imputed_data)
}


# Bayesian Network functions ----

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

plot.network <- function(structure, ht = "400px"){
  nodes.uniq <- unique(c(structure$arcs[,1], structure$arcs[,2]))
  nodes <- data.frame(id = nodes.uniq,
                      label = nodes.uniq,
                      color = "darkturquoise",
                      shadow = TRUE)
  edges <- data.frame(from = structure$arcs[,1],
                      to = structure$arcs[,2],
                      arrows = "to",
                      smooth = TRUE,
                      shadow = TRUE,
                      color = "black")
  return(visNetwork(nodes, edges, height = ht, width = "100%"))
}

drop_vars_df <- function(data_drop, drop_var){
  assertthat::assert_that(is.vector(drop_var))
  data_drop <- data_drop[ , !(names(data_drop) %in% drop_var)]
  return(data_drop)
}

# conditional probability function ----
cpq_effe_modif <- function(.data, vars, outcome, state, model, repeats = 500000) {
  all.levels <- if (any(length(vars) > 1)) {
    lapply(.data[, (names(.data) %in% vars)], levels)
  } else {
    all.levels <- .data %>%
      select(all_of(vars)) %>%
      sapply(levels) %>%
      as_tibble()
  }
  combos <- do.call("expand.grid", c(all.levels, list(stringsAsFactors = FALSE))) # al combiations

  # generate character strings for all combinations
  str1 <- ""
  for (i in seq(nrow(combos))) {
    str1[i] <- paste(combos %>% names(), " = '",
                     combos[i, ] %>% sapply(as.character), "'",
                     sep = "", collapse = ", "
    )
  }

  # repeat the string for more than one outcome
  str1 <- rep(str1, times = length(outcome))
  str1 <- paste("list(", str1, ")", sep = "")

  # repeat loop for outcome variables (can have more than one outcome)
  all.levels.outcome <- if (any(length(outcome) > 1)) {
    lapply(.data[, (names(.data) %in% outcome)], levels)
  } else {
    all.levels <- .data %>%
      select(all_of(outcome)) %>%
      sapply(levels) %>%
      as_tibble()
  }
  combos.outcome <- do.call("expand.grid", c(all.levels.outcome))

  # repeat each outcome for the length of combos
  str3 <- rep(paste("(", outcome, " == '", state, "')", sep = ""), each = length(str1) / length(outcome))

  # fit the model
  fitted <- bn.fit(cextend(model), .data)

  # join all elements of string together
  cmd <- paste("replicate(200,cpquery(fitted, ", str3, ", ", str1, ", method = 'lw', n = ", repeats, "))", sep = "")

  prob <- rep(0, length(str1)) # empty vector for probabilities
  q05 <- rep(0,length(str1))
  q975 <- rep(0,length(str1))
  for (i in seq(length(cmd))) {
    prop_vec = eval(parse(text =  cmd[i]))
    q05[i] = quantile(prop_vec, 0.05)
    q975[i] = quantile(prop_vec,0.975)
    prob[i] <- mean(prop_vec)
  } # for each combination of strings, what is the probability of outcome
  res <- cbind(combos, prob, q05, q975)

  return(res)
}
# get conditional probability plot for one variable
get_cpq_plot_one_var <- function(res_data, effe_modif_vars, original_raw_data, final_data){
  # get the var cases
  var1Cases <- unique(final_data[,effe_modif_vars[1]])
  #var2Cases <- unique(final_data[,effe_modif_vars[2]])
  # get the var labels
  Var1Labels <- names(attributes(unique(original_raw_data[,effe_modif_vars[1]]))$labels)
  #
  var_1_sym <- sym(effe_modif_vars[1])
  # define conditions
  conditions_1 <- purrr::map2(var1Cases, Var1Labels, ~quo( !!var_1_sym == !!.x ~ !!.y))

  # mutate the new coloumns
  res_data <- res_data %>% mutate("{effe_modif_vars[1]}_Cases" := case_when(!!!conditions_1))
  res_colnames <- colnames(res_data)
  var1 <- sym(res_colnames[grepl("_Cases",res_colnames)])
  #prop_p <- ggplot(res_data, aes(x = !!var1, y = prob))  + geom_errorbar( aes(ymin = q05, ymax = q975, color = !!var1), position = position_dodge(0.3), width = 0.2) + geom_point(aes(color = !!var1), position = position_dodge(0.3)) + scale_color_brewer() + theme_classic() + scale_x_discrete(labels = function(x) {stringr::str_wrap(x, width = 16)})
  #prop_p <- ggplot(res_data, aes(x = !!var1, y = prob))  + geom_errorbar( aes(ymin = q05, ymax = q975, color = !!var1), position = position_dodge(0.3), width = 0.2) + geom_point(aes(color = !!var1), position = position_dodge(0.3)) + scale_color_viridis(discrete = TRUE, option = "D") + theme_classic() + scale_x_discrete(labels = function(x) {stringr::str_wrap(x, width = 16)})
  prop_p <- ggplot(res_data, aes(x = !!var1, y = prob))  + geom_errorbar( aes(ymin = q05, ymax = q975, color = !!var1), position = position_dodge(0.3), width = 0.2) + geom_point(aes(color = !!var1), position = position_dodge(0.3)) + theme_classic() + scale_x_discrete(labels = function(x) {stringr::str_wrap(x, width = 16)})
  res_list <- list()
  res_list[[1]] <- res_data
  res_list[[2]] <- prop_p
  return(res_list)
}
# get the result table and the conditional probability plots
get_cpq_plot <- function(res_data, effe_modif_vars, original_raw_data, final_data){
  # get the var cases
  var1Cases <- unique(final_data[,effe_modif_vars[1]])
  var2Cases <- unique(final_data[,effe_modif_vars[2]])
  # get the var labels
  Var1Labels <- names(attributes(unique(original_raw_data[,effe_modif_vars[1]]))$labels)
  Var2Labels <- names(attributes(unique(original_raw_data[,effe_modif_vars[2]]))$labels)
  #
  var_1_sym <- sym(effe_modif_vars[1])
  var_2_sym <- sym(effe_modif_vars[2])
  # define conditions
  conditions_1 <- purrr::map2(var1Cases, Var1Labels, ~quo( !!var_1_sym == !!.x ~ !!.y))
  conditions_2 <- purrr::map2(var2Cases, Var2Labels, ~quo( !!var_2_sym == !!.x ~ !!.y))
  # mutate the new coloumns
  res_data <- res_data %>% mutate("{effe_modif_vars[1]}_Cases" := case_when(!!!conditions_1)) %>% mutate("{effe_modif_vars[2]}cases" := case_when(!!!conditions_2))
  res_colnames <- colnames(res_data)
  var1 <- sym(res_colnames[grepl("_Cases",res_colnames)])
  var2 <- sym(res_colnames[grepl("cases",res_colnames)])
  prop_p <- ggplot(res_data, aes(x = !!var1, y = prob))  + geom_errorbar( aes(ymin = q05, ymax = q975, color = !!var2), position = position_dodge(0.3), width = 0.2) + geom_point(aes(color = !!var2), position = position_dodge(0.3)) + scale_color_manual(values = c("#00AFBB", "#E7B800",'#999999')) + theme_classic() + scale_x_discrete(labels = function(x) {stringr::str_wrap(x, width = 16)})
  res_list <- list()
  res_list[[1]] <- res_data
  res_list[[2]] <- prop_p
  return(res_list)
}
