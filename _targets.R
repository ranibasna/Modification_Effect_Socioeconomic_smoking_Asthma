# _targets.R
library(targets)
library(tarchetypes)
source("R/functions.R")
tar_option_set(packages = c("haven", "dplyr", "mice","tibble","VIM","rmarkdown","naniar","finalfit","ggplot2","bnlearn","visNetwork","cowplot"))
list(
  #record state of R environment
  tar_target(lockenv,renv::snapshot()),
  # set the data path
  tar_target(raw_data_file,"data/Analysis_dataset(60)25.2.21.dta.dta",format = "file"),
  # read the raw dta data
  tar_target(raw_data,read_dta(raw_data_file)),
  # set mice object path
  tar_target(mice_object, "Intermediate/outputs/final_mice_model.rds", format = "file"),
  # load the object
  tar_target(loaded_mice_object, readRDS(mice_object)),
  # select one sample form the loaded object
  tar_target(sample_selection, mice::complete(loaded_mice_object,20)),
  # render the report for the missing data
  tar_render(report, "Intermediate/Reports/miss_data_validation_report.Rmd"),
  # Prepare the data for the BN analysis
  tar_target(prepare_data_bn, Prepare_data_bn(raw_data, sample_selection)),
  # discretize the data
  tar_target(dis_data, bnlearn::discretize(data = prepare_data_bn %>% select(duration, startage), method = "hartemink", breaks = 3,
                                           ordered = TRUE, ibreaks=60, idisc="quantile")),
  # get the ready data for the fitting and Cross-Validation by binding
  tar_target(data_modeling, prepare_data_bn %>% dplyr::select(-c( duration, startage)) %>% cbind(dis_data)),
  # set the learned BN with Hill climbing file path
  tar_target(BN_HC_object, "Intermediate/outputs/boot_mice_final_hc.rds", format = "file"),
  # load the BN_HC_object
  tar_target(str_boot_HC, readRDS(BN_HC_object)),
  # Process the strength object
  tar_target(str_boot_HC_proc, str_boot_HC[with(str_boot_HC, strength >= 0.95 & direction >= 0.5), ]),
  # get the averaged network
  tar_target(avg.simpler_mice_hc, averaged.network(str_boot_HC_proc, threshold = 0.95)),
  # Cross_validation the netwrok
  # tar_target(cv.bn_hc, bnlearn::bn.cv(data = data_modeling, bn = cextend(avg.simpler_mice_hc), 
  #                                    runs = 4, method = "k-fold", loss = "pred-lw-cg", loss.args = list(target ="c_asthma"))),
  #tar_target(learn_bn_structure,  bnlearn::hc(prepare_data_bn, score = "bic-cg")),
  # render the report of r the BN analysis
  tar_render(report_bn, "Intermediate/Reports/BN_Analysis.Rmd")
)

