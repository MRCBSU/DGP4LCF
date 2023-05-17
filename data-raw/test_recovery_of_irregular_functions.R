library(Rcpp)

load("~/DGP4LCF/data/sim_fcs_truth.rda")
load("~/DGP4LCF/data/sim_fcs_init.rda")

sourceCpp("~/DGP4LCF/src/gibbs_within_mcem_irregular_time.cpp")
sourceCpp("~/DGP4LCF/src/gibbs_after_mcem_irregular_time.cpp")

source("~/DGP4LCF/R/mcem_parameters_setup_irregular_time.R")
source("~/DGP4LCF/R/subject_specific_objects.R")
source("~/DGP4LCF/R/table_generator.R")
source("~/DGP4LCF/R/mcem_algorithm_irregular_time.R")
source("~/DGP4LCF/R/mcem_cov_plot.R")
source("~/DGP4LCF/R/gibbs_after_mcem_diff_initials.R")
source("~/DGP4LCF/R/gibbs_after_mcem_algorithm_irregular_time.R")
source("~/DGP4LCF/R/gibbs_after_mcem_load_chains.R")
source("~/DGP4LCF/R/gibbs_after_mcem_combine_chains.R")
source("~/DGP4LCF/R/numerics_summary_need_alignment.R")
source("~/DGP4LCF/R/numerics_summary_do_not_need_alignment.R")
source("~/DGP4LCF/R/factor_score_trajectory.R")
source("~/DGP4LCF/R/factor_loading_heatmap.R")

# ################################# irregular_6_8 ################################################################
# set.seed(456)
#
# # setup: ok
# mcem_parameter_setup_irregular_time_result<-
#   mcem_parameter_setup_irregular_time(p = 100, k = 4, n = 17, q = 8,
#                                       obs_time_num = sim_fcs_truth$obs_time_num,
#                                       obs_time_index = sim_fcs_truth$obs_time_index,
#                                       a_person = sim_fcs_truth$a_person,
#                                       col_person_index = sim_fcs_truth$col_person_index,
#                                       y_init = sim_fcs_init$y_init_irregular,
#                                       a_init = sim_fcs_init$a_init_2,
#                                       z_init = sim_fcs_init$z_init_2,
#                                       phi_init = sim_fcs_init$phi_init_irregular,
#                                       a_full = sim_fcs_truth$a_full,
#                                       train_index = (1:8),
#                                       x = sim_fcs_truth$observed_x_train_irregular)
#
# # mcem: ok
# mcem_algorithm_irregular_time_result<-
#   mcem_algorithm_irregular_time(ind_x = 1,
#                                 x = sim_fcs_truth$observed_x_train_irregular,
#                                 mcem_parameter_setup_result = mcem_parameter_setup_irregular_time_result)
#
# par(mfrow = c(2,4))
# for (em_index in 1:mcem_algorithm_irregular_time_result$index_used){
#   mcem_cov_plot(mcem_algorithm_irregular_time_result$sigmay_record[em_index,,], k = 4, q = 8, title = paste0("MCEM Iteration ", em_index))
# }
#
# mcem_cov_plot(sim_fcs_truth$gp_sigmay_truth, k = 4, q = 10, title = "Truth: Correlated Factors")
#
# # original: used chol for sampling
#
# par(mfrow = c(2,4))
# for (em_index in 1:sim_fcs_results_under_init_2_use_all_fitting$mcem_algorithm_result$index_used){
#   mcem_cov_plot(sim_fcs_results_under_init_2_use_all_fitting$mcem_algorithm_result$sigmay_record[em_index,,], k = 4, q = 8, title = paste0("MCEM Iteration ", em_index))
# }
#
# mcem_cov_plot(sim_fcs_truth$gp_sigmay_truth, k = 4, q = 10, title = "Truth: Correlated Factors")
#
# # initial: ok
# gibbs_after_mcem_diff_initials_irregular_time_result<-
#   gibbs_after_mcem_diff_initials(ind_x = 1,
#                                  tot_chain = 5,
#                                  mcem_parameter_setup_result = mcem_parameter_setup_irregular_time_result,
#                                  mcem_algorithm_result = mcem_algorithm_irregular_time_result)
#
# # gibbs-after-mcem:
#
# tot_chain<- 5
#
# for (chain_index in 1:tot_chain){
#   gibbs_after_mcem_algorithm_irregular_time(chain_index = chain_index,
#                                             mc_num = 10000,
#                                             burnin = 3000,
#                                             thin_step = 10 ,
#                                             pathname = "~",
#                                             test_index = (9:10),
#                                             x = sim_fcs_truth$observed_x_train_irregular,
#                                             gibbs_after_mcem_diff_initials_result = gibbs_after_mcem_diff_initials_irregular_time_result,
#                                             mcem_algorithm_result = mcem_algorithm_irregular_time_result,
#                                             mcem_parameter_setup_result =  mcem_parameter_setup_irregular_time_result)
# }
#
#
# constant_list<- list(num_time_test = 2,
#                      mc_num = 10000,
#                      thin_step = 10,
#                      burnin = 3000,
#                      pathname = "/home/jc2312",
#                      p = 100,
#                      k = 4,
#                      n = 17,
#                      q = 8,
#                      ind_x = 1)
#
# for (chain_index in 1:tot_chain){
#
#   gibbs_after_mcem_load_chains_result<- gibbs_after_mcem_load_chains(chain_index = chain_index,
#                                                                      gibbs_after_mcem_algorithm_result = constant_list)
#
#   save(gibbs_after_mcem_load_chains_result,
#        file = paste0("/home/jc2312/chain_", chain_index,"_result.RData"))
# }
#
# gibbs_after_mcem_combine_chains_irregular_time_result<- gibbs_after_mcem_combine_chains(tot_chain = 5,
#                                                                                         gibbs_after_mcem_algorithm_result = constant_list)
#
# library(coda)
#
# # compare with original:
# numerics_summary_need_alignment_irregular_time_result<-
#   numerics_summary_need_alignment(gibbs_after_mcem_combine_chains_result =  gibbs_after_mcem_combine_chains_irregular_time_result)
#
# numerics_summary_need_alignment_irregular_time_result$convergence_summary
#
# factor_score_trajectory(numerics_summary_need_alignment_irregular_time_result$reordered_summary$latent_y,
#                         factor_index = 1,
#                         person_index = person_index,
#                         trajectory_title = paste0("Irregular-6-8: Estimated Trajectory of Factor 1 for Person ", person_index),
#                         cex_main = 0.8)
#
# mae_irregular_6_8 <- mean(abs(numerics_summary_need_alignment_irregular_time_result$reordered_summary$latent_y[1:8,,] - fcs_real_y_rescaled))
#
# # compare with original: ok
# numerics_summary_do_not_need_alignment_irregular_time_result<-
#   numerics_summary_do_not_need_alignment(pred_x_truth =  sim_fcs_truth$observed_x_pred_reformat,
#                                          gibbs_after_mcem_combine_chains_result =  gibbs_after_mcem_combine_chains_irregular_time_result)
#
# # finally, save the new results for use in the vignette
# sim_fcs_results_under_init_2_use_all_fitting_updated<-
#   list(mcem_algorithm_irregular_time_result = mcem_algorithm_irregular_time_result,
#        gibbs_after_mcem_combine_chains_irregular_time_result = gibbs_after_mcem_combine_chains_irregular_time_result,
#        numerics_summary_do_not_need_alignment_irregular_time_result = numerics_summary_do_not_need_alignment_irregular_time_result,
#        numerics_summary_need_alignment_irregular_time_result = numerics_summary_need_alignment_irregular_time_result)
#
# usethis::use_data(sim_fcs_results_under_init_2_use_all_fitting_updated, overwrite = TRUE)

# ############################################ regular_8_model ############################################
# set.seed(456)
#
# q<- 8
# n<- 17
#
# obs_time_num<- rep(q, times = n)
# obs_time_index<- list()
# a_person<- list()
# col_person_index<- list()
# observed_x_train_regular_8<- list()
#
# for (person_index in 1:n){
#   obs_time_index[[person_index]]<- (1:q)
#   a_person[[person_index]]<- sim_fcs_truth$a_full[(1:q)]
#   col_person_index[[person_index]]<- ((person_index-1)*q + 1):(person_index*q)
#   observed_x_train_regular_8[[person_index]]<- sim_fcs_truth$observed_x_train[,,person_index]
# }
#
# # setup: ok
# mcem_parameter_setup_result<-
#   mcem_parameter_setup_irregular_time(p = 100, k = 4, n = 17, q = 8,
#                        obs_time_num = obs_time_num,
#                        obs_time_index = obs_time_index,
#                        a_person = a_person,
#                        col_person_index = col_person_index,
#                        y_init = sim_fcs_init$y_init,
#                        a_init = sim_fcs_init$a_init_2,
#                        z_init = sim_fcs_init$z_init_2,
#                        phi_init = sim_fcs_init$phi_init,
#                        a_full = sim_fcs_truth$a_full,
#                        x = observed_x_train_regular_8,
#                        train_index = (1:8))
#
# # mcem: ok
#
# mcem_algorithm_result<-
#   mcem_algorithm_irregular_time(ind_x = 1,
#                                 x = observed_x_train_regular_8,
#                                 mcem_parameter_setup_result = mcem_parameter_setup_result)
#
# par(mfrow = c(2,4))
# for (em_index in 1:mcem_algorithm_result$index_used){
#   mcem_cov_plot(mcem_algorithm_result$sigmay_record[em_index,,], k = 4, q = 8, title = paste0("MCEM Iteration ", em_index))
# }
#
# mcem_cov_plot(sim_fcs_truth$gp_sigmay_truth, k = 4, q = 10, title = "Truth: Correlated Factors")
#
# # initial:
# gibbs_after_mcem_diff_initials_result<-
#   gibbs_after_mcem_diff_initials(ind_x = 1,
#                                  tot_chain = 5,
#                                  mcem_parameter_setup_result = mcem_parameter_setup_result,
#                                  mcem_algorithm_result = mcem_algorithm_result)
#
# # gibbs-after-mcem:
#
# tot_chain<- 5
#
# for (chain_index in 1:tot_chain){
#
#   gibbs_after_mcem_algorithm_irregular_time(chain_index = chain_index,
#                                             mc_num = 10000,
#                                             burnin = 3000,
#                                             thin_step = 10 ,
#                                             pathname = "~",
#                                             test_index = (9:10),
#                                             x = observed_x_train_regular_8,
#                                             gibbs_after_mcem_diff_initials_result = gibbs_after_mcem_diff_initials_result,
#                                             mcem_algorithm_result = mcem_algorithm_result,
#                                             mcem_parameter_setup_result = mcem_parameter_setup_result)
# }
#
#
# constant_list<- list(num_time_test = 2,
#                      mc_num = 10000,
#                      thin_step = 10,
#                      burnin = 3000,
#                      pathname = "/home/jc2312",
#                      p = 100,
#                      k = 4,
#                      n = 17,
#                      q = 8,
#                      ind_x = 1)
#
# for (chain_index in 1:tot_chain){
#
#   gibbs_after_mcem_load_chains_result<- gibbs_after_mcem_load_chains(chain_index = chain_index,
#                                                                      gibbs_after_mcem_algorithm_result = constant_list)
#
#   save(gibbs_after_mcem_load_chains_result,
#        file = paste0("/home/jc2312/chain_", chain_index,"_result.RData"))
# }
#
# gibbs_after_mcem_combine_chains_result<- gibbs_after_mcem_combine_chains(tot_chain = 5,
#                                                                          gibbs_after_mcem_algorithm_result = constant_list)
#
# library(coda)
#
# # compare with original:
# numerics_summary_need_alignment_result<-
#   numerics_summary_need_alignment(gibbs_after_mcem_combine_chains_result =  gibbs_after_mcem_combine_chains_result)
#
# numerics_summary_need_alignment_result$convergence_summary
#
# factor_score_trajectory(numerics_summary_need_alignment_result$reordered_summary$latent_y,
#                         factor_index = 1,
#                         person_index = person_index,
#                         trajectory_title = paste0("Regular-8: Estimated Trajectory of Factor 1 for Person ", person_index),
#                         cex_main = 0.8)
#
# ### rescale true factor scores to deal with identifiability issue of the covariance matrix
# # for training data
# q<- 8
# k<- 4
# n<- 17
#
# a_train<-  sim_fcs_truth$a_full[(1:q)]
#
# h3n2_data<- list()
#
# list_temp <- vector("list", k)
# for (list_index in 1:k){
#   list_temp[[list_index]]<- a_train
# }
# h3n2_data$input<- list_temp
#
# fcs_sigma_y_init_truth_for_train_data<- GPFDA::mgpCovMat(Data=h3n2_data, hp=sim_fcs_truth$gp_hp_truth)
#
# # rescale truth
# d_matrix<- diag(sqrt(diag(fcs_sigma_y_init_truth_for_train_data)))
# d_matrix_inv<- solve(d_matrix)
#
# fcs_real_y_rescaled<- array(0, dim = c(q,k,n))
#
# for (person_index in 1:n){
#   fcs_real_y_rescaled[,,person_index]<- matrix(d_matrix_inv%*%as.numeric(sim_fcs_truth$real_y[1:q,,person_index]),
#                                                nrow = q,
#                                                ncol = k)
# }
#
# mae_irregular_8 <- mean(abs(numerics_summary_need_alignment_result$reordered_summary$latent_y[1:8,,] - fcs_real_y_rescaled))
#
# # compare with original: ok
# numerics_summary_do_not_need_alignment_result<-
#   numerics_summary_do_not_need_alignment(pred_x_truth =  sim_fcs_truth$observed_x_pred_reformat,
#                                          gibbs_after_mcem_combine_chains_result =  gibbs_after_mcem_combine_chains_result)
#
# # finally, save the new results for use in the vignette
# sim_fcs_results_under_init_2_updated<-
#   list(mcem_algorithm_result = mcem_algorithm_result,
#        gibbs_after_mcem_combine_chains_result = gibbs_after_mcem_combine_chains_result,
#        numerics_summary_do_not_need_alignment_result = numerics_summary_do_not_need_alignment_result,
#        numerics_summary_need_alignment_result = numerics_summary_need_alignment_result)
#
# setwd("~/DGP4LCF")
# usethis::use_data(sim_fcs_results_under_init_2_updated, overwrite = TRUE)
