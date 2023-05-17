#' Generating posterior samples for parameters (other than DGP parameters) in the model and predicted gene expression for one chain.
#'
#' @details This function corresponds to Algorithm 2: Step 1 in the main manuscript; therefore reader can consult the paper for more explanations.
#'
#' @param chain_index A numeric scalar. Index of the chain.
#' @param mc_num A numeric scalar. Number of iterations in the Gibbs sampler.
#' @param burnin A numeric scalar. Number of iterations to be discarded as 'burn-in'.
#' @param thin_step A numeric scalar. This function will only save every 'thin_step'th iteration results in the specified directory to reduce storage space needed. Note that this number can be different from that used in the function 'mcem_algorithm'.
#' @param pathname A character. The directory where the saved Gibbs samplers are stored.
#' @param test_index Index of test time points in the full time vector.
#' @param x  A list of n elements. Each element is a matrix of dimension (p, q_i), storing the gene expression observed at q_i time points for the ith subject.
#' @param mcem_parameter_setup_result A list of objects returned from the function 'mcem_parameter_setup'.
#' @param mcem_algorithm_result A list of objects returned from the function 'mcem_algorithm'.
#' @param gibbs_after_mcem_diff_initials_result A list of objects returned from the function 'gibbs_after_mcem_diff_initials'.
#'
#' @examples
#' # See examples in vignette
#' vignette("bsfadgp_regular_data_example", package = "DGP4LCF")
#' vignette("bsfadgp_irregular_data_example", package = "DGP4LCF")
#'
#' @return Posterior samples for parameters (other than DGP parameters) in the model and predicted gene expression for one chain.
#' @export
gibbs_after_mcem_algorithm<- function(chain_index,
                                      mc_num,
                                      burnin,
                                      thin_step,
                                      pathname,
                                      test_index,
                                      x,
                                      gibbs_after_mcem_diff_initials_result,
                                      mcem_algorithm_result,
                                      mcem_parameter_setup_result){

  # assign results to objects
  p<- mcem_parameter_setup_result$p
  k<- mcem_parameter_setup_result$k
  n<- mcem_parameter_setup_result$n
  q<- mcem_parameter_setup_result$q
  model_dgp<- mcem_parameter_setup_result$model_dgp

  a_full<- mcem_parameter_setup_result$a_full
  a_person = mcem_parameter_setup_result$a_person
  obs_time_index = mcem_parameter_setup_result$obs_time_index
  obs_time_num = mcem_parameter_setup_result$obs_time_num

  ind_x<- gibbs_after_mcem_diff_initials_result$ind_x

  if (ind_x){
    mu_g <- mcem_parameter_setup_result$mu_g
    gene_person_mean_init <- gibbs_after_mcem_diff_initials_result$gene_person_mean_init_multiple_chains
  }

  index_used<- mcem_algorithm_result$index_used
  hyper_record<-  mcem_algorithm_result$hyper_record
  prior_sparsity<- mcem_algorithm_result$prior_sparsity
  ig_parameter<- mcem_algorithm_result$ig_parameter
  obs_gene_vector<- mcem_algorithm_result$obs_gene_vector
  obs_time_vector<- mcem_algorithm_result$obs_time_vector

  y_init <- gibbs_after_mcem_diff_initials_result$y_init_multiple_chains
  a_init <- gibbs_after_mcem_diff_initials_result$a_init_multiple_chains
  z_init <- gibbs_after_mcem_diff_initials_result$z_init_multiple_chains
  phi_init <- gibbs_after_mcem_diff_initials_result$phi_init_multiple_chains
  pai_init <- gibbs_after_mcem_diff_initials_result$pai_init_multiple_chains
  beta_init <- gibbs_after_mcem_diff_initials_result$beta_init_multiple_chains

  # remove
  rm(mcem_parameter_setup_result,
     mcem_algorithm_result,
     gibbs_after_mcem_diff_initials_result)

  ###############################################################################################################################################
  ################################################################# create gibbs parameter objects ##############################################
  ###############################################################################################################################################
  set.seed(chain_index)

  big_a<- a_init[chain_index,,]
  big_z<- z_init[chain_index,,]
  phi<- phi_init[,chain_index]
  pai<- rep(pai_init[[chain_index]],times=k)
  beta<- rep(beta_init[[chain_index]], times = k)

  if (ind_x){

    individual_mean<- gene_person_mean_init[chain_index,,]

    variance_g<- apply(gene_person_mean_init[chain_index,,], 1, var)

  } else {

      individual_mean<- matrix(0, nrow = p, ncol = n)

      variance_g<- rep(0, times = p)

      mu_g<- rep(0, times = p)
  }


    num_time_test<- length(test_index)

    pred_y<- array(0, dim=c(num_time_test, k, n))

    pred_x<- array(0, dim=c(num_time_test, p, n))

  num_time_all<- (q + num_time_test)

  # if do not need to predict: number of time points in latent y is the same as that in MCEM, which is q
  # if need to predict: number of time points in latent y should be q + num_time_test

  latent_y<- array(0, dim = c(num_time_all, k, n))

  latent_y[(1:q),,]<- y_init[chain_index,,,] # q*k*n array, initials generated by estimated covariance matrix

  ###############################################################################################################################################
  ################################# GP property to obtain predicted y based on estimated y and GP parameters ####################################
  ###############################################################################################################################################

  if(model_dgp){

    final_gp_hp<- hyper_record[index_used,]

    h3n2_data_all<- list()
    list_temp <- vector("list", k)
    for (list_index in 1:k){
      list_temp[[list_index]]<- a_full
    }
    h3n2_data_all$input<- list_temp

    # cov_all consists of all time points: time points with observed gene expression and without (q + num_time_test)
    cov_all<- GPFDA::mgpCovMat(Data=h3n2_data_all, hp=final_gp_hp)

    cov_all<- cov2cor(cov_all)

  } else {

    final_gp_hp<- hyper_record[index_used,,]

    cov_all<- matrix(0,nrow=(k*num_time_all),ncol=(k*num_time_all))

    h3n2_data_full_generate<- list()
    h3n2_data_full_generate$input<- list(a_full)

    for (j in 1:k){
      # for the jth process
      cov_all[((j-1)*num_time_all+1):(j*num_time_all),((j-1)*num_time_all+1):(j*num_time_all)]<- mgpCovMat(Data=h3n2_data_full_generate,
                                                                                                           hp=final_gp_hp[j,]) # hps for the jth process
    }

    cov_all<- cov2cor(cov_all)

  }

  cov_all_inv<- solve(cov_all)

    missing_time_num<- rep(0, times =  n)
    missing_time_index <- vector("list",  n)

    prod_covnew_covestinv<- vector("list",  n)
    cov_cond_dist<- vector("list",  n)
    sigmay_inv<- vector("list", n)

    # a_person: subject-specific time points in the training data
    # a_full: complete time points in the training data + new time points in the test data

    for (list_index in 1:n){

      result<- subject_specific_objects(k, num_time_all, a_full, a_person[[list_index]], cov_all)

      # results returned from the algorithm
      missing_time_num[list_index]<- result$missing_time_num
      missing_time_index[[list_index]]<- result$missing_time_index
      prod_covnew_covestinv[[list_index]]<- result$prod_covnew_covestinv
      cov_cond_dist[[list_index]]<- result$cov_cond_dist
      sigmay_inv[[list_index]]<- result$cov_est_inv

    }

  ##################################################################################################################################
  ################################################################# Constants ######################################################
  ##################################################################################################################################

  big_z_table<- as.matrix(table_generator(k))

  # Rcpp::sourceCpp('~/bsfadgp/src/gibbs_after_mcem_irregular_time.cpp')

  directory_name = paste0(pathname, "/gibbs_after_mcem_chain_", chain_index)
  dir.create(directory_name)
  setwd(directory_name)

  c0<- ig_parameter
  d0<- ig_parameter
  c1<- ig_parameter
  d1<- ig_parameter
  c2<- ig_parameter
  d2<- ig_parameter

  e0<- prior_sparsity*p
  f0<- (1-prior_sparsity)*p

  #######################################################################################################################################################
  ########################################### Final Gibbs Sampler ###########################################################################
  #######################################################################################################################################################

  out <- gibbs_after_mcem_irregular_time(latent_y,
                                         x,
                                         big_a,
                                         big_z,
                                         phi,
                                         pai,
                                         beta,
                                         k, n, p, q, num_time_test, c0, c1, d0, d1, e0, f0, 1, mc_num,
                                          TRUE, pred_y, pred_x, big_z_table, ind_x, individual_mean, mu_g, variance_g, c2, d2,
                                          obs_time_num,
                                          obs_time_index,
                                          missing_time_num,
                                          missing_time_index,
                                          prod_covnew_covestinv,
                                          cov_cond_dist,
                                          sigmay_inv,
                                          thin_step,
                                          burnin)

  # print("Gibbs Sampler Finished.")

}
