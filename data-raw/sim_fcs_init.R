## code to prepare `sim_fcs_init` dataset goes here

############################################# data preparation #############################################
setwd("~/rds/hpc-work/corrected_code/sim/truth_correlated_factors/generated_var_not_1/fitting_4_factors")

k<- 4
p<- 100
q<- 8
n<- 17

regression_matrix<- read.delim("mA.txt",header = FALSE) # this is only a, need to combine this with mPostPib.txt (z) to obtain the final factor loading matrix
regression_matrix<- regression_matrix[,-(k+1)]
regression_matrix<- as.matrix(regression_matrix)

cont_matrix<- read.delim("mPostPib.txt",header = FALSE)
cont_matrix<- cont_matrix[,-(k+1)]
cont_matrix<- as.matrix(cont_matrix)

binary_matrix<- cont_matrix

# turn the above posterior probability matrix into a binary matrix z
for (i in 1:nrow(binary_matrix)){
  for (j in 1:ncol(binary_matrix)){
    #if (cont_matrix[i,j]>0.5){
    if (cont_matrix[i,j]>0.99){
      binary_matrix[i,j]<- 1
    } else {
      binary_matrix[i,j]<- 0
    }
  }
}

# latent factor scores
latent_factor_score<- read.delim("mF.txt",header = FALSE)
latent_factor_score<- latent_factor_score[,-((n*q)+1)]

# regression variance
regression_variance<- read.delim("mPsi.txt",header = FALSE) # variance in the regression model
regression_variance<- regression_variance[,-(p+1)]
regression_variance<- as.numeric(regression_variance)

# the first set of initials: a and z set as estimation from BFRM
y_init = latent_factor_score
phi_init = regression_variance
a_init_1 = regression_matrix
z_init_1 = binary_matrix

# the second set of initials: a and z set as truth
a_init_2 = sim_fcs_truth$real_a
z_init_2 = sim_fcs_truth$real_z

################################################ import results for irregular data from BFRM #################################################################################
setwd("/home/jc2312/use_all_fitting/")

regression_matrix<- read.delim("mA.txt",header = FALSE) # this is only a, need to combine this with mPostPib.txt (z) to obtain the final factor loading matrix
regression_matrix<- regression_matrix[,-(k+1)]
regression_matrix<- as.matrix(regression_matrix)

cont_matrix<- read.delim("mPostPib.txt",header = FALSE)
cont_matrix<- cont_matrix[,-(k+1)]
cont_matrix<- as.matrix(cont_matrix)

binary_matrix<- cont_matrix

# turn the above posterior probability matrix into a binary matrix z
for (i in 1:nrow(binary_matrix)){
  for (j in 1:ncol(binary_matrix)){
    #if (cont_matrix[i,j]>0.5){
    if (cont_matrix[i,j]>0.99){
      binary_matrix[i,j]<- 1
    } else {
      binary_matrix[i,j]<- 0
    }
  }
}

# latent factor scores
latent_factor_score<- read.delim("mF.txt", header = FALSE)
latent_factor_score<- latent_factor_score[,-(sum(obs_time_num)+1)]

# latent_factor_score<- latent_factor_score[,-((n*q_common)+1)]

# regression variance
regression_variance<- read.delim("mPsi.txt",header = FALSE) # variance in the regression model
regression_variance<- regression_variance[,-(p+1)]
regression_variance<- as.numeric(regression_variance)

# irregular data: use all observations
y_init_irregular = latent_factor_score
phi_init_irregular = regression_variance
a_init_irregular = regression_matrix
z_init_irregular = binary_matrix

################################################### import results for six common time points from BFRM ###############################################################
setwd("/home/jc2312/use_common_fitting/")

q_common<- 6
k<- 4
p<- 100
n<- 17

regression_matrix<- read.delim("mA.txt",header = FALSE) # this is only a, need to combine this with mPostPib.txt (z) to obtain the final factor loading matrix
regression_matrix<- regression_matrix[,-(k+1)]
regression_matrix<- as.matrix(regression_matrix)

cont_matrix<- read.delim("mPostPib.txt",header = FALSE)
cont_matrix<- cont_matrix[,-(k+1)]
cont_matrix<- as.matrix(cont_matrix)

binary_matrix<- cont_matrix

# turn the above posterior probability matrix into a binary matrix z
for (i in 1:nrow(binary_matrix)){
  for (j in 1:ncol(binary_matrix)){
    #if (cont_matrix[i,j]>0.5){
    if (cont_matrix[i,j]>0.99){
      binary_matrix[i,j]<- 1
    } else {
      binary_matrix[i,j]<- 0
    }
  }
}

# latent factor scores
latent_factor_score<- read.delim("mF.txt", header = FALSE)
latent_factor_score<- latent_factor_score[,-((n*q_common)+1)]

# regression variance
regression_variance<- read.delim("mPsi.txt",header = FALSE) # variance in the regression model
regression_variance<- regression_variance[,-(p+1)]
regression_variance<- as.numeric(regression_variance)

# irregular data: use all observations
y_init_six_common = latent_factor_score
phi_init_six_common = regression_variance
a_init_six_common = regression_matrix
z_init_six_common = binary_matrix

# store all data want to be used further in a single R object - the name should be the same as the R file
sim_fcs_init<- list(y_init = y_init,
                    phi_init = phi_init,
                    a_init_1 = a_init_1,
                    z_init_1 = z_init_1,
                    a_init_2 = a_init_2,
                    z_init_2 = z_init_2,
                    y_init_irregular =  y_init_irregular,
                    phi_init_irregular = phi_init_irregular,
                    a_init_irregular = a_init_irregular,
                    z_init_irregular = z_init_irregular,
                    y_init_six_common = y_init_six_common,
                    phi_init_six_common = phi_init_six_common,
                    a_init_six_common = a_init_six_common,
                    z_init_six_common = z_init_six_common)

# add prepared data to package
usethis::use_data(sim_fcs_init, overwrite = TRUE)
