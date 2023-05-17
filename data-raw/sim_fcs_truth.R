## code to prepare `sim_fcs_truth` dataset goes here

### set seed for reproducbility
set.seed(456)

### constants
p<- 100
k<- 4
n<- 17

# training data + test data
q<- 8
num_time_test<- 2
num_time_all<- (q + num_time_test)

sparsity_e<- 0.1*p
sparsity_f<- 0.9*p

a_full<- c(0.0, 5.0, 12.0, 45.5, 53.0, 60.0, 69.5, 77.0, 84.0, 93.5)
a_train<- a_full[(1:q)]
a_test<- a_full[(q+1):(q+num_time_test)]

# for all data
h3n2_data_all<- list()
list_temp <- vector("list", k)
for (list_index in 1:k){
  list_temp[[list_index]]<- a_full
}
h3n2_data_all$input<- list_temp

# for training data
h3n2_data<- list()
list_temp <- vector("list", k)
for (list_index in 1:k){
  list_temp[[list_index]]<- a_train
}
h3n2_data$input<- list_temp

# for testing data
h3n2_data_test<- list()
list_temp <- vector("list", k)
for (list_index in 1:k){
  list_temp[[list_index]]<- a_test
}
h3n2_data_test$input<- list_temp

# load the source of covariance matrix
load("~/rds/hpc-work/corrected_code/sim/truth_correlated_factors/generated_var_is_1/generated_real_data.RData")

### generate real data
real_y<- array(0,dim=c(num_time_all, k, n)) # time_all*k*n

# for small sd + correlated factors (already tried)
sigma_y_init_truth<- GPFDA::mgpCovMat(Data=h3n2_data_all, hp=final_gp_hp_truth)

for (n_index in 1:n){
  temp_vec<- mvtnorm::rmvnorm(1, mean=rep(0,times=(num_time_all*k)), sigma= sigma_y_init_truth)
  real_y[,,n_index]<- temp_vec
}

real_a<- matrix(rnorm((p*k), mean = 4, sd = 1), nrow=p, ncol=k) # loading ranges from 2 to 5 (if loaded on one factor)

pai_init<- rbeta(1, sparsity_e, sparsity_f) # for the same chain, initial values for different pi is the same

real_z<- matrix(0, nrow = p, ncol = k)

for (k_index in 1:k){
  real_z[,k_index]<- as.integer(purrr::rbernoulli(p, pai_init))
}

real_l<- real_a*real_z

# subject-gene mean
gene_mean<- runif(p, 4, 16) # gene's mean ranges from 4~16
gene_person_mean<- matrix(0, nrow = p, ncol =n)

for (gene_index in 1:p){
  gene_person_mean[gene_index, ]<- rnorm(n, mean = gene_mean[gene_index], sd = 0.5)
}

# random error: noise ranges from -1 to 1
real_e<- array(0, dim = c(p,num_time_all,n))
real_e_sd<- 0.5

for (person_index in 1:n){
  real_e[,,person_index]<-
    matrix(rnorm(n=(p*num_time_all), mean = 0, sd=real_e_sd),
           nrow= p,
           ncol= num_time_all)
}

# construct x: for each person at each time point, the expression of biomarker is the sum of subject-gene mean, y's expression, and random noise
observed_x<- array(0, dim = c(p,num_time_all,n))

for (person_index in 1:n){
  for (time_index in 1:num_time_all){
     observed_x[,time_index,person_index]<-
      gene_person_mean[,person_index]+
      (real_l%*%(real_y[time_index,,person_index])) +
      real_e[,time_index,person_index]
  }
}

# all data split into training data and test data
observed_x_train<- observed_x[,(1:q),] # p*q*n
observed_x_pred<- observed_x[,((q+1):(q+num_time_test)),]

observed_x_pred_reformat<- array(0, dim = c(n, p, num_time_test))
for (person_index in 1:n){
  observed_x_pred_reformat[person_index,,]<- observed_x_pred[,,person_index]
}

# all data as training data
observed_x_train<- observed_x

############################ irregular time points: if succeed, save this part as well ################################
### 1. regular missing - for people with odd number, only have 6 observed time points
# create irregular observations
observed_x_train_irregular<- list()

obs_time_num<- rep(0, times = n)

obs_time_index<- list()

a_person<- list()

for (person_index in 1:n){
  if (person_index %%2 ==1){
    observed_x_train_irregular[[person_index]]<- observed_x_train[,-c(3,6),person_index]
    a_person[[person_index]]<- a_train[-c(3,6)]
    obs_time_index[[person_index]]<- c(1,2,4,5,7,8)
  } else {
    observed_x_train_irregular[[person_index]]<- observed_x_train[,,person_index]
    a_person[[person_index]]<- a_train
    obs_time_index[[person_index]]<- (1:8)
  }

  obs_time_num[person_index]<- ncol(observed_x_train_irregular[[person_index]])

}

# reformat the observations into a matrix ready to be input to BFRM software: p*sum(n_i) matrix for returning a factor score matrix of same dimension
observed_x_train_irregular_bfrm<- matrix(0, nrow = p, ncol = sum(obs_time_num))

observed_x_train_irregular_bfrm_center_within_individual<- matrix(0, nrow = p, ncol = sum(obs_time_num))

gene_individual_mean_est<- matrix(0, nrow = p, ncol = n)

col_person_index<- list()

for (gene_index in 1:p){

  col_start_index<- 0

  for (person_index in 1:n){

  observed_x_train_irregular_bfrm[gene_index, ((col_start_index+1):(col_start_index + obs_time_num[person_index]))]<- observed_x_train_irregular[[person_index]][gene_index,] # n_i dimensional vector

  # center within individual

  gene_individual_mean_est[gene_index, person_index]<- mean(observed_x_train_irregular_bfrm[gene_index, ((col_start_index+1):(col_start_index + obs_time_num[person_index]))])

  observed_x_train_irregular_bfrm_center_within_individual[gene_index, ((col_start_index+1):(col_start_index + obs_time_num[person_index]))]<-
    observed_x_train_irregular_bfrm[gene_index, ((col_start_index+1):(col_start_index + obs_time_num[person_index]))] - gene_individual_mean_est[gene_index, person_index]

  if (gene_index == 1){
    col_person_index[[person_index]]<- ((col_start_index+1):(col_start_index + obs_time_num[person_index]))
  }

  # restarting index for the next person
  col_start_index<- col_start_index +  obs_time_num[person_index]

  }
}

# store all data want to be used further in a single R object - the name should be the same as the R file
sim_fcs_truth<- list(a_full = a_full,
                     gp_hp_truth = final_gp_hp_truth,
                     gp_sigmay_truth = sigma_y_init_truth,
                     sparsity_e_truth = sparsity_e,
                     sparsity_f_truth = sparsity_f,
                     real_y = real_y,
                     real_a = real_a,
                     real_z = real_z,
                     gene_mean = gene_mean, gene_person_mean = gene_person_mean,
                     observed_x = observed_x,
                     observed_x_train = observed_x_train,
                     observed_x_pred = observed_x_pred,
                     observed_x_pred_reformat = observed_x_pred_reformat,
                     observed_x_train_irregular = observed_x_train_irregular,
                     obs_time_num = obs_time_num,
                     obs_time_index = obs_time_index,
                     a_person = a_person,
                     col_person_index = col_person_index)

# add prepared data to package
usethis::use_data(sim_fcs_truth, overwrite = TRUE)
