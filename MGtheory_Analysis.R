library(tidyverse)
#Set directory to location where data and source file are saved.
source("mGtheory_source.R")
#Read in data. Data is in long format
text_dat1 <- read.csv("text_dat_long.csv", header = T)

###### G STUDY
PARE_g <- calcMGtheory_Gstudy(text_dat1$id, text_dat1$rater, text_dat1$scale, text_dat1$TextQual)
PARE_g


######
######
###### D STUDY

###Get D study number of raters
n_r_prime <- 4
###Get covariance matrices from G study
sigma_p <- PARE_g[[1]][[1]]
sigma_r <- PARE_g[[2]][[1]]
sigma_pr <- PARE_g[[3]][[1]]
###Compute absolute error covariance components
delta_pI <- sigma_r/n_r_prime + sigma_pr/n_r_prime
###Define weights for composite
#weights1 <- rep(1, times = 8)
weights_vec <- c(1, 1, 1, 1, 1, 1, 0.5, 0.5)
weights <- weights_vec/(sum(weights_vec))

##Compute weights matrix
w_mat <- CalcWeightMat(weights)
##Composite reliability
sum(w_mat*sigma_p)/sum(w_mat*(sigma_p + delta_pI))


