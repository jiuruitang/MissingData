pop_dat_small <- gen_pop_Poisson(N = 200)
sub_dat_small <- sub_dat <- getSample_Poisson(population = pop_dat_small, N=200)

sampleTotal <- replicate(200,getSampleStatsPoisson(population = pop_dat_small,N=200))
pop_sd_HT <- sd(sampleTotal)
pop_mean_HT <- mean(sampleTotal)

sub_dat_obs <- sub_dat_small %>% filter(Rx == 0)
sub_dat_mis <- sub_dat_small %>% filter(Rx == 1)

# total of observed X1
sum(sub_dat_obs$X1*sub_dat_obs$W)

n_mis = dim(sub_dat_mis)[1]
mis_comb <- expand.grid(replicate(n_mis, 0:1, simplify = FALSE))
mis_comb_mat <- as.matrix(mis_comb)
mis_total <- mis_comb_mat %*% sub_dat_mis$W
total_comb <- mis_total + sum(sub_dat_obs$X1*sub_dat_obs$W)
hist(total_comb)
abline(v=sum(pop_dat_small$X1),col="red")
abline(v=sum(pop_dat_small$X1)-3*pop_sd_HT, col = "green")

HT_var <- sum((sub_dat_small$X1)^2/(sub_dat_small$pi)^2*(1-sub_dat_small$pi))

