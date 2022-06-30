#-------------------------------------------------------------------------
#Packages
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)
library(sandwich)
library(here)

#Set path to main folder
main_file_path <- "..."

ref_tables <- paste0(main_file_path,'/raw_data')
libin <- paste0(main_file_path,'/processed_data')
libout <- paste0(main_file_path,'/final_data')

#-------------------------------------------------------------------------
#Utility Functions and External Macors
any_in <- function(egg, nest){
  ifelse(sum(egg %in% nest) > 0, 1, 0)
  }

#Misc Functions
source(here("code/ackley_funs.R")) #paste0(main_file_path,'/macros/ackley_funs.R'))

#This is a modification of the add_pi function in ciTools which uses robust standard errors
#See https://cran.r-project.org/web/packages/ciTools/vignettes/ciTools-glm-vignette.html
source(here("code/sim_glm_robust_fun.R")) #paste0(main_file_path,'/macros/sim_glm_robust_fun.R')) 


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#Chunk 0: Load input data and conduct minor edits

#-------------------------------------------------------------------------
#Import analysis data
#setwd(libin)

dat <- read.csv(here('raw_data/acm_2011-2019.csv'))
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

base_year <- 1998

dat_edit <- dat %>% filter(year >= 2011, !is.na(county)) %>%
  mutate(time = year - base_year, #Normalize time to 1999 = 1
         time_orig_vals = time)

#dat_est <- dat_edit %>% filter(year < 2020) NOTE: maybe i need 2022 in this dataset

#take the mean population over the 8 yrs represented 
dat_500 <- dat_edit %>% 
  group_by(county_code) %>% 
  summarize(mean_pop = mean(population)) %>% 
  ungroup()

#subset the top 500 most populous counties 
dat_500 <- dat_500[order(dat_500$mean_pop, decreasing = TRUE),]
dat_500 <- dat_500[1:500,]
dat_edit <- subset(dat_edit, county_code %in% dat_500$county_code)

#check to see that 500 unique counties are represented
num_counties_dat_500 <- length(unique(dat_edit$county_code))
num_counties_dat_500

#arrange/order counties by name and year
dat_est <- dat_edit %>% arrange(county, year)

#dat_2020 <- dat_edit %>% filter(year == 2020)

fit <- glm(deaths ~ death_rate_lag1 +
             time +
             time*factor(county_code),
           family = quasipoisson(link = "log"),
           data = dat_est)


# if(1==1){
#   setwd(libin)
#   save(fit, file = 'estimated_poisson_glm_parameters.rda')
# }

#Compute cluster-robust errors
vmat_rob <- vcovCL(fit, cluster = ~ county_code)
clust_ses <- sqrt(diag(vmat_rob)) #produces warning: NAs produced 

# if(1==1){
#   save(clust_ses, file = 'estimated_poisson_glm_ses.rda')
# }

#Add fitted values for all years and predicted values for 2022
dat_est <- dat_est %>%
  mutate(fitted_deaths_all_yrs = exp(predict.glm(fit, dat_est)),
         #fitted_death_rate_all_yrs = fitted_deaths_all_yrs/death_offset
         )

#read in 2022 provisional data 
dat_2021_2022 <- read_csv("raw_data/acm_2021-2022_provisional.csv")

#set time var
dat_2021_2022_edit <- dat_2021_2022 %>%
  filter(year >= 2011, !is.na(county)) %>%
  mutate(time = year - base_year, #Normalize time to 1999 = 1
         time_orig_vals = time)

#subset to include the 500 most populous counties
dat_2021_2022_edit <- subset(dat_2021_2022_edit, county_code %in% dat_500$county_code)

#check to see that 500 unique counties are represented
num_counties_dat_500_2022 <- length(unique(dat_2021_2022_edit$county_code))
num_counties_dat_500_2022

dat_2021_2022 <- dat_2021_2022_edit %>% 
  mutate(fitted_deaths_2020 = exp(predict.glm(fit, dat_2021_2022_edit)),
         #fitted_death_rate_2020 = fitted_deaths_2020/death_offset
         )

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#Chunk 2: Estimate prediction 95% interval for all years 
#-------------------------------------------------------------------------
library(arm)
#Set parameters
set.seed(79)
mod_mat <- model.matrix(fit)
alpha <- .05
nsims <- 1000
npreds <- nrow(dat_est)
overdisp <- summary(fit)$dispersion
sim_response_mat <- matrix(NA, ncol = nsims, nrow = npreds)

#Sample from coef distributions
#This is a modification of the add_pi function in ciTools which uses robust standard errors
#See https://cran.r-project.org/web/packages/ciTools/vignettes/ciTools-glm-vignette.html for complete details
sim_coefs <- sim_glm_robust(fit, n.sims = nsims, ses = clust_ses)
#save(sim_coefs, file = 'estimated_poisson_glm_simulated_coefs.rda')


for(i in 1:nsims){
  #Fitted value with new coef draw
  yhat <- dat_est$death_offset * exp (mod_mat %*% sim_coefs@coef[i,])
  
  disp_parm <- yhat / (overdisp - 1) #Set new dispersion parameter
  
  #Draw new death count and fill in matrix. Each row is now a sample of size nsims for each county-year. Each col is a draw.
  sim_response_mat[,i] <- rnegbin(n = npreds,
                                  mu = yhat,
                                  theta = disp_parm)
}



#Gather statistics from simulated distributions
sds <- sqrt(apply(sim_response_mat,1,var))
lwr <- apply(sim_response_mat, 1, FUN = quantile, probs = alpha/2)
upr <- apply(sim_response_mat, 1, FUN = quantile, probs = 1 - alpha / 2)

#Add the computed predictions intervals to main table
dat_est <- dat_est %>% 
  mutate(pred_deaths_lwr_ci = lwr,
         pred_deaths_upr_ci = upr,
         pred_death_rate_lwr_ci = lwr/death_offset,
         pred_death_rate_upr_ci = upr/death_offset,
         pred_death_std_err = sds,
         pred_death_rate_std_err = sds/death_offset
         
  )


#-------------------------------------------------------------------------

#Produce plot to make sure eveything looks fine
check_dat <- dat_est %>% filter(cs_code == '01CS031')

ggplot(check_dat, aes(x = year, y = total_deaths)) +
  ggtitle("Quasipoisson Regression", subtitle = "Model fit (black line), with Prediction intervals (gray), Confidence intervals (dark gray)") +
  geom_point(size = 2) +
  geom_line(aes(x = year, y = fitted_deaths_all_yrs), size = 1.2) +
  geom_ribbon(aes(ymin = pred_deaths_lwr_ci , ymax = pred_deaths_upr_ci), alpha = 0.2)
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Estimate prediction interval for 2022. This repeats everything from prior chunk just using 2020 data
#Need this just to get model matrix for 2020
dat_2022 <- semi_join(dat_2022, dat_est, by = 'cs_code')
dat_est2 <- dat_est %>% filter(time == 19) %>% bind_rows(dat_2022)


fit2 <- glm(total_deaths ~ offset(log(death_offset)) + death_rate_lag1 +
              time +
              time*factor(cs_code) , family = quasipoisson(link = "log") ,data = dat_est2)


mod_mat <- model.matrix(fit2)
mod_mat <- mod_mat[which(mod_mat[,'time'] %in% c(22)),]
npreds <- nrow(dat_2020)
sim_response_mat <- matrix(NA, ncol = nsims, nrow = npreds)


for(i in 1:nsims){
  #Fitted value with new coef draw
  yhat <- dat_2020$death_offset* exp (mod_mat %*% sim_coefs@coef[i,])
  
  disp_parm <- yhat / (overdisp - 1) #Set new dispersion parameter
  
  #Draw new death count and fill in matrix. Each row is now a sample of size nsims for each county-year. Each col is a draw.
  sim_response_mat[,i] <- rnegbin(n = npreds,
                                  mu = yhat,
                                  theta = disp_parm)
}


if(1==1){
  boot_samp_out <- bind_cols(as.data.frame(sim_response_mat),dat_2020[,c('cs_code','death_rate','death_offset')]) %>% 
    rename(all_cause_death_rate_2020 = death_rate, death_offset_2020 = death_offset)
  saveRDS(boot_samp_out, file = 'full_bootstrap_sample_2020.rds')
}

sds <- sqrt(apply(sim_response_mat,1,var))
lwr <- apply(sim_response_mat, 1, FUN = quantile, probs = alpha/2)
upr <- apply(sim_response_mat, 1, FUN = quantile, probs = 1 - alpha / 2)

dat_2020 <- dat_2020 %>% 
  mutate(pred_deaths_lwr_ci_2020 = lwr,
         pred_deaths_upr_ci_2020 = upr,
         pred_death_rate_lwr_ci_2020 = lwr/death_offset,
         pred_death_rate_upr_ci_2020 = upr/death_offset,
         pred_death_std_err = sds,
         pred_death_rate_std_err = sds/death_offset
  )

#NOTE MUST UNLOAD THESE LIBRARIES OR CODE BELOW MAY CRASH. THERE IS A CONFLICT WITH THE USE OF SELECT IN MASS PACKAGE AND DPLYR!
#SOMETIMES THIS UNLOAD FAILS IN A CONTINUOUS RUN OF THIS FILE
unloadNamespace("arm")
unloadNamespace("lme4")
unloadNamespace("MASS")
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#Chunk 3: Construct Final Output Data for 2020
#-------------------------------------------------------------------------


dat_2020_out <- dat_2020 %>% 
  rename(
    covid_deaths_2020 = covid_deaths,
    death_offset_2020 = death_offset,
    all_cause_deaths_2020 = total_deaths,
    all_cause_death_rate_2020 = death_rate
  ) %>% 
  mutate(
    covid_death_rate_2020 = covid_deaths_2020/death_offset_2020,
    excess_deaths_2020 = all_cause_deaths_2020 - fitted_deaths_2020,
    excess_assigned_covid_deaths = covid_deaths_2020,
    excess_unassigned_covid_deaths = excess_deaths_2020 - covid_deaths_2020,
    prop_deaths_unassigned = excess_unassigned_covid_deaths/excess_deaths_2020,
    prop_deaths_unassigned = ifelse(prop_deaths_unassigned < 0 | excess_unassigned_covid_deaths < 0,NA, prop_deaths_unassigned),
    excess_death_rate_2020 = all_cause_death_rate_2020 - fitted_death_rate_2020,
    excess_assigned_covid_death_rate = covid_death_rate_2020,
    excess_unassigned_covid_death_rate = excess_death_rate_2020 - covid_death_rate_2020,
    prop_death_rate_unassigned = excess_unassigned_covid_death_rate/excess_death_rate_2020,
    prop_death_rate_unassigned = ifelse(prop_death_rate_unassigned < 0 | excess_unassigned_covid_death_rate < 0,NA, prop_death_rate_unassigned)
  )


setwd(libout)
write.csv(dat_2020_out, file = 'fitted_and_actual_deaths_county_sets_2020_W2020_wash_6_3.csv')


dat_newrows_2020 <- dat_2020 %>% 
  rename(
    fitted_deaths_all_yrs = fitted_deaths_2020,
    fitted_death_rate_all_yrs = fitted_death_rate_2020,
    pred_deaths_upr_ci = pred_deaths_upr_ci_2020,
    pred_deaths_lwr_ci = pred_deaths_lwr_ci_2020,
    pred_death_rate_upr_ci = pred_death_rate_upr_ci_2020,
    pred_death_rate_lwr_ci = pred_death_rate_lwr_ci_2020)

dat_pan <- dat_est %>%
  bind_rows(.,dat_newrows_2020) %>% arrange(cs_code, time)


dat_pan_out <- dat_pan %>% 
  rename(
    total_death_rate = death_rate,
    fitted_deaths = fitted_deaths_all_yrs,
    fitted_death_rate = fitted_death_rate_all_yrs,
    population = cens_pop_est) %>%
  mutate(death_offset = (population)/1000)


#Output primary panel data from 2011-2020 used to make tables and figures
setwd(libout)
write.csv(dat_pan_out, file = 'fitted_and_actual_deaths_county_sets_2011_2020_W2020_wash_6_3.csv', row.names = FALSE)



