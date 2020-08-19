library(mgcv);
library(refund);
library(boot);
source("simulate_functional_mediation_example.R");
source("funreg_mediation.R");
source("tvem.R");
source("select_tvem.R");
source("plot_funreg_mediation.R");
source("print_funreg_mediation.R");
source("plot_tvem.R");
source("print_tvem.R");
set.seed(the_seed);
nboot <- 99;
answers <- NULL;
start_time <- Sys.time(); 
for (this_sim in 1:nsim) { 
  simulated_data <- simulate_functional_mediation_example(nsub=nsub,
                                                          simulate_binary_Y=TRUE);
  true_indirect <- simulated_data$true_indirect;
  true_alpha_int <- simulated_data$true_alpha_int;
  true_alpha_X <- simulated_data$true_alpha_X;
  true_beta_int <- simulated_data$true_beta_int;
  true_beta_X <- simulated_data$true_beta_X;
  true_beta_M <- simulated_data$true_beta_M;
  long_data <- simulated_data$dataset;
  fun_med_results <- funreg_mediation(data=long_data,
                                      treatment=X,
                                      mediator=M,
                                      outcome=Y,
                                      id=subject_id,
                                      time=t,
                                      logistic=TRUE,
                                      nboot=nboot);
  tvem_model_summary <- summary(fun_med_results$original_results$tvem_XM_details$back_end_model);
  alpha_int_residuals <- fun_med_results$original_results$alpha_int_estimate - true_alpha_int;
  alpha_int_se <- fun_med_results$original_results$alpha_int_se;
  alpha_X_residuals <- fun_med_results$original_results$alpha_X_estimate - true_alpha_X;
  alpha_X_se <- fun_med_results$original_results$alpha_X_se;
  beta_M_residuals <- fun_med_results$original_results$beta_M_estimate - true_beta_M;
  beta_M_se <- fun_med_results$original_results$beta_M_se;
  norm_cover <- (fun_med_results$bootstrap_results$indirect_effect_boot_norm_lower < true_indirect) & 
    (fun_med_results$bootstrap_results$indirect_effect_boot_norm_upper > true_indirect);
  norm_power <- (sign(fun_med_results$bootstrap_results$indirect_effect_boot_norm_lower)==
                   sign(fun_med_results$bootstrap_results$indirect_effect_boot_norm_upper));
  basic_cover <- (fun_med_results$bootstrap_results$indirect_effect_boot_basic_lower < true_indirect) & 
    (fun_med_results$bootstrap_results$indirect_effect_boot_basic_upper > true_indirect);
  basic_power <- (sign(fun_med_results$bootstrap_results$indirect_effect_boot_basic_lower)==
                    sign(fun_med_results$bootstrap_results$indirect_effect_boot_basic_upper));
  bca_cover <- (fun_med_results$bootstrap_results$indirect_effect_boot_bca_lower < true_indirect) & 
    (fun_med_results$bootstrap_results$indirect_effect_boot_bca_upper > true_indirect);
  bca_power <- (sign(fun_med_results$bootstrap_results$indirect_effect_boot_bca_lower)==
                  sign(fun_med_results$bootstrap_results$indirect_effect_boot_bca_upper));
  perc_cover <- (fun_med_results$bootstrap_results$indirect_effect_boot_perc_lower < true_indirect) & 
    (fun_med_results$bootstrap_results$indirect_effect_boot_perc_upper > true_indirect);
  perc_power <- (sign(fun_med_results$bootstrap_results$indirect_effect_boot_perc_lower)==
                   sign(fun_med_results$bootstrap_results$indirect_effect_boot_perc_upper));
  current_time <- Sys.time();
  answers <- rbind(answers, 
                   # Create structures to hold results:  
                   # ... for the effect of X on M: 
                   c(nsub=nsub, 
                     this_sim=this_sim,
                     alpha_int_estimate_bias =  mean(alpha_int_residuals), 
                     alpha_int_estimate_mse =  mean(alpha_int_residuals^2),
                     alpha_int_estimate_mean_std_err = mean(alpha_int_se),
                     alpha_int_estimate_coverage =  mean(abs(alpha_int_residuals)<1.96*alpha_int_se), 
                     alpha_int_estimate_family_coverage =  all(abs(alpha_int_residuals)<1.96*alpha_int_se), 
                     alpha_int_estimate_pvalue_overall =  as.numeric(tvem_model_summary$p.coeff["(Intercept)"]), 
                     alpha_int_estimate_pvalue_varying =  tvem_model_summary$s.table["s(time)","p-value"],
                     alpha_X_estimate_bias =  mean(alpha_X_residuals), 
                     alpha_X_estimate_mse =  mean(alpha_X_residuals^2),
                     alpha_X_estimate_mean_std_err = mean(alpha_X_se),
                     alpha_X_estimate_coverage =  mean(abs(alpha_X_residuals)<1.96*alpha_X_se), 
                     alpha_X_estimate_family_coverage =  all(abs(alpha_X_residuals)<1.96*alpha_X_se), 
                     alpha_X_estimate_pvalue_overall =  as.numeric(tvem_model_summary$p.coeff["(Intercept)"]),   
                     alpha_X_estimate_pvalue_varying =  tvem_model_summary$s.table["s(time):treatment","p-value"],  
                     # ... for the joint effect of X and M on Y: 
                     beta_int_bias =  fun_med_results$original_results$beta_int_estimate - true_beta_int,
                     beta_int_mse = (fun_med_results$original_results$beta_int_estimate - true_beta_int)^2,
                     beta_int_se =  fun_med_results$original_results$beta_int_se,
                     beta_int_coverage = abs(fun_med_results$original_results$beta_int_estimate - true_beta_int) < 1.96*fun_med_results$original_results$beta_int_se,
                     beta_X_bias =  fun_med_results$original_results$beta_X_estimate - true_beta_X,
                     beta_X_mse = (fun_med_results$original_results$beta_X_estimate - true_beta_X)^2,
                     beta_X_se =  fun_med_results$original_results$beta_X_se,
                     beta_X_coverage = abs(fun_med_results$original_results$beta_X_estimate - true_beta_X) < 1.96*fun_med_results$original_results$beta_X_se,
                     beta_M_estimate_bias =  mean(beta_M_residuals), 
                     beta_M_estimate_mse =  mean(beta_M_residuals^2), 
                     beta_M_estimate_coverage =  mean(abs(beta_M_residuals)<1.96*beta_M_se), 
                     beta_M_estimate_family_coverage =  all(abs(beta_M_residuals)<1.96*beta_M_se), 
                     beta_M_pvalue =   fun_med_results$original_results$beta_M_pvalue,  
                     # ... for the direct effect of X on Y: 
                     tau_int_estimate =  fun_med_results$original_results$tau_int_estimate, 
                     tau_int_se =  fun_med_results$original_results$tau_int_se, 
                     tau_X_estimate =  fun_med_results$original_results$tau_X_estimate, 
                     tau_X_se =  fun_med_results$original_results$tau_X_se,  
                     # ... for the mediation effect of X on Y through M:  
                     indirect_effect_bias = fun_med_results$bootstrap_results$indirect_effect_boot_estimate - true_indirect,
                     indirect_effect_mse = (fun_med_results$bootstrap_results$indirect_effect_boot_estimate - true_indirect)^2,
                     indirect_effect_boot_se = fun_med_results$bootstrap_results$indirect_effect_boot_se, 
                     indirect_effect_boot_norm_lower = fun_med_results$bootstrap_results$indirect_effect_boot_norm_lower, 
                     indirect_effect_boot_norm_upper = fun_med_results$bootstrap_results$indirect_effect_boot_norm_upper, 
                     indirect_effect_boot_norm_coverage = norm_cover,
                     indirect_effect_boot_norm_power = norm_power,
                     indirect_effect_boot_basic_lower = fun_med_results$bootstrap_results$indirect_effect_boot_basic_lower, 
                     indirect_effect_boot_basic_upper = fun_med_results$bootstrap_results$indirect_effect_boot_basic_upper, 
                     indirect_effect_boot_basic_coverage = basic_cover, 
                     indirect_effect_boot_basic_power = basic_power,
                     indirect_effect_boot_bca_lower = fun_med_results$bootstrap_results$indirect_effect_boot_bca_lower, 
                     indirect_effect_boot_bca_upper = fun_med_results$bootstrap_results$indirect_effect_boot_bca_upper, 
                     indirect_effect_boot_bca_coverage = bca_cover, 
                     indirect_effect_boot_bca_power = bca_power,
                     indirect_effect_boot_perc_lower = fun_med_results$bootstrap_results$indirect_effect_boot_perc_lower,
                     indirect_effect_boot_perc_upper = fun_med_results$bootstrap_results$indirect_effect_boot_perc_upper,
                     indirect_effect_boot_perc_coverage = perc_cover,
                     indirect_effect_boot_perc_power = perc_power
                   ));
  save.image("working.rdata");
}
finish_time <- Sys.time();
print(c(nsim, nboot));
time_required <- difftime(finish_time,start_time);
print(time_required);
mean_answers <- apply(answers,2,mean);
print(round(mean_answers,3));
save.image(paste("answers",nsub,"_",the_seed,".rdata",sep=""));

