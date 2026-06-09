##########################
##########################
##   Resistance model   ##
## Sensitivity analysis ##
##########################
##########################
setwd("~/Desktop/Repos/resistanceEvolution_Model/src")
source("modelFunc.R")

## default default main parameters
R0 = 5 # expected number of offspring / survival time
b = 1 - (1 / (1+R0)) # individual per-generation birth probability
  R0_range = c(1,10)
  lh_pars_range = data.frame(R0 = R0_range)

mu = 1e-5 # per-base-per-replication mutation rate
  mu_range = c(1e-10, 1e-3)
  lh_pars_range = cbind(lh_pars_range, data.frame(mu = mu_range))
  
N_max = 1e12 # max population at start of treatment
  N_max_range = c(1e8, 1e13)
  lh_pars_range = cbind(lh_pars_range, data.frame(N_max = N_max_range))
  
w = 0.9 # omega - drug susceptible fitness cost - can be a vector of length K
  w_range = c(0.9, 1)
  lh_pars_range = cbind(lh_pars_range, data.frame(w = w_range))
  
  fit_cost = 0.1 # per-mutation fitness cost
  fit_cost_range = c(0,0.2)
  lh_pars_range = cbind(lh_pars_range, data.frame(fit_cost = fit_cost_range))
  
n_drugs = 5 # number of drugs used -- start at 0 for range
  n_drugs_range = c(0,5)
  lh_pars_range = cbind(lh_pars_range, data.frame(n_drugs = n_drugs_range))

G = 1e6 # genome size
  G_range = c(1e3, 1e8)
  lh_pars_range = cbind(lh_pars_range, data.frame(G = G_range))

frac_del = 0.1 # fraction of lethal deleterious sites
Z = G*frac_del

rownames(lh_pars_range) <- c("lwr", "upr")

## initial population size
N_init <- data.frame(WT = 1) ## initial population size

## estimate default upper bound on generation counts to prevent near-infinite simulations
T_max_val <- ceiling(ceiling(log(base = 2*b, N_max)) * 1.2)
T_maxTreat_val <- max(ceiling(ceiling(log(base = 2*b, N_max * 2)) * 1.2), 100)

## all parameters with defaults
pars <- data.frame(
  n_drugs = n_drugs, 
  T_max = T_max_val, T_maxTreat = T_maxTreat_val, 
  N_max = N_max, N_maxTreat = N_max * 2,
  b = b, fit_cost = fit_cost, mu = mu, frac_del = frac_del, Z = Z,
  
  n_reps = 10, n_retry = 10, resCrit = 1e8, 
  verbose = 0, data_out = 1, keepReps = 20, ## data_out = 3 gives minimal output for sensitivity analysis
  batchname = paste0("fullSimRes"),
  par_dat = paste0("fullSim_pars"),
  exe_path = "./multiType_FullSim_p.exe",
  saveFiles = F
)

## parameters to use in sensitivity 
#### only built to work over parameters 
#### specified in get_lhs_pars function
sensitivity_vars = c("R0", "mu", "N_max", "w", "fit_cost", "n_drugs", "G")
log_vars = c("mu", "N_max", "G")


lh_pars_valid = get_lh_pars(N_LH_sets = 1000, lh_pars = lh_pars_range, par_names = sensitivity_vars, sampling_scheme = "unif", vary_single = 100)
lh.valid.summary <- as.data.frame(round(apply(lh_pars_valid[,c(sensitivity_vars)], 2, function(x) summary(x)),4)) #parameter summary file


## run with lhs parameters -- valid sets only
## Set up parallel backend 
ncores <- floor( parallel::detectCores() - 2 )
cl <- makeCluster(ncores) 
registerDoParallel(cl)

## evaluate model over parameter grid
TIC <- Sys.time()
sens_results <- foreach(i = 1:nrow(lh_pars_valid), .combine = rbind, 
                   .packages = c("dplyr"), .options.snow = list(preschedule = FALSE)) %dopar% {
                     
                     ## extract the i-th parameter set
                     vals <- lh_pars_valid[i, , drop = FALSE]        # data.frame 1-row
                     param_names <- names(vals)
                     val_list <- as.list(vals[1, ])                  # named scalars for easy access
                     
                     ## === Assign only those params that exist in pars ===
                     pars_cols_to_set <- intersect(param_names, names(pars))
                     if (length(pars_cols_to_set) > 0) {
                       # create a 1-row data.frame with the same column order as pars_cols_to_set
                       tmp <- vals[ , pars_cols_to_set, drop = FALSE]
                       # assign into pars (recycling/scalar assignment is fine since tmp is 1-row)
                       pars[ , pars_cols_to_set] <- tmp[rep(1, nrow(pars)), , drop = FALSE]
                     }
                     
                     ## update batchname and par_dat using the grid values (stringify safely)
                     batch_suffix <- paste0("sensRun_", i)
                     pars$batchname <- paste0("fullSimRes_", batch_suffix)
                     pars$par_dat   <- paste0("fullSim_", batch_suffix, "_pars")
                     
                     
                     
                     ## mu: if present in the grid, capture for site_pars; do NOT assign into pars directly
                     mu_val <- if ("mu" %in% param_names) {
                       as.numeric(val_list[["mu"]])
                     } else {
                       mu    # fallback to global/default mu
                     }
                     
                     
                     ## k_drug: if provided in grid, use it for site_pars; else fallback to existing k_drug variable
                     w_val <- if ("w" %in% param_names) {
                       as.numeric(val_list[["w"]])
                     } else {
                       w # fallback to global default
                     }
                     

                     ## define parameters that depend on n_drugs or mu
                     #pars$mu_del = get_mu_del(mu = mu_val, frac_del = pars$frac_del) ## **fix so other parameters can be specified easily**
                     pars$mu_del <- if("Z" %in% param_names){
                       1 - (1-mu_val)^as.numeric(val_list[["Z"]])
                     } else {
                       1 - (1-mu_val)^Z # fallback to global default
                     }

                     ## number of sites that give resistance to a given drug -- here, assume only 1 site
                     N_B = rep(1, pars$n_drugs) ## **fix so this need not be hard coded**
                     
                     ## === Build site_pars (mu repeated per-site, k repeated per-site) ===
                     site_pars <- data.frame(
                       mu = get_mu_ben(mu_val, N_B), ## updated to define beneficial mutation prob based on per-site values
                       w = rep(w_val, pars$n_drugs)
                     ) %>% t() %>% as.data.frame()
                     
                     
                     ## define parameters for each class type
                     type_pars <- make_initial_pop_p(
                       K = pars$n_drugs, init_size = N_init
                     )
                     
                     
                     ## Run the model
                     tic <- Sys.time()
                     out <- multiType_FullSim(
                       n_reps = pars$n_reps, 
                       type_pars = type_pars, site_pars = site_pars, 
                       parameters = pars, 
                       exe_path = pars$exe_path,
                       run_dir = NULL, 
                       saveFiles = pars$saveFiles
                     )
                     toc <- Sys.time()
                     tictoc <- toc - tic
                     
                     
                     ## remove replicates where extinction occurred before treatment
                     out <- out[!(out$ext == 1), ]
                     
                     
                     ## calculate simulation statistics
                     sims_used <- nrow(out)
                     frac_resolved  <- (sum(out$crit) + sum(out$treatSuccess)) / nrow(out)
                     frac_resistance <- (nrow(out) - sum(out$treatSuccess)) / nrow(out) # fraction of simulations where resistance emerged
                     
                     ## return data frame with grid parameters and the associated simulation statistics
                     data.frame(
                       vals,
                       simsUsed = sims_used,
                       fracResolved = frac_resolved,
                       fracResistance = frac_resistance,
                       resFilename = paste0("./simulation_files/", pars$batchname, "_results.txt"),
                       repsFilename = paste0("./simulation_files/", pars$batchname, "_full_results.txt")
                     )
                   }
TOC <- Sys.time()
print(TOC - TIC)
stopCluster(cl)


# remove any simulations where none resolved (i.e., extinction conditions)
sum(sens_results$simsUsed == 0)
sens_results <- sens_results %>%
  filter(sens_results$simsUsed != 0)

# split out baseline parameters
baseline_vals <- sens_results %>%
  filter(base_par %in% "baseline") %>%
  select(all_of(c(sensitivity_vars, "base_id", "base_par", "frac_del", "fracResistance", "simsUsed")))

sens_results_trim <- sens_results %>%
  filter(base_par != "baseline") %>%
  select(all_of(c(sensitivity_vars, "base_id", "base_par", "frac_del", "fracResistance", "simsUsed")))



## loop over each parameter
mode = 5; model = "log"; minDat = 5; sigChange = 0.5
results <- imap(
  setNames(unique(sens_results_trim$base_par),
           unique(sens_results_trim$base_par)),
  ~ process_param(.y, sens_results_trim,
                  mode = mode, model = model, 
                  minDat = minDat, sigChange = sigChange)
)


## combine results in singel dataframe
slope_df <- bind_rows( map(results, "reg"), .id = "parameter" )
regDat <- bind_rows( map(results, "data"), .id = "parameter" )
rm(results)

names(slope_df) <- c("base_par", "base_id", "slope", "diff", "diff_sign", "val_min", "val_max", "mean", "ps", "fs")

slope_df$logRange = if(model == "log"){ slope_df$val_max - slope_df$val_min }else{ log10(slope_df$val_max) - log10(slope_df$val_min) }
slope_df_trim <- slope_df %>%
  filter(!is.infinite(slope)) %>%
  filter(!is.na(slope))

sens_results_filtered <- sens_results_trim %>%
  semi_join(slope_df_trim, by = c("base_id", "base_par"))

## convert all parameter values to log scale
slope_df_trim <- slope_df_trim %>%
  mutate(
    val_min = if_else( !(base_par %in% log_vars), 
                       log10(val_min), val_min ),
    
    val_max = if_else( !(base_par %in% log_vars), 
                       log10(val_max), val_max )
  )


## negative correlations in mu
badID = slope_df$base_id[which(slope_df$base_par == "mu" & slope_df$diff_sign < -.5)]
max(sens_results_trim$n_drugs[sens_results_trim$base_id %in% badID & sens_results_trim$base_par %in% "mu"])
length(badID)
length(badID)/sum(slope_df$base_par == "mu") ## fraction of parameter sets with non-monotonic or degenerate edge behavior
id = badID[1];
if(length(badID) > 0 & F){
  plot_monotonicity( 
    sens_results = sens_results_trim[(sens_results_trim$base_id == id & sens_results_trim$base_par == "mu") ,], 
    sens_vars = "mu", # sensitivity_vars, 
    output_var = "fracResistance",
    save_plots = F, out_dir = "./figures/", tag = "_unifGrid",
    width = 840, height = 840, res = 100
  )
}
bad.rm = T
if(bad.rm){ 
  slope_df_trim <- slope_df_trim %>% 
    filter(!(base_par == "mu" & base_id %in% badID))
}
  



## sensitivity/importance metric -- how often there was a shift from 0-1
alpha_min = 0.1
alpha = 0.05
slope_df_trim %>% 
  group_by(base_par) %>%
  summarize(
    detectChange = sum(diff > alpha_min),
    detectChangePct = sum(diff > alpha_min) / length(diff),
    sigChange = sum(diff > 1-alpha),
    sigChangePct = sum(diff > 1-alpha) / length(diff) ) %>%
  ungroup()

## log-unit parameter ranges where change was detected
slope_df_trim %>% 
  group_by(base_par) %>%
  summarize(
    rangeSize = mean(logRange, na.rm = T) ) %>%
  ungroup()



###########
##  PRCC ##
###########
library(epiR)

fracRes <- regDat[ c(sensitivity_vars, 'fracResistance') ]
prcc.fracRes <- epi.prcc(fracRes, sided.test = 2)


p_vals <- data.frame(parameters = sensitivity_vars)
p_vals$fracRes_pval <- prcc.fracRes$p.value

prcc.fracRes <- prcc.fracRes[order(prcc.fracRes$est),] ## order variables by PRCC value
p_vals <- p_vals[order(match(p_vals$parameters, prcc.fracRes$var)),] ## order p values based on PRCC ordering
prcc.fracRes <- prcc.fracRes %>%
  select(var, est) %>%
  rename(parameters = var, fracRes.est = est
  )

prcc.trans = prcc.fracRes

## only necessary if more than one  metric is  being observed (eg. TruPrev, Fadeout, etc. in this example)
# prcc.trans <- cbind(TruPrev = prcc.TruPrev, Fadeout = prcc.Fadeout$est, FadeoutTime = prcc.FadeoutTime$est, HuntPrev = prcc.HuntPrev$est)
# prcc.trans <- prcc.trans[order(prcc.trans$TruPrev.est),]
# p_vals <- p_vals[order(match(p_vals$parameters, prcc.TruPrev$parameters)),]

#### select metric columns and convert to long format
# prcc.trans <- prcc.trans[,c("resFrac.est", "parameters", "Fadeout.est", "FadeoutTime.est", "HuntPrev.est")]

p_vals.long <- p_vals %>%
  pivot_longer( cols = -parameters,
                names_to = "metric", values_to = "value" )


prcc.trans.long <- prcc.trans %>%
  pivot_longer( cols = -parameters, 
                names_to = "metric", values_to = "value" )


#Renaming facets -- currently only 1
variable_names <- list(
  # "Fadeout" = "Fadeout",
  # "FadeoutTime" = "Fadeout Time" ,
  # "HuntPrev" = "Observed Prevalence",
  "fracRes.est" = "Probability of resistance"
)

## variable_labeller <- function(variable,value){ return(variable_names[value]) }
variable_labeller <- as_labeller(variable_names) ## ai suggestion
significant = F
if(significant){
  prcc.trans.long$value[p_vals.long$value >= .05] <- NA
}


blue<-c('#2171b5')
#plot PRCC results

prcc.trans.long <- prcc.trans.long %>%
  mutate(parameters = factor(parameters, levels = unique(parameters)))



param_labels <- c(
  mu       = "atop('Mutation', 'rate ('*mu*')')",
  G        = "atop('Genome', 'size ('*italic(G)*')')",
  w        = "atop('Drug killing', 'rate ('*kappa*')')",
  R0       = "atop('Reproductive',  'number ('*R[0]*')')",
  N_max    = "atop('Population', 'size ('*N[max]*')')",
  fit_cost = "atop('Fitness', 'cost (c)')",
  n_drugs  = "'# Drugs'"
)


sens <- ggplot() + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.y = element_text(size = 14),  axis.title.x = element_text(size = 14)) +
  geom_col(data = prcc.trans.long, aes(y = value, x = parameters), fill = blue) +
  #scale_x_discrete(position = "top") +
  #scale_y_discrete(limits=c(0, 0.25, 0.5), expand=c(0.02,0.02)) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Parameter", y = "Rank correlation coefficient") + 
  scale_x_discrete( labels = function(x) parse(text = param_labels[x]) ) + 
  
  coord_flip() + 
  facet_grid(~metric, labeller = variable_labeller) +
  scale_fill_manual(values = blue) + 
  theme(panel.spacing.y = unit(1.5, "lines")) + 
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16,),
    
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
  )
#sens <- annotate_figure(sens, top = text_grob(paste0("PRCC Analysis  DS.mode = ", mode), color = "black", face = "bold", size = 18))
sens

# jpeg(filename = paste0("PRCC_Trans_fracRes.jpeg"), width = 1440, height = 840, units = 'px', res = 100)
# print(sens)
# dev.off()

###########

#########################
## sensitivity metrics ##
#########################
# ggplot( data = slope_df_trim, aes(x = base_par, y = log10(abs(slope)), color = base_par) ) +
#   geom_boxplot()
slope_df_trim$effectSpan = slope_df_trim$slope*slope_df_trim$logRange
effSpan_box <- ggplot( data = slope_df_trim %>% filter(!(base_par == "mu" & slope < 0)), 
                       aes(x = base_par, y = effectSpan, color = base_par) ) +
  geom_boxplot() +
  labs( x = "Parameter", y = "effectSpan Metric", color = "Parameter")
effSpan_box <- annotate_figure(effSpan_box, top = text_grob(paste0("Regression Analysis  DS.mode = ", mode), 
                                                            color = "black", face = "bold", size = 18))
effSpan_box

## 
effect_summary <- slope_df_trim %>%
  group_by(base_par) %>%
  summarise(
    mean = mean(effectSpan, na.rm = TRUE),
    median = median(effectSpan, na.rm = TRUE),
    sd = sd(effectSpan, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(parameters = base_par) %>%
  mutate(metric = "fracResistance", cv = sd/mean) %>%
  mutate(parameters = factor(parameters, levels = prcc.trans.long$parameters))
x=effect_summary$mean / effect_summary$median; names(x)=effect_summary$parameters; x


blue = '#2171b5'
effSpan <- ggplot(data = effect_summary, aes(y = mean, x = parameters, fill = blue)) + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.y = element_text(size = 14),  axis.title.x = element_text(size = 14)) +
  geom_col() + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),width = 0.2) + 
  #scale_y_continuous(limits=c(-1, 1)) +
  coord_flip() + 
  facet_grid(~metric) +
  scale_fill_manual(values = blue) + 
  labs( x = "Parameter", y = "Standardized slope" ) +
  scale_x_discrete( labels = function(x) parse(text = param_labels[x]) ) + 
  theme(panel.spacing.y = unit(1.5, "lines")) + 
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16,),
    
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
  )
# effSpan <- annotate_figure(effSpan, top = text_grob(paste0("Regression Analysis  DS.mode = ", mode), 
#                                                     color = "black", face = "bold", size = 18))
effSpan



#########################




















# ## distribution of pct. change in resistance by parameter
# ggplot(data = slope_df_trim) +
#   geom_histogram( aes(x = diff_sign, y = after_stat(density)), 
#                   binwidth = 0.05, boundary = -1 ) +
#   facet_wrap( ~base_par, scales = "free" ) +
#   theme_bw( base_size = 12 ) +
#   xlim(c(-1,1)) +
#   labs( x = "Max Change (all)", y = "Density" ) +
#   theme( legend.position = "none" )
# 
# ## remove all-or-nothing cases
# ggplot(data = slope_df_trim[!(slope_df_trim$diff %in% c(0,1)),]) +
#   geom_histogram( aes(x = diff_sign, y = after_stat(density)), 
#                   binwidth = 0.05, boundary = -1 ) +
#   facet_wrap( ~base_par, scales = "free" ) +
#   theme_bw( base_size = 12 ) +
#   labs( x = "Max Change (intermediate)", y = "Density" ) +
#   theme( legend.position = "none" )


####################################
## magnitude of range over which ##
##   detectable change occurs    ##
####################################

fullRange <- sens_results[,sensitivity_vars] %>%
  pivot_longer(
    cols = everything(),
    names_to = "base_par",
    values_to = "value"
  ) %>% 
  group_by(base_par) %>%
  summarise(
    logRange_Full = log10( max(value, na.rm = TRUE) ) - log10( min(value, na.rm = TRUE) ),
    .groups = "drop"
  )

slope_df_trim <- slope_df_trim %>%
  left_join(fullRange, by = "base_par") %>%
  mutate(range_frac = logRange / logRange_Full)


## raw range values
ggplot(data = slope_df_trim) +
  geom_histogram(aes(x = logRange, y = ..density..)) +
  facet_wrap( ~base_par, scales = "free" ) +
  theme_bw( base_size = 12 ) +
  labs( x = "log-range size", y = "Density" ) +
  theme( legend.position = "none" )
## normalized range values by max possible range
ggplot(data = slope_df_trim) +
  geom_histogram(aes(x = range_frac, y = ..density..)) +
  geom_density( aes(x = range_frac), linewidth = 1, color = "blue") +
  facet_wrap( ~base_par, scales = "free" ) +
  theme_bw( base_size = 12 ) +
  labs( x = "log-range size", y = "Density" ) +
  theme( legend.position = "none" )

####################################



## convert relevant params to  log scale if not already
slope_df_trim <- slope_df_trim %>% mutate(
  val_min = if_else((base_par %in% c("mu", "N_max", "G") & model != "log"), log10(val_min), val_min),
  val_max = if_else((base_par %in% c("mu", "N_max", "G") & model != "log"), log10(val_max), val_max),
  mean = if_else((base_par %in% c("mu", "N_max", "G") & model != "log"), log10(mean), mean)
)



plot_df <- slope_df_trim %>%
  select(base_par, val_min, val_max) %>%
  pivot_longer(cols = c(val_min, val_max), names_to = "bound_type", values_to = "value")
ggplot(plot_df, aes(x = value, fill = bound_type)) +
  geom_histogram(bins = 30, alpha = 0.5, position = "identity") +
  facet_wrap(~base_par, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Parameter value",fill = "Bound")
rm(plot_df)
slope_df_trim <- slope_df_trim %>% mutate(
  val_min = if_else((base_par %in% c("mu", "N_max", "G") & model != "log"), 10^(val_min), val_min),
  val_max = if_else((base_par %in% c("mu", "N_max", "G") & model != "log"), 10^(val_max), val_max)
)




## view particular regression ranges
var = "w"; id = 604
sens_results_trim$fracResistance.logit = logit(sens_results_trim$fracResistance*0.99 + 0.005)
plot_monotonicity( 
  sens_results = sens_results_trim[(sens_results_trim$base_id == id & sens_results_trim$base_par == var) ,], 
  sens_vars = var, # sensitivity_vars, 
  output_var = "fracResistance.logit",
  save_plots = F, out_dir = "./figures/", tag = "_unifGrid",
  width = 840, height = 840, res = 100
)

hist_df <- sens_results_trim %>%
  filter( base_par %in% sensitivity_vars ) %>%
  pivot_longer( cols = all_of(sensitivity_vars),
                names_to = "variable", values_to = "value_raw" ) %>%
  filter(variable == base_par) %>%
  mutate(
    value = if_else( variable %in% log_vars, log10(value_raw), value_raw),
    scale_type = if_else( variable %in% log_vars,"log10", "linear")
  )


ggplot(hist_df, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  facet_wrap(~variable, scales = "free") +
  theme_bw(base_size = 12) +
  labs(x = "Parameter value", y = "Density")
rm(hist_df)


## -----------------------------
## Monotonicity plot generator
## -----------------------------
library(ggplot2)
library(gridExtra)


# plot_monotonicity( 
#   sens_results = regDat, 
#   sens_vars = "mu", # sensitivity_vars, 
#   output_var = "fracResistance.logit", #"fracResistance"
#   save_plots = F, out_dir = "./figures/", tag = "_unifGrid",
#   width = 840, height = 840, res = 100
# )







####################
## Logistic Model ##
####################
library(sjPlot)
library(sjmisc)
library(sjlabelled)


sensitivity_vars_trans = paste0("log10(", sensitivity_vars, ")")
#single.model = as.formula( paste( 'RESPONSE~', paste(sensitivity_vars, collapse = '+') ) )
single.model = as.formula( paste( 'RESPONSE~', paste(sensitivity_vars_trans, collapse = '+') ) )

p_vals_lm <- data.frame(parameters = sensitivity_vars)
#TruPrev: No interaction regression 
single.model.fracRes = update.formula(single.model, `fracResistance.trans` ~ .)

regDat$fracResistance.trans = logit(regDat$fracResistance*0.99 + 0.005)
lm.fracRes = lm(single.model.fracRes, data = regDat)
#write.csv(lm.fracRes.PRCC, file=paste0("lmOutput_fracRes.csv"), row.names = F)
summary(lm.fracRes)
tab_model(lm.fracRes)

values.in.PRCC <- lm.fracRes$coefficients
values.in.PRCC <- as.data.frame(values.in.PRCC)
values.in.PRCC$par <- rownames(values.in.PRCC)
values.in.PRCC <- values.in.PRCC[-1,]
values.in.PRCC.scaled <- values.in.PRCC$values.in.PRCC/max(abs(values.in.PRCC$values.in.PRCC)) #regression results scaled by the largest abs regression value

lm.fracRes.scaled <- data.frame(values.in.PRCC.scaled, row.names = sensitivity_vars)

p_vals_lm$fracRes <- summary(lm.fracRes)$coefficients[-1,4]

lm.fracRes.scaled$parameters <- rownames(lm.fracRes.scaled)
lm.fracRes.scaled <- lm.fracRes.scaled[order(match(lm.fracRes.scaled$parameters, prcc.fracRes$parameters)),]
p_vals_lm <- p_vals_lm[order(match(p_vals_lm$parameters, prcc.fracRes$parameters)),]

p_vals_lm.long <- p_vals_lm %>%
  pivot_longer(cols = -parameters, 
               names_to = "variable", values_to = "value")

lm.fracRes.scaled.long <- lm.fracRes.scaled %>%
  pivot_longer( cols = -parameters, 
                names_to = "metric", values_to = "value" )
lm.fracRes.scaled.long$metric = "fracRes.est" ## only works for single metric plots

if(significant){
  lm.fracRes.scaled.long$value[p_vals_lm.long$value >= .05] <- NA
}

#plot single model results
lm.fracRes.scaled.long <- lm.fracRes.scaled.long %>%
  mutate(parameters = as.factor(parameters))
reg <- ggplot() + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.y = element_text(size = 14),  axis.title.x = element_text(size = 14)) +
  geom_col(data = lm.fracRes.scaled.long, aes(y = value, x = parameters, fill = blue)) + 
  scale_y_continuous(limits=c(-1, 1)) +
  coord_flip() + 
  facet_grid(~metric, labeller = variable_labeller) +
  scale_fill_manual (values = blue) + 
  theme(panel.spacing.y = unit(1.5, "lines")) + 
  theme(strip.text = element_text(size = 12))
reg <- annotate_figure(reg, top = text_grob("Regression Analysis", 
                                            color = "black", face = "bold", size = 18))
reg

jpeg(filename = paste0("SingleModel_Trans_", infType, ".jpeg"), width = 1440, height = 840, units = 'px', res = 100)
reg
dev.off()
