################################################################################
################################################################################
###############  Pre - Post treatment extinction validation Runs ###############
################################################################################
################################################################################
setwd("~/Desktop/Repos/resistanceEvolution_Model/src")
source("modelFunc.R")

##* note that N_max is limited to 4.29 billion (4.29 x 10^9) as the upper limit on unsigned int values *##
##* note that N_max of 1.8e19 is a theoretical max using C++ long long, but 1e15 is likely a rational upper limit *##
##* long long values of population sizes are acceptable in the simulation, but dramatically increase computational overhead *##
k_drug = 10 ## fold- increase in death rate due to drug susceptibility
N_max_val = 1e9 ## 1e9 takes ~3.1 min to run 500000 sims
b_fact = 3 ## scalar increase to birth rate
mu = 1e-6 ## per-base-pair mutation rate
N_init <- data.frame(WT = 1) ## initial population size
K_max = 5
N_max = 1e10

## calculate birth / death probabilities
b <- 1.01*b_fact; d <- 1; b_prob <- b/(b+d)
## estimate reasonable upper bounds for T_max 
T_max_val <- ceiling(ceiling(log(base = 2*b_prob, N_max_val)) * 1.2)
T_maxTreat_val <- max(ceiling(ceiling(log(base = 2*b_prob, N_max_val * 2)) * 1.2), 100)

## parameter data frame with default values
pars <- data.frame(
  n_reps = 1000000, n_drugs = 5, T_max = T_max_val, n_retry = 20,
  N_max = N_max_val, T_maxTreat = T_maxTreat_val, N_maxTreat = N_max_val * 2,
  fit_cost = 0.01, frac_del = 0.1,
  R = 1000, k = 0.1,
  resCrit = 1e8, verbose = 1, data_out = 2, keepReps = 50,
  batchname = paste0("fullSimRes"),
  par_dat = paste0("fullSim_pars"),
  exe_path = "./multiType_FullSim.exe",
  saveFiles = F
)


## define parameter grid values
## use value of NA for default ranges [specified in mapVals function]
## use value of NULL to ignore in parameter grid
extinctionMap <- mapVals(list(
  n_drugs = c(1:K_max),
  N_max = exp_range(8,13, decreasing = T),
  mu = exp_range(-3,-10, decreasing = T),
  b_fact = NULL,
  k_drug = NULL,
  #fit_cost = c(.5, exp_range(-1,-3, decreasing = T)[-1], 0)
  fit_cost = c(0, 0.01, 0.05, 0.1)
))


## Set up parallel backend 
ncores <- floor( parallel::detectCores() - 2 )
cl <- makeCluster(ncores) 
registerDoParallel(cl)

## evaluate model over parameter grid
TIC <- Sys.time()
results <- foreach(i = 1:nrow(extinctionMap), .combine = rbind, 
                   .packages = c("dplyr"), .options.snow = list(preschedule = FALSE)) %dopar% {
                     
                     ## extract the i-th parameter set
                     vals <- extinctionMap[i, , drop = FALSE]        # data.frame 1-row
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
                     batch_suffix <- paste0(param_names, "_", sapply(val_list, function(x) paste0(x, collapse = ",")), collapse = "_")
                     pars$batchname <- paste0("fullSimRes_", batch_suffix)
                     pars$par_dat   <- paste0("fullSim_", batch_suffix, "_pars")
                     
                     ## === Handle special parameters that are not plain scalar columns in pars ===
                     ## (e.g., mu might be per-site vector, n_drugs affects dimensions, k_drug used later)
                     # n_drugs override if present in grid
                     if ("n_drugs" %in% param_names) {
                       pars$n_drugs <- as.integer(val_list[["n_drugs"]])
                     }
                     
                     # N_max override if present (used for T_max recomputation, but keep local copy)
                     if ("N_max" %in% param_names) {
                       N_max_val <- as.numeric(val_list[["N_max"]])
                       pars$N_maxTreat = 2*N_max_val
                     }
                     
                     # b_fact or d overrides (used for recomputing b,d and T_max)
                     if (any(c("b_fact","d") %in% param_names)) {
                       # use values either from val_list or fallback to existing variables
                       b_fact_val <- if ("b_fact" %in% param_names) as.numeric(val_list[["b_fact"]]) else b_fact
                       d_val      <- if ("d"      %in% param_names) as.numeric(val_list[["d"]])      else d
                       
                       b <- 3 * b_fact_val
                       d <- d_val
                       b_prob <- b / (b + d)
                       # recompute T_max values using the (possibly updated) N_max_val
                       T_max_val <- ceiling(ceiling(log(N_max_val, base = 2*b_prob)) * 1.2)
                       T_maxTreat_val <- max(ceiling(ceiling(log(N_max_val * 1.2, base = 2*b_prob)) * 1.2), 100)
                       # update pars T_max fields if you want
                       pars$T_max <- T_max_val
                       pars$T_maxTreat <- T_maxTreat_val
                     }
                     
                     ## mu: if present in the grid, capture for site_pars; do NOT assign into pars directly
                     mu_val <- if ("mu" %in% param_names) {
                       as.numeric(val_list[["mu"]])
                     } else {
                       mu    # fallback to global/default mu
                     }
                     
                     ## k_drug: if provided in grid, use it for site_pars; else fallback to existing k_drug variable
                     k_drug_val <- if ("k_drug" %in% param_names) {
                       as.numeric(val_list[["k_drug"]])
                     } else {
                       k_drug
                     }
                     
                     ## define parameters that depend on n_drugs or mu
                     N_B = rep(1, pars$n_drugs) ## **fix so this need not be hard coded**
                     pars$mu_del = get_mu_del(mu = mu_val, frac_del = pars$frac_del) ## **fix so other parameters can be specified easily**
                     
                     ## === Build site_pars (mu repeated per-site, k repeated per-site) ===
                     site_pars <- data.frame(
                       # mu = rep(mu_val, pars$n_drugs), ## old version -- technically equivalent to N_B = 1
                       mu = get_mu_ben(mu_val, N_B), ## updated to define beneficial mutation prob based on per-site values
                       k  = rep(k_drug_val, pars$n_drugs)
                     ) %>% t() %>% as.data.frame()
                     
                     
                     ## define parameters for each class type
                     b_vec <- makeLabels_vec(pars$n_drugs) + b
                     d_vec <- makeLabels_vec(pars$n_drugs) + d
                     type_pars <- make_initial_pop(
                       K = pars$n_drugs, init_size = N_init,
                       b_vec = b_vec, d_vec = d_vec
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
                     num_success <- sum(out$treatSuccess) # process went extinct after introduction of treatment
                     num_resCrit <- sum(out$preResistance > 0 & out$crit > 0) # any resistant genotypes in the population
                     num_FullResCrit <- sum(out$preResistance == 2 & out$crit > 0) # de-novo full mutant
                     num_res <- sum(out$preResistance > 0 ) # any resistant genotypes in the population
                     num_FullRes <- sum(out$preResistance == 2 ) # de-novo full mutant
                     
                     num_preTreatCrit <- sum(out$preTreatCrit) # count of simulations where full mutants at high prevalence before treatment
                     
                     ## return data frame with grid parameters and the associated simulation statistics
                     data.frame(
                       vals,
                       simsUsed = sims_used,
                       fracResolved = frac_resolved,
                       numSuccess = num_success,
                       numResCrit = num_resCrit,
                       numFullResCrit = num_FullResCrit,
                       numRes = num_res,
                       numFullRes = num_FullRes,
                       numPreTreatCrit = num_preTreatCrit,
                       resFilename = paste0("./simulation_files/", pars$batchname, "_results.txt"),
                       repsFilename = paste0("./simulation_files/", pars$batchname, "_full_results.txt")
                     )
                   }
TOC <- Sys.time()
print(TOC - TIC)
stopCluster(cl)

# Merge results back into extinctionMap
extinctionMap <- results; rm(results)
options(scipen = 0)

## Option to read in output of previous simulation from file
# write.table(extinctionMap, file = paste0("./extinctionMap_Nmax6", Sys.Date(), ".txt"), sep = ",", row.names = FALSE, col.names = TRUE)
# extinctionMap = read.table(file.choose(), sep = ',')
# out = read.table(extinctionMap$resFilename[9], sep = ';', header = T); out <- out[!(out$ext == 1), ] 

## --- probability of successful treatment ---
extinctionMap$successFrac <- extinctionMap$numSuccess / extinctionMap$simsUsed

## --- probability of resistance --- (1-P(success))
extinctionMap$resProb <- 1 - extinctionMap$successFrac

## --- Probability of pre-treatment resistance ---
extinctionMap$resFrac <- (extinctionMap$numRes / extinctionMap$simsUsed)

## --- Probability of full pre-treatment resistance ---
### note that this is specifically the number of simulations where a full mutant
### was detected prior to treatment divided by the number of critical simulations.
### Thus values greater than one indicates stochastic extinction, and values less than one
### indicate cases where resistance emerged after starting treatment and became critical
extinctionMap$fullResCritFrac <- (extinctionMap$numFullRes / (extinctionMap$numFullResCrit))

extinctionMap$fracPreResistance <- (extinctionMap$numFullResCrit / extinctionMap$numResCrit)
# safe numeric conversion (handles factors and characters)
extinctionMap <- extinctionMap %>%
  mutate(
    mu_num    = as.numeric(as.character(mu)),
    N_max_num = as.numeric(as.character(N_max))
  ) %>%
  filter(!is.na(mu_num), !is.na(N_max_num))

# unique sorted values for ordering
mu_vals     <- sort(unique(extinctionMap$mu_num), decreasing = TRUE)  # largest -> smallest
Nmax_vals   <- sort(unique(extinctionMap$N_max_num), decreasing = FALSE) # ascending

# formatted strings for parsed labels (scientific notation)
fmt_mu     <- function(x) sprintf("%.0e", x)        # "1e-03" etc
fmt_Nmax   <- function(x) sprintf("%.0e", x)

# construct label strings and factors with correct explicit levels
extinctionMap <- extinctionMap %>%
  mutate(
    mu_label_str   = paste0("mu == ", fmt_mu(mu_num)),
    Nmax_label_str = paste0("N[max] == ", fmt_Nmax(N_max_num)),
    mu_lab = factor(mu_label_str, levels = paste0("mu == ", fmt_mu(mu_vals))),
    N_max_lab = factor(Nmax_label_str, levels = paste0("N[max] == ", fmt_Nmax(Nmax_vals)))
  )

# Map color to the same discrete factor so manual palette lines up
extinctionMap <- extinctionMap %>%
  mutate(mu_fact = factor(mu_num, levels = mu_vals))

# build palette of the right length
n_mu_levels <- length(mu_vals)
pal <- rev(sequential_hcl(n_mu_levels+2, palette = "Heat"))[2:(n_mu_levels+1)]


# -----------------------------
# Facet variables
# -----------------------------
facetVar_list <- c(
  "N_max",
  "mu",
  "fit_cost",
  "b_fact",
  "k_drug"
)

# -----------------------------
# Choose axes
# -----------------------------
x_var <- facetVar_list[2]
y_var <- facetVar_list[1]

# Variables that should use log scale
log_vars <- c("mu", "N_max")

is_log_var <- function(v) v %in% log_vars
log_x <- is_log_var(x_var)
log_y <- is_log_var(y_var)


other_vars <- setdiff(facetVar_list, c(x_var, y_var))
other_vars <- intersect(other_vars, names(extinctionMap))  # only keep those that exist

index_list <- list(
  fit_cost = 1
)
fixed_vals <- lapply(other_vars, function(v) {
  u <- sort(unique(extinctionMap[[v]]))
  
  idx <- index_list[[v]]
  if (is.null(idx)) idx <- 1  # default if not specified
  
  u[idx]
})
names(fixed_vals) <- other_vars


# -----------------------------
# Threshold surface
# -----------------------------
p_crit <- 0.001

threshold <- extinctionMap %>%
  reduce(
    .x = intersect(names(fixed_vals), names(.)),
    .init = .,
    .f = function(df, v) {
      if (v %in% c(x_var, y_var)) return(df)
      df %>% filter(.data[[v]] == fixed_vals[[v]])
    }
  ) %>%
  select(all_of(c("n_drugs", x_var, y_var, other_vars, "resProb", "fracResolved"))) %>%
  filter(fracResolved > 0) %>%
  group_by(across(all_of(c(y_var, x_var, other_vars)))) %>%  
  arrange(n_drugs) %>%
  summarize(
    drugs_needed = {
      vals <- n_drugs[resProb < p_crit]
      if (length(vals) == 0) Inf else min(vals)
    },
    .groups = "drop"
  ) %>%
  filter(is.finite(drugs_needed))

# -----------------------------
# Axis limits
# -----------------------------
x_max <- max(extinctionMap[[x_var]], na.rm = TRUE)
y_max <- max(extinctionMap[[y_var]], na.rm = TRUE)

x_min <- min(threshold[[x_var]])
y_min <- min(threshold[[y_var]])

break_vals <- sort(unique(threshold$drugs_needed))
legend_labels <- as.character(break_vals)

# -----------------------------
# Regression setup
# -----------------------------
x_term <- if (log_x) paste0("log(", x_var, ")") else x_var
y_term <- if (log_y) paste0("log(", y_var, ")") else y_var

gam_formula <- as.formula(
  paste0("drugs_needed ~ s(", x_term, ", ", y_term, ")")
)

gam_fit <- gam(gam_formula, data = threshold)

scam_formula <- as.formula(
  paste0(
    "drugs_needed ~ ",
    "s(", x_term, ", bs='mpi') + ",
    "s(", y_term, ", bs='mpi')"
  )
)

scam_fit <- scam(scam_formula, data = threshold)

# -----------------------------
# Prediction grid
# -----------------------------
x_seq <- seq(min(threshold[[x_var]]), max(threshold[[x_var]]), length.out = 200)
y_seq <- seq(min(threshold[[y_var]]), max(threshold[[y_var]]), length.out = 200)

grid <- expand.grid(x = x_seq, y = y_seq)
names(grid) <- c(x_var, y_var)

grid$gam_pred <- predict(gam_fit, newdata = grid)
grid$scam_pred <- predict(scam_fit, newdata = grid)

# -----------------------------
# Read in "true" pathogen values
# -----------------------------
path_dat <- read.csv("est_vals.csv")

path_dat_wide <- path_dat %>%
  mutate(variable = tolower(variable)) %>%  # ensure consistency
  select(pathogen, variable, value, lwr, upr) %>%
  pivot_wider(
    names_from = variable,
    values_from = c(value, lwr, upr),
    names_glue = "{variable}_{.value}"
  )
path_dat_wide <- path_dat_wide %>%
  rename(
    N_max = n_value,
    N_max_lwr = n_lwr,
    N_max_upr = n_upr,
    mu = mu_value,
    mu_lwr = mu_lwr,
    mu_upr = mu_upr
  )


set.seed(42)


path_dat_bounds <- path_dat_wide %>%
  rowwise() %>%
  mutate(
    N_bounds = list(fill_bounds(N_max, N_max_lwr, N_max_upr, "N_max")),
    mu_bounds = list(fill_bounds(mu, mu_lwr, mu_upr, "mu"))
  ) %>%
  mutate(
    N_max_lwr = N_bounds$lwr,
    N_max_upr = N_bounds$upr,
    mu_lwr = mu_bounds$lwr,
    mu_upr = mu_bounds$upr
  ) %>%
  ungroup() %>%
  select(-N_bounds, -mu_bounds)

path_dat_ellipse <- path_dat_bounds %>%
  mutate(id = row_number()) %>%
  rowwise() %>%
  mutate(
    ellipse = list(
      make_ellipse(
        cx = .data[[x_var]],
        cy = .data[[y_var]],
        lx = .data[[paste0(x_var, "_lwr")]],
        ux = .data[[paste0(x_var, "_upr")]],
        ly = .data[[paste0(y_var, "_lwr")]],
        uy = .data[[paste0(y_var, "_upr")]],
        log_x = log_x,
        log_y = log_y
      )
    )
  ) %>%
  unnest(ellipse)

# -----------------------------
# Generate hypothetical dataset
# -----------------------------
n_points <- 1
set.seed(42)

df_bounds <- tibble(
  name = paste0("Pathogen_", seq_len(n_points)),
  !!x_var := sample_param(x_var, n_points),
  !!y_var := sample_param(y_var, n_points)
)

# -----------------------------
# Construct uncertainty ranges
# -----------------------------
df_bounds <- df_bounds %>%
  rowwise() %>%
  mutate(
    
    # span of ellipse in x
    x_span = if (log_x) {
      # span orders of magnitude
      runif(1, 0.6, 1.3)  # orders of magnitude
    } else {
      # proportional width
      .data[[x_var]] * runif(1, 0.05, 0.2)
    },
    
    # span of ellipse in y
    y_span = if (log_y) {
      runif(1, 0.6, 1.3)
    } else {
      .data[[y_var]] * runif(1, 0.45, 0.70)
    }
  ) %>%
  mutate(
    
    # Convert spans into bounds
    x_lwr = if (log_x) {
      10^(log10(.data[[x_var]]) - x_span/2)
    } else {
      .data[[x_var]] - x_span
    },
    
    x_upr = if (log_x) {
      10^(log10(.data[[x_var]]) + x_span/2)
    } else {
      .data[[x_var]] + x_span
    },
    
    y_lwr = if (log_y) {
      10^(log10(.data[[y_var]]) - y_span/2)
    } else {
      .data[[y_var]] - y_span
    },
    
    y_upr = if (log_y) {
      10^(log10(.data[[y_var]]) + y_span/2)
    } else {
      .data[[y_var]] + y_span
    }
    
  ) %>%
  ungroup()

# -----------------------------
# Build ellipse dataframe
# -----------------------------
ellipse_df <- df_bounds %>%
  mutate(id = row_number()) %>%
  rowwise() %>%
  mutate(
    ellipse = list(
      make_ellipse(
        cx = .data[[x_var]],
        cy = .data[[y_var]],
        lx = x_lwr,
        ux = x_upr,
        ly = y_lwr,
        uy = y_upr,
        log_x = log_x,
        log_y = log_y
      )
    )
  ) %>%
  unnest(ellipse) %>%
  ungroup()


# -----------------------------
# Generate plots
# -----------------------------
p <- plot_surface("gam_pred", "GAM smooth")
q <- plot_sim_surface()
r <- plot_surface("scam_pred", "Monotone SCAM")

p
q
r

combined_plot <- q + r + 
  plot_layout(ncol = 3) &   # 3 columns
  plot_annotation(
    title = "Comparison of Drug Resistance Thresholds",
    subtitle = paste0("Axes: ", x_var, " vs ", y_var),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

combined_plot

combined_plot <- q + p +
  plot_layout(ncol = 3) &   # 3 columns
  plot_annotation(
    title = "Comparison of Drug Resistance Thresholds",
    subtitle = paste0("Axes: ", x_var, " vs ", y_var),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

combined_plot

combined_plot <- p + q + r +
  plot_layout(ncol = 3) &   # 3 columns
  plot_annotation(
    title = "Comparison of Drug Resistance Thresholds",
    subtitle = paste0("Axes: ", x_var, " vs ", y_var),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

combined_plot

x_var = "n_drugs"
row_var = facetVar_list[2]
col_var = facetVar_list[1]
color_var = facetVar_list[3]


n_Nmax <- 6   # number you want to keep
n_mu   <- 6

N_vals <- sort(unique(extinctionMap[,row_var]))
mu_vals <- sort(unique(extinctionMap[,col_var]))

N_keep <- N_vals[round(seq(1, length(N_vals), length.out = min(n_Nmax, length(N_vals))))]
mu_keep <- mu_vals[round(seq(1, length(mu_vals), length.out = min(n_mu, length(mu_vals))))]

extinctionMap_trim <- extinctionMap %>%
  filter(.data[[row_var]] %in% N_keep, .data[[col_var]] %in% mu_keep) 

p1 <- plot_metric_grid(extinctionMap_trim, 
                       y_var = "resProb", 
                       x_var = x_var,
                       facet_row = row_var, 
                       facet_col = col_var, 
                       const_vals = list(b_fact = b_fact, k_drug = k_drug),
                       color_var = color_var,
                       transform = "log")
p1
