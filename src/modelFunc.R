####################
## Load necessary ##
##    packages    ##
####################
# install.packages("ggplot2")
# install.packages("viridis")
# install.packages("tidyverse")
# install.packages("patchwork")
# install.packages("ggpubr")
# install.packages("colorspace")
# install.packages("data.table")
# install.packages("doParallel")
# install.packages("foreach")

library(ggplot2)
library(viridis)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(GGally)
library(colorspace)
library(data.table)
library(doParallel)
library(foreach)
library(readxl)
library(mgcv)
library(scam)
library(gridExtra)

library(dplyr)
library(purrr)
library(tidyr)


multiType_FullSim <- function(n_reps, type_pars, site_pars, parameters, 
                              exe_path = "./multiType_FullSim.exe", run_dir = NULL, 
                              saveFiles = TRUE) {
  
  parvec <- c(format(n_reps, scientific = FALSE),
              parameters$n_drugs, # number of independent drugs/mutation sites
              parameters$T_max, # max generation number
              format(parameters$N_max, scientific = FALSE), # critical population size
              parameters$T_maxTreat, # max generation number during treatment
              format(parameters$N_maxTreat, scientific = FALSE), # critical population size during treatment
              format(parameters$resCrit, scientific = FALSE), # critical resistant pop size
              parameters$verbose, # debugging error output
              parameters$data_out, # print to console: 0; otherwise print to file
              parameters$keepReps, # number of full replicate trajectories to save to file
              parameters$batchname, # name element for output / error files
              parameters$par_dat, # specified name for parameter file
              parameters$fit_cost,
              parameters$n_retry,
              parameters$mu_del, 
              parameters$b
  ) 
  strvec <- format(parvec, digits = 5)
  
  
  # Optionally change working directory for the run
  oldwd <- getwd()
  if (!is.null(run_dir)) setwd(run_dir)
  
  options(scipen = 999)
  ## generate parameter file for vectorized parameters
  write.table(type_pars, file = paste0("./simulation_files/", parameters$par_dat, ".txt"), sep = ",", row.names = FALSE, col.names = FALSE)
  write.table(site_pars, file = paste0("./simulation_files/", parameters$par_dat, ".txt"), sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  
  ## print command line statement for comparison
  if(parameters$verbose > 1){ print( paste(shQuote(exe_path), paste(strvec, collapse = " ")) ) }
  
  res <- system2(exe_path, args = strvec, stdout = FALSE)
  
  out <- read.table(paste0("./simulation_files/", parameters$batchname, "_results.txt"), header = TRUE, sep = ";", check.names = FALSE) %>%
    mutate_all(as.numeric)
  
  ## remove all generated files
  
  if(!saveFiles){
    files = c(paste0("./simulation_files/", parameters$batchname, "_results.txt"), ## results file
              paste0("./simulation_files/", parameters$batchname, "_full_results.txt"), ## replicate trajectory file
              paste0("./simulation_files/", parameters$batchname, "_err.txt"), ## error file
              paste0("./simulation_files/", parameters$par_dat, ".txt") ## parameter file
    )
    file.remove(files)
  }
  
  
  options(scipen = 0)
  if (!is.null(run_dir)) setwd(oldwd)
  return( out )
}

# ---------------------------
# Function to generate initial state vector, 
## returned using bitmap ordering
## for consistency with simulation code
# K defines the number of drugs (so # classes = 2^K)
# init_size is a data frame for known initial states {WT, M1, ... MXX}
## values not included are assumed to be 0
# ---------------------------
make_initial_pop <- function(K, init_size = data.frame(), b_vec, d_vec) {
  # number of genotypes
  n_geno <- 2^K
  
  # generate bitmask labels
  bitmask_labels <- sapply(0:(n_geno - 1), function(mask) {
    if (mask == 0) {
      "WT"
    } else {
      mutated_sites <- which(intToBits(mask)[1:K] == 1)
      paste0("M", paste(mutated_sites, collapse = ""))
    }
  })
  
  # initialize all counts as zero
  counts <- setNames(rep(0, n_geno), bitmask_labels)
  
  # check that init_size has correct column names
  if (!missing(init_size) && ncol(init_size) > 0) {
    provided_names <- colnames(init_size)
    bad_names <- setdiff(provided_names, bitmask_labels)
    
    # warn if names don't match
    if (length(bad_names) > 0) {
      warning("Ignoring invalid initial states: ", paste(bad_names, collapse = ", "))
    }
    
    # assign matching values
    good_names <- intersect(provided_names, bitmask_labels)
    counts[good_names] <- init_size[1, good_names]
  }
  
  # return as a one-row data frame
  # reorder b_vec and d_vec to match counts
  b_vec_ord <- b_vec[names(counts)]
  d_vec_ord <- d_vec[names(counts)]
  
  # construct final 3-row data frame
  df <- rbind(
    N = counts,
    b = b_vec_ord,
    d = d_vec_ord
  )
  
  return(df)
}

make_initial_pop_p <- function(K, init_size = data.frame()) {
  # number of genotypes
  n_geno <- 2^K
  
  # generate bitmask labels
  bitmask_labels <- sapply(0:(n_geno - 1), function(mask) {
    if (mask == 0) {
      "WT"
    } else {
      mutated_sites <- which(intToBits(mask)[1:K] == 1)
      paste0("M", paste(mutated_sites, collapse = ""))
    }
  })
  
  # initialize all counts as zero
  counts <- setNames(rep(0, n_geno), bitmask_labels)
  
  # check that init_size has correct column names
  if (!missing(init_size) && ncol(init_size) > 0) {
    provided_names <- colnames(init_size)
    bad_names <- setdiff(provided_names, bitmask_labels)
    
    # warn if names don't match
    if (length(bad_names) > 0) {
      warning("Ignoring invalid initial states: ", paste(bad_names, collapse = ", "))
    }
    
    # assign matching values
    good_names <- intersect(provided_names, bitmask_labels)
    counts[good_names] <- init_size[1, good_names]
  }
  
  # construct final 3-row data frame
  df <- rbind(
    N = counts
  )
  
  return(df)
}

## function to generate a default named vector 
#### for use in make_initial_pop (b_vec and d_vec)
makeLabels_vec <- function(K, bitwise = FALSE) {
  # number of genotypes
  n_geno <- 2^K
  
  # generate all masks
  masks <- 0:(n_geno - 1)
  
  # turn masks into labels
  labels <- sapply(masks, function(mask) {
    if (mask == 0) {
      "WT"
    } else {
      mutated_sites <- which(as.integer(intToBits(mask))[1:K] == 1)
      paste0("M", paste(mutated_sites, collapse = ""))
    }
  })
  
  if(bitwise){
    return( setNames(rep(0, length(labels)), labels) )
  }else{
    
    # compute Hamming weights
    hamming_weights <- sapply(masks, function(mask) sum(as.integer(intToBits(mask))[1:K]))
    
    # order by weight, then lexicographically
    order_idx <- order(hamming_weights, labels)
    
    labels <- labels[order_idx]
    
    return(setNames(rep(0, length(labels)), labels))
  }
  
}


## simple expression for mu_del -- lethal deleterious
get_mu_del <- function(mu, genomeSize = NA, frac_del = 0.01){
  
  # if not provided, calculate appropriate genome size based on correlation w/ mu
  if(is.na(genomeSize)){
    # dat = read_excel("genome_mutation_corr.xlsx") # read in data
    # dat$genomeSize = dat$`Genome size (kb)`*1000 # convert to units bp
    # coef = lm(data = dat, formula = log10(genomeSize) ~ log10(mu)) 
    ## above regression used to obtain correlation model parameters
    #10^((log10(G) - 2.728987)/ (-0.2676165))
    genomeSize = floor( 10^(2.728987 - 0.2676165*log10(mu)) )
  }
  
  1 - ( 1 - mu )^( genomeSize * frac_del )
}

## simple expression for mu_ben -- resistance granting
#### N_B should be a vector of length K where the i-th element
#### gives the number of sites that grant resistance to drug i
get_mu_ben <- function(mu, N_B){
  1 - ( 1 - mu )^( N_B )
}

# --- helper: exponential range generator ---
exp_range <- function(start_exp, end_exp, decreasing = FALSE) {
  exps <- seq(start_exp, end_exp)
  vals <- as.vector(rbind(1*10^exps, 5*10^exps))
  if (decreasing) sort(vals, decreasing = TRUE) else sort(vals)
}

inv.logit = function(x){exp(x)/(1+exp(x))}
logit = function(x){log(x/(1-x))}

# --- helper to build parameter grid data frame ---
## as long as names match the model parameter names, any parameter combination can be used
mapVals <- function(param_list) {
  
  # --- helper to generate exponential ranges ---
  exp_range <- function(start_exp, end_exp, decreasing = FALSE) {
    exps <- seq(start_exp, end_exp)
    vals <- as.vector(rbind(1 * 10^exps, 5 * 10^exps))
    vals <- if (decreasing) sort(vals, decreasing = TRUE) else sort(vals, decreasing = FALSE)
    return(vals)
  }
  
  # --- default parameter ranges ---
  defaults <- list(
    n_drugs = 1:7,
    N_max   = exp_range(5, 10),
    k_drug  = c(1, 2, seq(5, 15, by = 5)),
    mu      = exp_range(-3, -6, decreasing = TRUE),
    b_fact  = c(1, seq(2, 8, by = 2))
  )
  
  # --- sanity checks ---
  stopifnot(is.list(param_list))
  if (is.null(names(param_list)) || any(names(param_list) == "")) {
    stop("All elements of param_list must be named.")
  }
  if (length(param_list) < 1) {
    stop("param_list must include at least one named parameter.")
  }
  
  # --- process parameters ---
  final_params <- list()
  for (nm in names(param_list)) {
    val <- param_list[[nm]]
    
    if (is.null(val)) {
      # ignore NULL parameters entirely
      next
    } else if (all(is.na(val))) {
      # replace NA parameters with defaults
      if (!is.null(defaults[[nm]])) {
        final_params[[nm]] <- defaults[[nm]]
      } else {
        stop(sprintf("No default defined for parameter '%s'.", nm))
      }
    } else {
      # use provided values as-is
      final_params[[nm]] <- val
    }
  }
  
  # --- expand over the parameters that remain ---
  if (length(final_params) == 0) {
    stop("No valid parameters supplied after ignoring NULL entries.")
  }
  
  df <- do.call(expand.grid, final_params)
  return(df)
}


fill_bounds <- function(val, lwr, upr, var) {
  if (!is.na(lwr) & !is.na(upr)) {
    return(list(lwr = lwr, upr = upr))
  }
  
  if (var %in% log_vars) {
    factor <- runif(1, 2, 4)  # span orders of magnitude
    return(list(
      lwr = val / factor,
      upr = val * factor
    ))
  } else {
    spread <- val * runif(1, 0.2, 0.5)
    return(list(
      lwr = val - spread,
      upr = val + spread
    ))
  }
}


sample_param <- function(var, n) {
  vmin <- min(extinctionMap[[var]], na.rm = TRUE)
  vmax <- max(extinctionMap[[var]], na.rm = TRUE)
  
  if (var %in% log_vars) {
    # uniform in log space
    return(10^runif(n, log10(vmin), log10(vmax)))
  } else {
    # uniform in linear space
    return(runif(n, vmin, vmax))
  }
}



## -----------------------------
## Latin hypercube generator
## -----------------------------
get_lh_pars <- function(N_LH_sets = 100, lh_pars, par_names = NULL, sampling_scheme = "unif", vary_single = 1){
  
  library(dplyr)
  library(lhs)
  library(purrr)
  
  # -----------------------------
  # 1. Valid parameters only
  # -----------------------------
  valid_names <- c("R0", "mu", "N_max", "w", "fit_cost", "n_drugs", "G")
  
  # if user specifies subset
  if (!is.null(par_names)) {
    valid_names <- intersect(valid_names, par_names)
  }
  
  # trim to valid columns only
  lh_pars <- lh_pars[, colnames(lh_pars) %in% valid_names, drop = FALSE]
  
  # safety checks
  if (nrow(lh_pars) == 0) {
    stop("No valid parameters found in lh_pars.")
  }
  
  if (any(is.na(lh_pars['lwr',])) | any(is.na(lh_pars['upr',]))) {
    stop("All parameters must have lwr and upr values.")
  }
  
  if (any(lh_pars['lwr',] >= lh_pars['upr',])) {
    stop("Each lwr must be strictly less than upr.")
  }
  
  # -----------------------------
  # 2. Generate normalized LHS
  # -----------------------------
  if(sampling_scheme == "lh"){
    
    lhs_raw <- optimumLHS(
      n = N_LH_sets,
      k = ncol(lh_pars)
    )
    
  }else{
    
    lhs_raw <- matrix(
      runif(N_LH_sets * ncol(lh_pars), min = 0, max = 1),
      nrow = N_LH_sets,
      ncol = ncol(lh_pars),
      dimnames = list(
        NULL,
        colnames(lh_pars)
      )
    )
    
  }
  colnames(lhs_raw) <- names(lh_pars)
  
  # -----------------------------
  # 3. Scale each parameter
  # -----------------------------
  scale_param <- function(x, lwr, upr, pname) {
    
    # log-scale biologically wide parameters
    if (pname %in% c("mu", "N_max", "G")) {
      return(10^(log10(lwr) + x * (log10(upr) - log10(lwr))))
    }
    
    # integer parameters
    if (pname %in% c("n_drugs")) {
      return(ceiling(lwr + x * (upr - lwr)))
    }
    
    # linear-scale parameters
    return(lwr + x * (upr - lwr))
  }
  
  lhs_scaled <- as.data.frame(lhs_raw)
  
  ## apply scaling to each parameter column
  for (i in 1:dim(lh_pars)[2]) {
    
    pname = names(lh_pars)[i]
    lhs_scaled[[pname]] <- scale_param(
      lhs_scaled[[pname]],
      lh_pars["lwr", pname],
      lh_pars["upr", pname],
      pname
    )
    
  }
  
  
  # -----------------------------
  # Filter / count biologically invalid sets
  # -----------------------------
  
  hard_valid <- lhs_scaled %>%
    rowwise() %>%
    mutate(
      hard_valid = valid_pars_lim(cur_data(), frac_del)
    ) %>%
    ungroup() %>%
    pull(hard_valid)
  
  print( paste0(sum(hard_valid), " of ", N_LH_sets, " parameter sets valid" ))
  
  ## keep only parameter sets that satisfy model conditions
  lhs_scaled <- as.data.frame(lhs_scaled[hard_valid,])
  lhs_raw <- as.data.frame(lhs_raw[hard_valid,])
  
  ## for each parameter set, generate 'vary_single' new sets for each parameter, one at a time
  if(vary_single > 1){
    oat_list <- vector("list", nrow(lhs_raw))
    for( i in 1:nrow(lhs_raw) ){
      
      base_row <- lhs_raw[i, , drop = FALSE]
      
      ## Generate OAT sets for each parameter
      row_expansions <- map_dfr(par_names, function(focal_par) {
        
        ## replicate base row
        raw_grid <- seq(from = round(min(lhs_raw[,focal_par]), digits = 3), 
                        to = round(max(lhs_raw[,focal_par]), digits = 3), 
                        by = 0.01)
        
        new_block <- as.data.frame(base_row[rep(1, length(raw_grid)), , drop = FALSE])
        
        ## vary only focal parameter
        #new_block[, focal_par] <- runif(vary_single)
        new_block[, focal_par] <- raw_grid
        
        # metadata
        new_block$base_id <- i
        new_block$base_par <- focal_par
        
        return(new_block)
      })
      
      # optional: include original baseline point
      base_with_meta <- base_row %>%
        mutate(
          base_id = i,
          base_par = "baseline"
        )
      
      oat_list[[i]] <- bind_rows(base_with_meta, row_expansions)
    }
    
    lhs_raw = bind_rows(oat_list) 
    rownames(lhs_raw) = NULL
    
    ## scale updated raw parameters if 'vary_single was used'
    lhs_scaled <- as.data.frame(lhs_raw)
    
    for (i in 1:dim(lh_pars)[2]) {
      
      pname = names(lh_pars)[i]
      lhs_scaled[[pname]] <- scale_param(
        lhs_scaled[[pname]],
        lh_pars["lwr", pname],
        lh_pars["upr", pname],
        pname
      )
      
    }
    
    # -----------------------------
    # Filter / count biologically invalid sets
    # -----------------------------
    
    hard_valid <- lhs_scaled %>%
      rowwise() %>%
      mutate(
        hard_valid = valid_pars_lim(cur_data(), frac_del)
      ) %>%
      ungroup() %>%
      pull(hard_valid)
    
    print( paste0(sum(hard_valid), " of ", length(hard_valid), " parameter sets valid" ))
    
    ## keep only parameter sets that satisfy model conditions
    lhs_scaled <- as.data.frame(lhs_scaled[hard_valid,])
    lhs_raw <- as.data.frame(t(lhs_raw[hard_valid,]))
    
  }
  
  
  # -----------------------------
  # Derived parameters
  # -----------------------------
  
  ## birth probability from R0
  lhs_scaled$b <- 1 - (1 / (1 + lhs_scaled$R0))
  
  ## lethal deleterious sites (Z)
  lhs_scaled$frac_del <- frac_del
  lhs_scaled$Z <- floor(lhs_scaled$G * lhs_scaled$frac_del)
  
  
  ## treatment population size
  lhs_scaled$N_maxTreat <- lhs_scaled$N_max * 2
  
  ## generation limits
  lhs_scaled$T_max <- ceiling(
    ceiling(log(lhs_scaled$N_max, base = 2 * lhs_scaled$b)) * 1.2
  )
  
  lhs_scaled$T_maxTreat <- max(
    ceiling(
      ceiling(log(lhs_scaled$N_maxTreat, base = 2 * lhs_scaled$b)) * 1.2
    ),
    100
  )
  
  return(lhs_scaled)
  
}


## valid parameter threshold functions
## ensures full resistant has positive growth and the drug is effective against WT
valid_pars_lim <- function(pars, frac_del){
  cond_res_grow = pars$R0 * ( (1-pars$mu)^(pars$G*frac_del) ) * (1-pars$fit_cost)^(pars$n_drugs) > 1 ## ensure full resistant genotype can grow
  cond_drug_kil =  pars$R0*(1-pars$w) < 1 ## ensure drug is effective in killing resistant strains
  return(cond_res_grow & cond_drug_kil) ## conditions that ensure no guranteed extinction / impossible resistance
}


## function for data downsampling and regression over y-logit space
process_param <- function(par_name, dat, mode = 0, model = "log",
                          minDat = 5, sigChange = 0.5) {
  
  tmp <- dat %>%
    filter( base_par %in% c(par_name, "baseline") ) %>%
    mutate( base_par = if_else(base_par == "baseline", par_name, base_par) )
  
  if (model == "log") { tmp[[par_name]] <- log10(tmp[[par_name]]) }
  
  #---- function to process each base_id ----#
  per_id <- function(i) {
    
    tmp_by_id <- tmp %>%
      filter(base_id == i) %>%
      mutate(pRes_trans = logit(fracResistance * 0.99 + 0.005))
    
    tmp_by_id <- switch(
      as.character(mode),
      "0" = { tmp_by_id },
      
      "1" = {
        drop_dat = tmp_by_id[!(tmp_by_id$fracResistance %in% c(0,1)),]
        if (nrow(drop_dat) >= minDat){ drop_dat }else{ NULL }
      },
      
      ## caution: can cause issues in rare cases of non-monotonicity in mu at boundary
      "2" = {
        drop_dat = tmp_by_id[!(tmp_by_id$fracResistance %in% c(0,1)),]
        x <- as.numeric(rownames(drop_dat))
        runs <- cumsum(c(TRUE, diff(x) != 1))
        x_keep <- x[runs %in% which(tabulate(runs) >= 3)]
        
        if (length(x_keep) == 0) {
          NULL
        } else {
          x_keep <- c(min(x_keep) - 1, x_keep, max(x_keep) + 1)
          
          # Keep only rows that actually exist
          x_keep <- intersect(x_keep, as.numeric(rownames(drop_dat)))
          
          drop_dat <- drop_dat[as.character(x_keep), , drop = FALSE]
          if (nrow(drop_dat) >= minDat){ drop_dat }else{ NULL }
        }
      },
      
      ## caution: can cause issues in rare cases of non-monotonicity in mu at boundary
      "3" = {
        drop_dat = tmp_by_id[!(tmp_by_id$fracResistance %in% c(0,1)),]
        x <- as.numeric(rownames(drop_dat))
        runs <- cumsum(c(TRUE, diff(x) != 1))
        x_keep <- x[runs %in% which(tabulate(runs) >= 3)]
        
        if (length(x_keep) == 0) {
          NULL
        } else {
          drop_dat <- drop_dat[as.character(x_keep), , drop = FALSE]
          if (nrow(drop_dat) >= minDat){ drop_dat }else{ NULL }
        }
      },
      
      "4" = {
        drop_dat = tmp_by_id[!(tmp_by_id$fracResistance %in% c(0,1)),]
        x <- as.numeric(rownames(drop_dat))
        runs <- cumsum(c(TRUE, diff(x) != 1))
        x_keep <- x[runs == which.max(tabulate(runs))]
        
        if (length(x_keep) == 0) {
          NULL
        } else {
          out <- drop_dat[as.character(x_keep), , drop = FALSE]
          if (nrow(out) >= minDat){ out }else{ NULL }
        }
      },
      
      "5" = {
        drop_dat = tmp_by_id
        
        if ( abs(max(drop_dat$fracResistance) - min(drop_dat$fracResistance)) < sigChange ){ 
          x_keep <- 1:dim(drop_dat)[1] ## keep all values if no strong change detected
        }else{
          x <- cumsum(sign(diff(drop_dat$fracResistance)))
          x_keep <- 1:( which(abs(x) == max(abs(x)))[1]+1 ) #cutoff at first stable value after change
        }
        
        ## trim data
        if (length(x_keep) == 0) {
          NULL
        } else {
          drop_dat <- drop_dat[as.character(x_keep), , drop = FALSE]
          drop_dat <- drop_dat[!(drop_dat$fracResistance %in% c(0,1)),]
          if (nrow(drop_dat) >= minDat){ drop_dat }else{ NULL }
        }
      },
      
      stop("Invalid mode")
    )
    
    if (is.null(tmp_by_id)) {
      return(list(
        reg = tibble( par_id = i, slope = NA, val_min = NA, val_max = NA, mean = NA),
        data = NULL
      ))
    }
    
    # Fit model
    fit <- lm( reformulate(par_name, response = "pRes_trans"), data = tmp_by_id )
    
    ## compute change and range metrics
    max_diff = (max(tmp_by_id$fracResistance) - min(tmp_by_id$fracResistance))
    if(max_diff != 0){
      alpha = 0.05
      partSwitch = tmp_by_id$fracResistance[1] %in% c(0,1) | tmp_by_id$fracResistance[dim(tmp_by_id)[1]] %in% c(0,1)
      fullSwitch = ( max(tmp_by_id$fracResistance) >= 1-alpha )  & ( min(tmp_by_id$fracResistance) <= alpha )
      run = cumsum(c(TRUE, diff(tmp_by_id$fracResistance) != 0))
      switchRange = tmp_by_id[c(max(which(run == min(run))), min(which(run == max(run)))), par_name]
      
    }else{ switchRange = c(NA,NA); partSwitch = F; fullSwitch = F }
    
    
    return(list(
      reg = tibble( par_id = i, slope = coef(fit)[[par_name]], diff = max_diff, diff_sign = max_diff*sign(coef(fit)[[par_name]]), val_min = switchRange[1], val_max = switchRange[2], mean = mean(tmp_by_id[,par_name]), ps = partSwitch, fs = fullSwitch),
      data = tmp_by_id
    ))
    
  }
  
  ## output parameter values should remain on natural scale
  if (model == "log") { tmp[[par_name]] <- 10^(tmp[[par_name]]) }
  
  res <- map(unique(tmp$base_id), per_id)
  
  list(
    reg = bind_rows(map(res, "reg")),
    data = bind_rows( map(res, "data"), .id = "base_id" )
  )
}

# -----------------------------
# Ellipse generator
# -----------------------------
make_ellipse <- function(cx, cy, lx, ux, ly, uy, n = 120,
                         log_x = FALSE, log_y = FALSE) {
  
  theta <- seq(0, 2*pi, length.out = n)
  
  if (log_x) {
    cx <- log10(cx)
    lx <- log10(lx)
    ux <- log10(ux)
  }
  
  if (log_y) {
    cy <- log10(cy)
    ly <- log10(ly)
    uy <- log10(uy)
  }
  
  rx <- (ux - lx) / 2
  ry <- (uy - ly) / 2
  
  x <- cx + rx * cos(theta)
  y <- cy + ry * sin(theta)
  
  if (log_x) x <- 10^x
  if (log_y) y <- 10^y
  
  tibble(x = x, y = y)
}

# Function to create ellipse in log space
ellipse_logspace <- function(N, mu, N_lwr, N_upr, mu_lwr, mu_upr, n = 150) {
  
  # transform to log10 space
  x0 <- log10(N)
  y0 <- log10(mu)
  
  rx <- max(abs(log10(N_lwr) - x0), abs(log10(N_upr) - x0))
  ry <- max(abs(log10(mu_lwr) - y0), abs(log10(mu_upr) - y0))
  
  t <- seq(0, 2 * pi, length.out = n)
  
  x_log <- x0 + rx * cos(t)
  y_log <- y0 + ry * sin(t)
  
  # convert back to original scale
  tibble(
    x = 10^x_log,
    y = 10^y_log
  )
}

# -----------------------------
# Plot function
# -----------------------------

K_max=max(threshold$drugs_needed)
my_colors <- setNames(
  viridis::viridis(K_max, option = "plasma"),
  as.character(1:K_max)
)

var_labels <- c(
  mu = expression("Mutation rate (" * mu * ")"),
  G = "Genome size (G)",
  w = expression("Drug killing rate (" * kappa * ")"),
  R0 = expression("Reproductive number (" * R[0] * ")"),
  N_max = expression("Population size ("*N[max]*")"),
  fit_cost = "Fitness cost (c)"
)


## Threshold heatmap plots
plot_surface <- function(pred_var, title_suffix) {
  
  # enforce discrete drug counts
  grid[[pred_var]] <- ceiling(grid[[pred_var]])
  
  # update breaks
  break_vals <- sort(unique(grid[[pred_var]]))
  threshold$drugs_needed_f <- factor(threshold$drugs_needed, levels = 1:K_max)
  
  legend_labels <- break_vals
  
  p <- ggplot() +
    geom_contour(
      data = grid,
      aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
        z = .data[[pred_var]]
      ),
      breaks = break_vals,
      color = "black",
      linewidth = 0.4
    ) +
    geom_contour_filled(
      data = grid,
      aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
        z = .data[[pred_var]]
      ),
      breaks = break_vals,
      alpha = 0.35
    ) +
    
    # random hypothetical data
    geom_polygon(
      data = ellipse_df,
      aes(x = x, y = y, group = id),
      alpha = 0.25,
      fill = "steelblue",
      color = "steelblue",
      inherit.aes = FALSE
    ) +
    geom_point(
      data = df_bounds,
      aes(x = .data[[x_var]], y = .data[[y_var]]),
      size = 3,
      color = "red"
    ) +
    geom_text(
      data = df_bounds,
      aes(x = .data[[x_var]], y = .data[[y_var]], label = name),
      vjust = -0.9,
      size = 3.8,
      color = "red"
    ) +
    
    # True pathogen data
    geom_polygon(
      data = path_dat_ellipse,
      aes(x = x, y = y, group = id),
      alpha = 0.25,
      fill = "steelblue",
      color = "steelblue",
      inherit.aes = FALSE
    ) +
    geom_point(
      data = path_dat_bounds,
      aes(x = .data[[x_var]], y = .data[[y_var]]),
      size = 3,
      color = "red"
    ) +
    geom_text(
      data = path_dat_bounds,
      aes(x = .data[[x_var]], y = .data[[y_var]], label = pathogen),
      vjust = -0.9,
      size = 3.8,
      color = "red"
    ) +
    
    scale_fill_manual(
      values = unname(my_colors[break_vals]),  
      labels = break_vals,
      name = "Drugs needed"
    ) +
    coord_cartesian(
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max),
      expand = FALSE
    ) +
    labs(
      x = x_var,
      y = y_var,
      title = paste0(
        "Drug threshold for resistance < ",
        p_crit,
        " (", title_suffix, ")"
      )
    ) +
    theme_bw(base_size = 12)
  
  if (log_x) p <- p + scale_x_log10()
  if (log_y) p <- p + scale_y_log10()
  
  return(p)
}

## alternative plot function using simulation data
plot_sim_surface <- function() {
  
  
  # plot
  threshold$drugs_needed_f <- factor(threshold$drugs_needed, levels = 1:max(threshold$drugs_needed))
  break_vals <- sort(unique(threshold$drugs_needed))
  p <- ggplot() + 
    geom_contour( 
      data = threshold, 
      aes( x = .data[[x_var]], y = .data[[y_var]], z = drugs_needed ), 
      breaks = break_vals, 
      color = "black", 
      linewidth = 0.5 ) + 
    
    geom_contour_filled( data = threshold, 
                         aes( x = .data[[x_var]], y = .data[[y_var]], z = drugs_needed ), 
                         breaks = break_vals, 
                         alpha = 0.35 ) +
    
    # # random hypothetical data
    # geom_polygon(
    #   data = ellipse_df,
    #   aes(x = x, y = y, group = id),
    #   alpha = 0.25,
    #   fill = "steelblue",
    #   color = "steelblue",
    #   inherit.aes = FALSE
    # ) +
    # geom_point(
    #   data = df_bounds,
    #   aes(x = .data[[x_var]], y = .data[[y_var]]),
    #   size = 3,
    #   color = "red"
    # ) +
    # geom_text(
    #   data = df_bounds,
    #   aes(x = .data[[x_var]], y = .data[[y_var]], label = name),
    #   vjust = -0.9,
    #   size = 3.8,
    #   color = "red"
    # ) +
    
    # True pathogen data
    geom_polygon(
      data = path_dat_ellipse,
      aes(x = x, y = y, group = id),
      alpha = 0.25,
      fill = "steelblue",
      color = "steelblue",
      inherit.aes = FALSE
    ) +
    geom_point(
      data = path_dat_bounds,
      aes(x = .data[[x_var]], y = .data[[y_var]]),
      size = 3,
      color = "red"
    ) +
    geom_text(
      data = path_dat_bounds,
      aes(x = .data[[x_var]], y = .data[[y_var]], label = pathogen),
      vjust = -0.9,
      size = 3.8,
      color = "red"
    ) +
    
    # discrete fill scale
    scale_fill_manual(
      values = unname(my_colors[break_vals]),  
      labels = break_vals,
      name = "Drugs needed"
    ) +
    
    # coordinate limits
    coord_cartesian(
      xlim = c(x_min, x_max),
      ylim = c(y_min, y_max),
      expand = FALSE
    ) +
    
    # labels and theme
    labs(
      x = var_labels[x_var],
      y = var_labels[y_var],
      title = paste0(
        "Drug threshold for resistance < ",
        p_crit,
        " (Simulation data)"
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      
      axis.text.x = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 14, color = "black"),
      
      legend.key.height = unit(0.6, "cm"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      plot.title = element_text(face = "bold")
    )
  
  # ---------------------------
  # 3. Log axes only for N_max and mu
  # ---------------------------
  if (log_x) p <- p + scale_x_log10()
  if (log_y) p <- p + scale_y_log10()
  p
  return(p)
}





## function to generate plot grids over parameter values
## y_var indicates the column name of the metric to plot
## facet row/column indicates variable to be displayed on the respective facets
## color var indicates variable name mapped to the color aesthetic 
plot_metric_grid <- function(df, y_var, x_var = "n_drugs",
                             facet_row = "N_max", 
                             facet_col = "mu", 
                             const_vals = NULL,
                             color_var = "mu",
                             transform = NULL) {
  
  # -- ensure all necessary columns exist --
  stopifnot(all(c(x_var, y_var, facet_row, facet_col, color_var) %in% names(df)))
  
  # --- Filter for constant variables if provided ---
  if (!is.null(const_vals)) {
    if (!is.list(const_vals)) {
      stop("`const_vals` must be a named list or named vector (e.g., list(alpha = 0.1, beta = 5)).")
    }
    
    for (var in names(const_vals)) {
      if (!var %in% names(df)) {
        warning(paste("Variable", var, "not found in data; skipping."))
        next
      }
      val <- const_vals[[var]]
      df <- df[df[[var]] == val, , drop = FALSE]
    }
    
    if (nrow(df) == 0) {
      stop("No rows found matching the provided const_vals filters.")
    }
  }
  
  # --- apply optional transformation to y_var ---
  if (!is.null(transform) && nzchar(transform)) {
    if (!transform %in% c("log", "sqrt")) {
      stop("Unsupported transform: use 'log', 'sqrt', or leave blank.")
    }
    yvals <- df[[y_var]]
    
    if (transform == "log") {
      # Apply natural log and replace illegal results with NA
      transformed <- suppressWarnings(log(yvals, base = 10))
      invalid <- !is.finite(transformed)
      if (any(invalid, na.rm = TRUE)) {
        warning(sum(invalid), " invalid y-values (<= 0 or NaN) replaced with NA during log-transform.")
        transformed[invalid] <- NA
      }
      df[[y_var]] <- transformed
    } else if (transform == "sqrt") {
      if (any(df[[y_var]] < 0, na.rm = TRUE)) {
        warning("Negative values found in y_var; sqrt-transform may produce NaN.")
      }
      df[[y_var]] <- sqrt(df[[y_var]])
    }
  }
  
  # --- convert and prepare numeric versions of the chosen parameters ---
  df <- df %>%
    mutate(
      across(all_of(c(facet_row, facet_col, color_var)), 
             ~ as.numeric(as.character(.x)), 
             .names = "{.col}_num")
    )
  
  # --- unique sorted values for ordering ---
  vals <- function(v, decreasing = FALSE) sort(unique(v), decreasing = decreasing)
  
  row_vals  <- vals(df[[paste0(facet_row, "_num")]], decreasing = FALSE)
  col_vals  <- vals(df[[paste0(facet_col, "_num")]], decreasing = TRUE)
  colr_vals <- vals(df[[paste0(color_var, "_num")]], decreasing = TRUE)
  
  # --- label formatting ---
  fmt <- function(x) sprintf("%.0e", x)
  
  # --- variable name → plotmath-compatible label mapping ---
  tex_labels <- c(
    k_drug = "kappa[drug]",
    N_max  = "N[T]",
    mu     = "mu",
    b_fact = "b"
  )
  
  get_tex_label <- function(var) {
    if (var %in% names(tex_labels)) tex_labels[[var]] else var
  }
  
  # --- formatted factor labels for facets ---
  df <- df %>%
    mutate(
      !!paste0(facet_row, "_lab") := factor(
        paste0(get_tex_label(facet_row), " == ", fmt(.data[[paste0(facet_row, "_num")]])),
        levels = paste0(get_tex_label(facet_row), " == ", fmt(row_vals))
      ),
      !!paste0(facet_col, "_lab") := factor(
        paste0(get_tex_label(facet_col), " == ", fmt(.data[[paste0(facet_col, "_num")]])),
        levels = paste0(get_tex_label(facet_col), " == ", fmt(col_vals))
      ),
      color_factor = factor(.data[[paste0(color_var, "_num")]], levels = colr_vals)
    )
  
  # --- palette ---
  n_col_levels <- length(colr_vals)
  pal <- rev(sequential_hcl(n_col_levels + 2, palette = "Heat"))[2:(n_col_levels + 1)]
  
  # --- axis labels ---
  y_lab <- gsub("_", " ", y_var)
  if (!is.null(transform) && nzchar(transform)) {
    y_lab <- paste0(y_lab, " (", transform, "-scale)")
  }
  x_lab <- expression("Number of Drugs (" * N[drugs] * ")") ## ***fix*** ##
  
  # --- build plot ---
  p <- ggplot(df,
              aes(x = .data[[x_var]],
                  y = .data[[y_var]],
                  color = color_factor,
                  group = color_factor)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = pal, name = color_var) +
    facet_grid(
      rows = vars(!!sym(paste0(facet_row, "_lab"))),
      cols = vars(!!sym(paste0(facet_col, "_lab"))),
      labeller = label_parsed
    ) +
    theme_bw() +
    labs(
      x = x_lab,
      y = y_lab,
      subtitle = paste(
        c(
          if (!is.null(const_vals))
            paste(paste0(names(const_vals), " = ", unlist(const_vals)), collapse = ", "),
          if (!is.null(transform) && nzchar(transform))
            paste0("Transformed: ", transform)
        ),
        collapse = "; "
      )
    ) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      #legend.position = "none",
      strip.text = element_text(size = 11, face = "bold"),
      panel.spacing = unit(0.8, "lines")
    )
  
  return(p)
}


# function to generate heatmap over parameter values
## fill_val indicates the column name of the metric to plot
## x/y_var indicates variables to be displayed on the respective axes
## const_vals is necessary for parameter grids over 3 or more parameters 
## param_names provide 
plot_metric_heatmap <- function(df, fill_val,
                                x_var = "param_x",
                                y_var = "param_y",
                                const_vals = NULL,
                                param_names = NULL) {
  # Convert names to symbols for tidy evaluation
  x_param <- sym(x_var)
  y_param <- sym(y_var)
  fill_param <- sym(fill_val)
  
  # --- Filter for constant variables if provided ---
  if (!is.null(const_vals)) {
    if (!is.list(const_vals)) {
      stop("`const_vals` must be a named list or named vector (e.g., list(alpha = 0.1, beta = 5)).")
    }
    
    for (var in names(const_vals)) {
      if (!var %in% names(df)) {
        warning(paste("Variable", var, "not found in data."))
        next
      }
      val <- const_vals[[var]]
      df <- df[df[[var]] == val, , drop = FALSE]
    }
    
    if (nrow(df) == 0) {
      stop("No rows found matching the provided const_vals filters.")
    }
  }
  
  # --- Use param_names for axes labels if provided ---
  if (is.null(param_names)) {param_names = c(x_var, y_var)}
  
  
  # --- Compute label color contrast ---
  df$label_color <- ifelse(df[[as.character(fill_param)]] < 0.5, "white", "black")
  
  # --- Build heatmap ---
  p <- ggplot(df, aes(x = !!x_param, y = !!y_param)) +
    geom_tile(aes(fill = !!fill_param), color = "grey80") +
    geom_text(aes(label = sprintf("%.2f", !!fill_param), color = label_color),
              size = 4, fontface = "bold") +
    scale_fill_viridis(option = "B", limits = c(0, 1), direction = 1) +
    scale_color_identity() +
    labs(
      x = param_names[1],
      y = param_names[2],
      fill = "Metric"
    ) +
    coord_cartesian(expand = FALSE) +
    theme_bw()
  
  print(p)
}




plot_monotonicity <- function(
    sens_results,
    sens_vars,
    transform = '',
    output_var = "fracResistance",
    save_plots = TRUE,
    out_dir = "./figures/",
    tag = "",
    width = 840,
    height = 840,
    res = 100
){
  
  # -----------------------------
  # Safety checks
  # -----------------------------
  if (!(output_var %in% names(sens_results))) {
    stop(paste("Output variable", output_var, "not found in sens_results"))
  }
  
  missing_vars <- sens_vars[!sens_vars %in% names(sens_results)]
  if (length(missing_vars) > 0) {
    stop(
      paste(
        "Missing sensitivity variables:",
        paste(missing_vars, collapse = ", ")
      )
    )
  }
  
  # -----------------------------
  # Helper for log-scale variables
  # -----------------------------
  log_vars <- c("mu", "N_max", "G")
  
  
  # -----------------------------
  # Storage for plots
  # -----------------------------
  plot_list <- list()
  
  # -----------------------------
  # Generate one plot per variable
  # -----------------------------
  for (var in sens_vars) {
    
    fit <- lm(reformulate(var, response = output_var), data = sens_results)
    slope <- coef(fit)[[var]]
    
    p <- ggplot( data = sens_results, aes(x = .data[[var]], y = .data[[output_var]]) ) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", formula = y ~ x, color = "blue") +
      labs(
        x = var,
        y = output_var,
        title = paste("Monotonicity:", output_var, "vs", var)
      ) +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = paste0("Slope = ", round(slope, 3)),
        hjust = 1.1,
        vjust = 1.5
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid = element_blank()
      )
    
    # apply log scale if appropriate
    if (var %in% log_vars | transform == 'log') {
      p <- p + scale_x_log10()
    }
    
    # save individually if requested
    if (save_plots) {
      jpeg(
        filename = file.path(
          out_dir,
          paste0("monotonic_", var, tag, ".jpeg")
        ),
        width = width,
        height = height,
        units = "px",
        res = res
      )
      print(p)
      dev.off()
    }
    
    plot_list[[var]] <- p
  }
  
  # -----------------------------
  # Combined panel plot -- 2 plots per row
  # -----------------------------
  row_tot = ceiling( length(plot_list)/2 )
  combined_plot <- marrangeGrob(
    grobs = plot_list,
    nrow = row_tot,
    ncol = 2,
    top = paste("Sensitivity monotonicity:", output_var)
  )
  
  # optional combined PDF
  if (save_plots) {
    ggsave(
      filename = file.path(
        out_dir,
        paste0("monotonic_ALL_", output_var, tag, ".pdf")
      ),
      combined_plot,
      width = 12,
      height = 10
    )
  }
  
  # -----------------------------
  # Return plots
  # -----------------------------
  return(plot_list)
}



