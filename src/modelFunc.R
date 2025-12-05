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
library(colorspace)
library(data.table)
library(doParallel)
library(foreach)


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
              parameters$n_retry
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


# --- helper: exponential range generator ---
exp_range <- function(start_exp, end_exp, decreasing = FALSE) {
  exps <- seq(start_exp, end_exp)
  vals <- as.vector(rbind(1*10^exps, 5*10^exps))
  if (decreasing) sort(vals, decreasing = TRUE) else sort(vals)
}


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
      legend.position = "none",
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

