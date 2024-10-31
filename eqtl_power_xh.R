#!/usr/bin/env Rscript

# eQTL Power Analysis Tool
# Author: XH and claud.ai
# Description: Calculate and visualize statistical power for eQTL studies
#              in both bulk tissue and cell-type specific contexts

# Function to print help message
print_help <- function() {
    cat("eQTL Power vs PVE Analysis Tool\n")
    cat("Usage: Rscript eqtl_power_pve_xh.R [parameters]\n\n")
    cat("Shared Parameters (optional):\n")
    cat("  pval_thr=<number>    P-value threshold (default: 1e-3)\n")
    cat("  n=<integer>          Sample size (default: 500)\n")
    cat("  output=<string>      Output figure filename (optional)\n\n")
    cat("Bulk Analysis Parameter:\n")
    cat("  pve=<string>         Single PVE value or comma-separated PVE values\n")
    cat("                       (default: sequence from 0.01 to 0.05 by 0.005)\n\n")
    cat("Cell-type Specific Analysis Parameters:\n")
    cat("  cell_prop=<string>   Comma-separated cell proportions\n")
    cat("                       Providing this triggers cell-type analysis\n")
    cat("                       (default: 0.3,0.2,0.2,0.1,0.1,0.05,0.05)\n")
    cat("  se=<string>          Comma-separated standard errors (default: 1 for each cell)\n")
    cat("  beta_focal=<string>  Comma-separated effect sizes for focal cells\n")
    cat("                       (default: sequence from sqrt(0.01) to sqrt(0.1))\n\n")
    cat("Examples:\n")
    cat("  # Bulk analysis with single PVE:\n")
    cat("  Rscript eqtl_power_pve_xh.R pval_thr=1e-4 n=1000 pve=0.03\n\n")
    cat("  # Bulk analysis with multiple PVE values and plot:\n")
    cat("  Rscript eqtl_power_pve_xh.R pval_thr=1e-4 n=1000 pve=0.01,0.02,0.03 output=power.pdf\n\n")
    cat("  # Cell-type specific analysis:\n")
    cat("  Rscript eqtl_power_pve_xh.R pval_thr=1e-4 n=1000 cell_prop=0.4,0.3,0.3 output=cell_power.pdf\n")
}

# Function to print summary results
print_summary <- function(PVE, power, power.reduced.50 = NULL, power.reduced.30 = NULL) {
    cat("\nPower Analysis Results:\n")
    cat("---------------------\n")
    for(i in 1:length(PVE)) {
        cat(sprintf("PVE = %.3f:\n", PVE[i]))
        cat(sprintf("  Power (pi = 1.0): %.3f\n", power[i]))
        if(!is.null(power.reduced.50)) {
            cat(sprintf("  Power (pi = 0.5): %.3f\n", power.reduced.50[i]))
        }
        if(!is.null(power.reduced.30)) {
            cat(sprintf("  Power (pi = 0.3): %.3f\n", power.reduced.30[i]))
        }
        cat("\n")
    }
}

# Function for basic power analysis
run_basic_power_analysis <- function(pval_thr, n, PVE) {
    z.thr <- qnorm(1-pval_thr/2)
    z <- sqrt(PVE * n)
    
    power <- 1 - pnorm(z.thr, z, 1)
    power.reduced.50 <- 1 - pnorm(z.thr, z * 0.5, 1)
    power.reduced.30 <- 1 - pnorm(z.thr, z * 0.3, 1)
    
    return(list(PVE = PVE, 
                power = power,
                power.reduced.50 = power.reduced.50,
                power.reduced.30 = power.reduced.30))
}

# Function for cell-type specific pure power analysis
est_power_pure <- function(beta.focal, p, sigma, n, pval_thr) {
    z.thr <- qnorm(1-pval_thr/2)
    PVE.pure <- rep(0, length(beta.focal))
    power.pure <- rep(0, length(beta.focal))
    
    for (i in 1:length(beta.focal)) {
        PVE.pure[i] <- beta.focal[i]^2 / (beta.focal[i]^2 + sigma[1]^2)
        z.pure <- sqrt(PVE.pure[i] * n)
        power.pure[i] <- 1 - pnorm(z.thr, z.pure, 1)
    }
    
    return(list(PVE = PVE.pure, power = power.pure))
}

# Function for cell-type specific bulk power analysis
est_power_bulk <- function(beta.focal, p, sigma, n, focal, pval_thr) {
    z.thr <- qnorm(1-pval_thr/2)
    beta <- rep(0, length(p))
    PVE.bulk <- rep(0, length(beta.focal))
    power.bulk <- rep(0, length(beta.focal))
    
    for (i in 1:length(beta.focal)) {
        PVE.bulk[i] <- p[focal]^2 * beta.focal[i]^2 / 
            (p[focal]^2 * beta[focal]^2 + sum(p^2 * sigma^2))
        z.bulk <- sqrt(PVE.bulk[i] * n)
        power.bulk[i] <- 1 - pnorm(z.thr, z.bulk, 1)
    }
    
    return(list(PVE = PVE.bulk, power = power.bulk))
}

# Parse comma-separated string to numeric vector
parse_vector <- function(str) {
    if (is.null(str)) return(NULL)
    as.numeric(unlist(strsplit(str, ",")))
}

# Reconcile this model:  https://bwhbioinfo.shinyapps.io/powerEQTL/
compute_pve_from_slope <- function(slope, MAF, sigma.y) {
    sigma2.x = 2 * MAF * (1 - MAF)
    pve = (slope^2 * sigma2.x) / (sigma.y^2)
    return(pve)
}

# Main execution function
main <- function() {
    # Parse command line arguments
    args <- commandArgs(TRUE)
    if (length(args) == 0) {
        print_help()
        return()
    }
    
    # Set default shared parameters
    params <- list(
        pval_thr = 1e-3,
        n = 500
    )
    
    # Parse parameters
    for (arg in args) {
        split_arg <- strsplit(arg, "=")[[1]]
        if (length(split_arg) == 2) {
            param_name <- split_arg[1]
            param_value <- split_arg[2]
            
            # Convert numeric parameters
            if (param_name %in% c("pval_thr", "n")) {
                param_value <- as.numeric(param_value)
            }
            
            params[[param_name]] <- param_value
        }
    }
    
    # Check if cell-type specific analysis is requested
    is_cell_specific <- !is.null(params$cell_prop)
    
    if (!is_cell_specific) {
        # Process PVE values
        if (!is.null(params$pve)) {
            PVE <- parse_vector(params$pve)
        } else {
            PVE <- seq(0.01, 0.05, 0.005)  # Default sequence
        }
        
        # Run bulk power analysis
        results <- run_basic_power_analysis(params$pval_thr, params$n, PVE)
        
        # Print summary
        cat(sprintf("\nParameters: n = %d, p-value threshold = %g\n", 
                   params$n, params$pval_thr))
        print_summary(results$PVE, results$power, results$power.reduced.50, 
                     results$power.reduced.30)
        
        # Generate plot if output is specified and there's more than one PVE value
        if (!is.null(params$output) && length(PVE) > 1) {
            pdf(file = params$output, height = 3, width = 4)
            plot(results$PVE, results$power, type = "l", col = "red",
                 xlab = "PVE", ylab = "Power", ylim = c(0, 1))
            lines(results$PVE, results$power.reduced.50, type = "l", col = "blue")
            lines(results$PVE, results$power.reduced.30, type = "l", col = "green")
            legend("topleft", c("pi = 1", "pi = 0.5", "pi = 0.3"),
                   col = c("red", "blue", "green"), lty = 1, cex = 1.0)
            dev.off()
            cat(sprintf("\nPlot saved to: %s\n", params$output))
        }
        
    } else {
        # Set up cell-type specific parameters
        p <- parse_vector(params$cell_prop)
        if (is.null(p)) {
            p <- c(0.3, 0.2, 0.2, 0.1, 0.1, 0.05, 0.05)
        }
        
        sigma <- if (!is.null(params$se)) {
            parse_vector(params$se)
        } else {
            rep(1, length(p))
        }
        
        beta.focal <- if (!is.null(params$beta_focal)) {
            parse_vector(params$beta_focal)
        } else {
            seq(sqrt(0.01), sqrt(0.1), 0.02)
        }
        
        # Calculate powers
        result.pure <- est_power_pure(beta.focal, p, sigma, params$n, params$pval_thr)
        result.bulk.high <- est_power_bulk(beta.focal, p, sigma, params$n, 1, params$pval_thr)
        result.bulk.low <- est_power_bulk(beta.focal, p, sigma, params$n, 2, params$pval_thr)
        
        # Print summary
        cat(sprintf("\nParameters: n = %d, p-value threshold = %g\n", 
                   params$n, params$pval_thr))
        cat(sprintf("Cell proportions: %s\n", 
                   paste(sprintf("%.3f", p), collapse = ", ")))
        cat("\nPower Analysis Results:\n")
        cat("---------------------\n")
        for(i in 1:length(beta.focal)) {
            cat(sprintf("Effect size (beta) = %.3f:\n", beta.focal[i]))
            cat(sprintf("  Power (pure): %.3f\n", result.pure$power[i]))
            cat(sprintf("  Power (bulk, %.1f%%): %.3f\n", 100*p[1], result.bulk.high$power[i]))
            cat(sprintf("  Power (bulk, %.1f%%): %.3f\n", 100*p[2], result.bulk.low$power[i]))
            cat("\n")
        }
        
        # Generate plot if output is specified
        if (!is.null(params$output)) {
            pdf(file = params$output, height = 5, width = 5)
            plot(result.pure$PVE, result.pure$power, type = "l", col = "red",
                 ylim = c(0,1), xlab = "Effect size (PVE scale)", ylab = "Power")
            lines(result.pure$PVE, result.bulk.high$power, type = "l", col = "blue")
            lines(result.pure$PVE, result.bulk.low$power, type = "l", col = "green")
            legend("topleft", c("pure", sprintf("%.1f%%", 100*p[1]), 
                               sprintf("%.1f%%", 100*p[2])),
                   col = c("red", "blue", "green"), lty = 1, cex = 1.0)
            dev.off()
            cat(sprintf("\nPlot saved to: %s\n", params$output))
        }
    }
}

# Run main function
main()
