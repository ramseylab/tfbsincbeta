## This file is part of the "tfbsincbeta.R" software program. "tfbsincbeta.R" is free
## software: you can redistribute it and/or modify it under the terms of
## Apache Software License version 2.0 (and incorporated into this package
## distribution in the file LICENSE). The R script reads a data file "Matrices.txt"
## which is in tab-delimited text format (see file ../data/Matrices.txt for an example).
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE.  See the Apache Software License version 2.0 for 
## details.
##
## Copyright Stephen Ramsey, Oregon State University
## 2016.09
##

library(parallel)

set.seed(12345)

## DEFINE GLOBAL VARIABLES:
g_fdtot <- rep(100000, 4)
g_fdtot_sum <- sum(g_fdtot)
g_rseq_len <- 20000
g_rseq_base <- sample(1:4,
                    prob=g_fdtot/(sum(g_fdtot)),
                    size=10000,
                    replace=TRUE)
g_lambda_max <- 0.001
g_alpha <- 1
g_nu <- 10000
g_mcmc_burn_in_steps <- 5000
g_mcmc_steps_per_sample <- 100
g_mcmc_num_samples <- 10
g_mcmc_weight_exponent <- 0.85

### read TFBS PCMs:
g_tfbs_pcms_all <- read.table("../data/Matrices.txt",
                              header=FALSE,
                              stringsAsFactors=FALSE,
                              sep="\t",
                              quote="",
                              comment.char="")
names(g_tfbs_pcms_all) <- c("PCM_NAME","A","C","G","T")

g_unique_tfbs_names <- unique(g_tfbs_pcms_all[,1])
g_num_tfbs_names_to_test <- 100
g_num_sites <- 1:10
g_num_trials <- 120

g_tfbs_pcms_list <- lapply(g_unique_tfbs_names,
                           function(tfbs_name) {
                               cmat <- g_tfbs_pcms_all[which(g_tfbs_pcms_all[,1] == tfbs_name), 2:5]
                               if (any(apply(cmat, 1, sum) <= 1)) {
                                   return(NULL)
                               }
                               cmat
                           })
names(g_tfbs_pcms_list) <- g_unique_tfbs_names
g_tfbs_pcms_list <- Filter(Negate(is.null), g_tfbs_pcms_list)

g_tfbs_pcms_try_list <- g_tfbs_pcms_list[sample(1:length(g_tfbs_pcms_list),
                                                replace=FALSE,
                                                size=g_num_tfbs_names_to_test)]

## logarithm of the incomplete beta function
##
## Test results:
##   > g_compute_log_incomplete_beta_function(0, 1, 1)
##   [1] -Inf
##   > g_compute_log_incomplete_beta_function(1, 1, 1)
##   [1] 0
##   > g_compute_log_incomplete_beta_function(1, 1, 10)
##   [1] -2.302585
##   > g_compute_log_incomplete_beta_function(1, 1, 20)
##   [1] -2.995732
##   > g_compute_log_incomplete_beta_function(1, 1, 50)
##   [1] -3.912023
##   > g_compute_log_incomplete_beta_function(1, 50, 50)
##   [1] -70.00272
g_compute_log_incomplete_beta_function <- function(x, p, w) {
  pbeta(x, p, w, log=TRUE) + lbeta(p, w)
}

## computes the log of the ratio of the TFBS model to the background model, for a sequence "lseq"
## given a PFM "cmat + sigmamat"

### lseq: integer l-vector of the sequence within the binding site
### cmat: lx4 position-count matrix for representative binding sites for the TF
### sigmamat: lx4 position-count matrix for putative binding sites for the TF
### fdtot:  integer 4-vector of nucleotide counts for sequence outside of binding sites
g_compute_log_prob_ratio_singlesite_seqdeppart <- function(lseq,
                                                           cmat,
                                                           sigmamat,
                                                           fdtot) {
    l <- nrow(cmat)
    seq_inds <- cbind(1:l, lseq)
    ndbl <- tabulate(lseq, 4)

    sum(log(cmat[seq_inds] + sigmamat[seq_inds] + 1)) - sum(ndbl*log(fdtot))
}

### cmat_sum: integer l-vector of the row-sums of the "cmat" matrix (should be 
###           less than or equal to the number of representative binding sites)
### total_length: the total length of sequence being analyzed
### beta_l1:  the number of putative binding sites
g_compute_log_prob_ratio_seq_indep_part_unif_lambda <- function(cmat_sum,
                                                              total_length,
                                                              beta_l1) {
    l <- length(cmat_sum) 

    -sum(log(cmat_sum + beta_l1 + 4)) + l*log(total_length + 3 - beta_l1*l)
}


g_compute_log_prob_ratio_change_num_binding_sites_tfbsincbeta <- function(l,
                                                                     beta_l1,
                                                                     lambda_max,
                                                                     p_alpha,
                                                                     p_nu,
                                                                     Lr) {
    -log(2) +
         g_compute_log_incomplete_beta_function(lambda_max,
                                                p_alpha + beta_l1 + 1,
                                                p_nu + Lr - (2*l - 1)*(beta_l1 + 1)) -
         g_compute_log_incomplete_beta_function(lambda_max,
                                                p_alpha + beta_l1,
                                                p_nu + Lr - (2*l - 1)*(beta_l1)) 
 
}

g_compute_log_prob_ratio_seq_indep_part <- function(cmat_sum,
                                                  total_length,
                                                  beta_l1,
                                                  Lr,
                                                  pf_compute_log_prob_ratio_change_num_binding_sites) {
    
   g_compute_log_prob_ratio_seq_indep_part_unif_lambda(cmat_sum,
                                                       total_length,
                                                       beta_l1) +
       pf_compute_log_prob_ratio_change_num_binding_sites(beta_l1, Lr)
}
                                                  
g_compute_mcmc_pos_weights_and_first_pass_seq <- function(cmat,
                                                        cmat_sum,
                                                        rseq,
                                                        fdtot,
                                                        total_length,
                                                        Lr,
                                                        mcmc_weight_exponent,
                                                        pf_compute_log_prob_ratio_change_num_binding_sites) {
    
    B <- length(rseq)
    l <- nrow(cmat)    

    beta_l1 <- 0
    log_prob_ratio_seq_indep_part <- g_compute_log_prob_ratio_seq_indep_part(cmat_sum,
                                                                           total_length,
                                                                           beta_l1,
                                                                           B,
                                                                           pf_compute_log_prob_ratio_change_num_binding_sites)
    
    sigmamat <- matrix(rep(0, l*4), ncol=4)
    
    log_prob_ratios_fwd <- log_prob_ratio_seq_indep_part +
                    c(sapply( 1:(B-l+1), function(b) {g_compute_log_prob_ratio_singlesite_seqdeppart(rseq[b:(b+l-1)],
                                                                    cmat,
                                                                    sigmamat,
                                                                    fdtot)} ),
                      rep(-Inf, l-1))

    log_prob_ratios_rev <- log_prob_ratio_seq_indep_part +
                    c(rep(-Inf, l-1),
                      sapply( l:B, function(b) {g_compute_log_prob_ratio_singlesite_seqdeppart(5 - rseq[b:(b-l+1)],
                                                                    cmat,
                                                                    sigmamat,
                                                                    fdtot)} ))
    p0 <- 1/(1 + exp(log_prob_ratios_fwd) + exp(log_prob_ratios_rev))
    pfwd <- p0*exp(log_prob_ratios_fwd)
    prev <- p0*exp(log_prob_ratios_rev)
    probmat <- cbind(prev, p0, pfwd)
    beta_vec_values <- c(-1, 0, 1)

    first_pass_seq <- order(apply(probmat, 1, function(prob_vec) {max(prob_vec[c(1,3)])}), decreasing=TRUE)
    
    entrop0 <- -p0*log(p0)
    entrop0[is.nan(entrop0)] <- 0
    entrop_fwd <- -pfwd*log(pfwd)
    entrop_fwd[is.nan(entrop_fwd)] <- 0
    entrop_rev <- -prev*log(prev)
    entrop_rev[is.nan(entrop_rev)] <- 0
    entropies <- entrop0 + entrop_fwd + entrop_rev
    weights <- entropies^mcmc_weight_exponent
    weights <- weights/sum(weights)

    list(pos_weights=weights,
         first_pass_seq=first_pass_seq)
}

g_pick_outcome_mcmc_step_given_b <- function(cmat,
                                           cmat_sum,
                                           sigmamat,
                                           rseq,
                                           fdtot,
                                           total_length,
                                           beta_vec,
                                           beta_l1,
                                           b,
                                           unif_sample,
                                           Lr,
                                           is_initialization_step,
                                           pf_compute_log_prob_ratio_change_num_binding_sites) {

    l <- nrow(cmat)
    B <- length(rseq)

    beta_vals <- -1:1
    
    curbeta_b <- beta_vec[b]
    
    log_prob_no_site_plus_seq_indep_part <- g_compute_log_prob_ratio_seq_indep_part(cmat_sum,
                                                                           total_length,
                                                                           beta_l1,
                                                                           B,
                                                                           pf_compute_log_prob_ratio_change_num_binding_sites)

    if (b < l) {
        allowed_rev <- FALSE
    }
    else {
        if (any(beta_vec[(b-l+1):(b-1)] != 0)) {
            allowed_rev <- FALSE
        }
        else {
            if (b > l && any(beta_vec[max(1, b - 2*l + 1):(b - l)] > 0)) {
                allowed_rev <- FALSE
            }
            else {
                allowed_rev <- TRUE
            }
        }
    }

    if (b > B - l + 1) {
        allowed_fwd <- FALSE
    }
    else {
        if (any(beta_vec[(b+1):(b+l-1)] != 0)) {
            allowed_fwd <- FALSE
        }
        else {
            if (b <= B - l && any(beta_vec[(b+l):min(B, b + 2*l - 1)] < 0)) {
                allowed_fwd <- FALSE
            }
            else {
                allowed_fwd <- TRUE
            }
        }        
    }

    ## there are two possible new states at this position b; choose one of them with equal probability
    newbeta_b <- sample(setdiff(beta_vals, curbeta_b), size=1)
    if (curbeta_b == 0) {
        adding <- TRUE
        orientation <- newbeta_b
        if (newbeta_b == -1) {
            allowed <- allowed_rev
        }
            else {
                allowed <- allowed_fwd
            }
    }
    else {
        adding <- FALSE
        allowed=TRUE
        orientation <- curbeta_b
    }
    
    if (! is_initialization_step && ! allowed) {
        newbeta_b <- curbeta_b
    }
    else {
        if (allowed_fwd) {
            log_prob_ratio_site_vs_none_fwd <- g_compute_log_prob_ratio_singlesite_seqdeppart(rseq[b:(b+l-1)],
                                                                                              cmat,
                                                                                              sigmamat,
                                                                                              fdtot) +
                                               log_prob_no_site_plus_seq_indep_part 
        }
        else {
            log_prob_ratio_site_vs_none_fwd <- -Inf
        }

        if (allowed_rev) {
            log_prob_ratio_site_vs_none_rev <- g_compute_log_prob_ratio_singlesite_seqdeppart(5 - rseq[b:(b-l+1)],
                                                                                              cmat,
                                                                                              sigmamat,
                                                                                              fdtot) +
                                               log_prob_no_site_plus_seq_indep_part 
        }
        else {
            log_prob_ratio_site_vs_none_rev <- -Inf
        }
        
        if (is_initialization_step) {
            # just choose the most probable state
            ef <- exp(log_prob_ratio_site_vs_none_fwd)
            er <- exp(log_prob_ratio_site_vs_none_rev)
            p0 <- 1/(1 + ef + er)
            pf <- p0 * ef
            pr <- p0 * er

            # just choose the most probable state, no need to randomize
            newbeta_b <- beta_vals[which.max(c(pr, p0, pf))]
        }
        else {
            if (orientation == 1) {
                log_prob_ratio_site_vs_none <- log_prob_ratio_site_vs_none_fwd
            }
            else {
                log_prob_ratio_site_vs_none <- log_prob_ratio_site_vs_none_rev
            }
            if (adding) {
                alpha <- exp(log_prob_ratio_site_vs_none)
            }
            else {
                alpha <- exp(-log_prob_ratio_site_vs_none)
            }
            if (unif_sample >= min(alpha, 1)) {
                newbeta_b <- curbeta_b
            }
        }
    }
    newbeta_b
}

g_sample_posterior_beta_using_mcmc <- function(cmat,
                                             rseq,
                                             fdtot,
                                             mcmc_weight_exponent,
                                             mcmc_burn_in_steps,
                                             mcmc_steps_per_sample,
                                             mcmc_num_samples,
                                             pf_compute_log_prob_ratio_change_num_binding_sites) {

    cmat_sum <- apply(cmat, 1, sum)
    
    B <- length(rseq)
    l <- nrow(cmat)
    sigmamat <- matrix(rep(0, l*4), ncol=4)

    beta_l1 <- 0

    total_length <- sum(fdtot) + length(rseq)
    
    ### compute weights for each of the base positions within the sequence
    ret_list <- g_compute_mcmc_pos_weights_and_first_pass_seq(cmat,
                                                            cmat_sum,
                                                            rseq,
                                                            fdtot,
                                                            total_length,
                                                            B,
                                                            mcmc_weight_exponent,
                                                            pf_compute_log_prob_ratio_change_num_binding_sites)
    pos_weights <- ret_list$pos_weights
    stopifnot(all(pos_weights > 0))
    first_pass_seq <- ret_list$first_pass_seq  # this is the vector of locations of most probable binding sites

    mcmc_total_steps <- B + mcmc_burn_in_steps + mcmc_num_samples*mcmc_steps_per_sample

    ### generate a random sequence of base positions to use in the MCMC, weighted by the per-position weights
    bselect_samples <- c(first_pass_seq,
                         sample(1:B,
                                prob=pos_weights,
                                size=mcmc_total_steps,
                                replace=TRUE))

    unif_samples <- runif(mcmc_total_steps)
    beta_vec <- rep(0, B)
    beta_l1_samples <- rep(0, mcmc_num_samples)
    sample_ctr = 1;
    num_tries <- 0
    num_accept <- 0
    
    for (si in 1:mcmc_total_steps) {
##        if (si %% 3000 == 0) {
##            print(sprintf("  have done %d steps", si))
##        }
        b <- bselect_samples[si]
        cur_beta_b <- beta_vec[b]

        if (si <= B) {
            is_initializing_step <- TRUE
        }
        else {
            is_initializing_step <- FALSE
        }

        new_beta_b <- g_pick_outcome_mcmc_step_given_b(cmat,
                                                    cmat_sum,
                                                    sigmamat,
                                                    rseq,
                                                    fdtot,
                                                    total_length,
                                                    beta_vec,
                                                    beta_l1,
                                                    b,
                                                    unif_samples[si],
                                                    B,
                                                    is_initializing_step,
                                                    pf_compute_log_prob_ratio_change_num_binding_sites)

        sigmamat_new <- sigmamat
        
        # if new_beta_b not same as cur_beta_b, update sigmamat
        if (new_beta_b != cur_beta_b) {
            
            if (b <= B - l + 1) {
                seq_fwd <- rseq[b:(b+l-1)]
            }
            else {
                seq_fwd <- NULL
            }
            
            if (b >= l) {
                seq_rev <- 5-rseq[b:(b-l+1)]
            }
            else {
                seq_rev <- NULL
            }
            
            if (new_beta_b != 0 && cur_beta_b == 0) {
                if (new_beta_b > 0) {
                    seq_use <- seq_fwd
                }
                else {
                    seq_use <- seq_rev
                }
                                        # ----------------- USED ONLY IN DEBUGGING ------------------
                                        #   cat("   new binding site: ")  # DEBUG
                                        #    cat(seq_use) # DEBUG
                                        #    cat(sprintf("  PWM score: %f", pwm.score(cmat, seq_use) - log.l.seqs.bg))
                                        #    cat("\n")    #DEBUG
                                        # ----------------- USED ONLY IN DEBUGGING ------------------
                
                                        # update the sigma matrix to reflect any changes in the number of binding sites
                sigmamat_new[cbind(1:l, seq_use)] <- sigmamat_new[cbind(1:l, seq_use)] + 1
            }
            else {
                if (new_beta_b == 0 && cur_beta_b != 0) {
                                        # we are getting rid of a binding sites
                    if (cur_beta_b > 0) {
                        seq_use <- seq_fwd
                    }
                    else   {
                        seq_use <- seq_rev
                    }
                    sigmamat_new[cbind(1:l, seq_use)] <- sigmamat_new[cbind(1:l, seq_use)] - 1    
                                        # ----------------- USED ONLY IN DEBUGGING ------------------
                                        #      stopifnot(all(sigmamat_new >= 0)) # DEBUG
                                        # ----------------- USED ONLY IN DEBUGGING ------------------
                }
                else {
                    if (new_beta_b != 0 && cur_beta_b != 0 && cur_beta_b != new_beta_b) {
                        if (new_beta_b > 0) {
                            seq_use_new <- seq_fwd
                        }
                    else {
                        seq_use_new <- seq_rev
                    }
                        if (cur_beta_b > 0) {
                            seq_use_old <- seq_fwd
                        }
                        else {
                            seq_use_old <- seq_rev
                        }
                    
                        sigmamat_new[cbind(1:l, seq_use_new)] <- sigmamat_new[cbind(1:l, seq_use_new)] + 1
                        sigmamat_new[cbind(1:l, seq_use_old)] <- sigmamat_new[cbind(1:l, seq_use_old)] - 1
                                        # ----------------- USED ONLY IN DEBUGGING ------------------
                                        #        stopifnot(all(sigmamat_new >= 0)) # DEBUG
                                        # ----------------- USED ONLY IN DEBUGGING ------------------
                    }
                }
            }

            if (new_beta_b != 0) {
                if (cur_beta_b == 0) {
                    beta_l1 <- beta_l1 + 1
                }
            }
            else {
                if (cur_beta_b != 0) {
                    beta_l1 <- beta_l1 - 1
                }
            }
            
            beta_vec[b] <- new_beta_b
            sigmamat <- sigmamat_new
        }

        if (si > mcmc_burn_in_steps + B) {

            num_tries <- num_tries + 1
            if (new_beta_b != cur_beta_b) {
                num_accept <- num_accept + 1
            }
            
            if ((si - mcmc_burn_in_steps - B) %% mcmc_steps_per_sample == 0) {
                beta_l1_samples[sample_ctr] <- beta_l1
                sample_ctr <- sample_ctr + 1
            }
        }
    }

    accept_ratio <- num_accept/num_tries

    list(samples=beta_l1_samples,
         accept=accept_ratio)
}

g_compute_log_prob_ratio_change_num_binding_sites_geometric_simple <- function(beta_l1, Lr) {
    -log(2) + log((beta_l1 + 1)/(Lr - beta_l1))
}

g_sample_posterior_beta_using_mcmc_driver <- function(cmat,
                                                      rseq,
                                                      fdtot,
                                                      mcmc_weight_exponent,
                                                      pf_compute_log_prob_ratio_change_num_binding_sites) {
    g_sample_posterior_beta_using_mcmc(cmat,
                                     rseq,
                                     fdtot,
                                     g_mcmc_weight_exponent,
                                     g_mcmc_burn_in_steps,
                                     g_mcmc_steps_per_sample,
                                     g_mcmc_num_samples,
                                     pf_compute_log_prob_ratio_change_num_binding_sites)
}

g_add_binding_site <- function(rseq, site_seq_list, site_ctr) {
    l <- length(site_seq_list[[1]])
    pos <- l*(site_ctr - 1) + 1
    if (site_ctr > 1) {
        rseq <- g_add_binding_site(rseq,
                                 site_seq_list,
                                 site_ctr - 1)
    }
    rseq[pos:(pos + l - 1 )] <- site_seq_list[[site_ctr]]
    rseq
}

g_compare_models_fixed_num_sites <- function(cmat, rseq_base, nsites) {
    l <- nrow(cmat)

    binding_site_seqs <- lapply(1:max(nsites), function(i) {
        apply(cmat, 1, function(cmat_row) {
                                  cmat_row_sum <- sum(cmat_row)
                                  sample(1:4, prob=cmat_row/cmat_row_sum, size=1)
                              }
              )})

    total_length <- sum(g_fdtot) + length(rseq_base)

    ret_list <- lapply(nsites, function(nsite) {
        rseq <- g_add_binding_site(rseq_base,
                                   binding_site_seqs,
                                   nsite)
        pf_compute_log_prob_ratio_change_num_binding_sites_tfbsincbeta_simple <- function(beta_l1, Lr) {
            g_compute_log_prob_ratio_change_num_binding_sites_tfbsincbeta(l,
                                                                     beta_l1,
                                                                     g_lambda_max,
                                                                     g_alpha,
                                                                     g_nu,
                                                                     Lr) }
        
        samples_tfbsincbeta_list <- g_sample_posterior_beta_using_mcmc(cmat,
                                                             rseq,
                                                             g_fdtot,
                                                             g_mcmc_weight_exponent,
                                                             g_mcmc_burn_in_steps,
                                                             g_mcmc_steps_per_sample,
                                                             g_mcmc_num_samples,
                                                             pf_compute_log_prob_ratio_change_num_binding_sites_tfbsincbeta_simple)

        samples_geometric_list <- g_sample_posterior_beta_using_mcmc(cmat,
                                                                 rseq,
                                                                 g_fdtot,
                                                                 g_mcmc_weight_exponent,
                                                                 g_mcmc_burn_in_steps,
                                                                 g_mcmc_steps_per_sample,
                                                                 g_mcmc_num_samples,
                                                                 g_compute_log_prob_ratio_change_num_binding_sites_geometric_simple)

        samples_tfbsincbeta <- samples_tfbsincbeta_list$samples
        samples_geometric <- samples_geometric_list$samples
        
        retmat <- cbind(samples_tfbsincbeta,
                        samples_geometric,
                        rep(nsite, length(samples_tfbsincbeta)),
                        rep(samples_tfbsincbeta_list$accept, length(samples_tfbsincbeta)),
                        rep(samples_geometric_list$accept, length(samples_tfbsincbeta)))


        colnames(retmat) <- c("tfbsincbeta", "geometric", "nsite", "accept_tfbsincbeta", "accept_geometric")

        retmat
    })

    ret_mat <- do.call(rbind, ret_list)

    ret_mat
}

g_compare_models_fixed_num_sites_multtfs <- function(cmat_list, rseq_base, nsites) {
    ret_list <- lapply(cmat_list,
                       function(mycmat) {
                           g_compare_models_fixed_num_sites(mycmat, rseq_base, nsites)
                       })
    ret_list_labeled <- lapply(names(ret_list),
                               function(tfbs_name) {
                                   ret_list[[tfbs_name]] <- data.frame(ret_list[[tfbs_name]],
                                                                       tfbs=tfbs_name,
                                                                       stringsAsFactors=FALSE)
                               })
    
    ret_data <- do.call(rbind, ret_list_labeled)

    ret_data
}


g_compare_models_fixed_num_sites_multtfs_rep <- function(cmat_list, nsites, N) {
    rseq_base <- sample(1:4,
                          prob=g_fdtot/g_fdtot_sum,
                          size=g_rseq_len,
                          replace=TRUE)
    ret_data <- do.call(rbind, mclapply(1:N,
                                        function(myarg) {
                                            g_compare_models_fixed_num_sites_multtfs(cmat_list, rseq_base, nsites)
                                        },
                                        mc.cores=60))                                       
    ret_data
}    

print(system.time(
    ret_data <- g_compare_models_fixed_num_sites_multtfs_rep(g_tfbs_pcms_try_list,
                                                             g_num_sites,
                                                             g_num_trials)))

write.table(ret_data,
            file="results.txt",
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=TRUE)

                                                         


