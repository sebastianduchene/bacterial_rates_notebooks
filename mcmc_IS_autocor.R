library(NELSI)
options(width = 150)

likelihood_function <- function(tree, is_observed, par){
    
    get_variances <- function(path_times, path_is, alpha, beta){
      # Note that times should start from zero and increase
        single_variances <- vector()
        for(i in 2:length(path_times)){
           var_time <- (path_times[i] - path_times[i-1]) * alpha
           var_is <- abs(path_is[i-1] * beta)
           single_variances[i-1] <- var_time + var_is
        }
        return(single_variances)
    }
    
    alpha <- par[1]
    beta <- par[2]
    times <- allnode.times(tree)
    times <- abs(times - max(times))
    likelihoods <- vector()
    nodes_count <- vector()
    # Each iteration is a path
    for(tip in 1:length(tree$tip.label)){
        path <- get.ancestor.nodes.branches(tree, tip)$ancestor.nodes
        path_times <- times[names(times) %in% path]
        path_times <- c(path_times[2:length(path_times)], path_times[1])
        path_is <- is_observed[names(is_observed) %in% path]
        path_is <- c(path_is[2:length(path_is)], path_is[1])
        ## And only for nodes for which the likelihood has not been calculated. 
        path_times <- path_times[!(names(path_times) %in% nodes_count)]
        path_is <- path_is[!(names(path_is) %in% nodes_count)]
        nodes_count <- c(nodes_count, names(path_times))
        
        variances <- get_variances(path_times, path_is, alpha, beta)
        single_likelihoods <- log(dnorm(path_is[-1], 0, variances))
        likelihoods[tip] <- sum(single_likelihoods)
    }
    sum_likelihoods <- sum(likelihoods)
    if(sum_likelihoods == Inf | sum_likelihoods == -Inf) stop('Error. Likelihood overflow or underflow')
    return(sum_likelihoods)
}

prior_function <- function(par){
    if(any(par < 0)) stop('Error computing prior: alpha and beta cannot be <0')
    prior_alpha <- log(dunif(par[1], 0, 100))
    prior_beta <- log(dunif(par[2], 0, 100))
    return(prior_alpha+prior_beta)
}


run_mcmc <- function(tree, is_observed, start_par, n_steps){
    start_likelihood <- likelihood_function(tree, is_observed, start_par)
    start_prior <- prior_function(start_par)
    start_posterior <- start_likelihood + start_prior

    output_log <- matrix(NA, n_steps, 6)
    colnames(output_log) <- c('step', 'likelihood', 'prior', 'posterior', 'alpha', 'beta')
    output_log[1, ] <- c(1, start_likelihood, start_prior, start_posterior, start_par)
    
    proposal_function <- function(par){
        alpha_proposed <- abs(par[1] + rnorm(1, 0, 0.01))
        beta_propsed <- abs(par[2] + rnorm(1, 0, 0.01))
        return(c(alpha_proposed, beta_propsed))
    }
    
    for(i in 2:n_steps){
        proposed_par <- proposal_function(output_log[i-1, c(5, 6)])
        proposed_prior <- prior_function(proposed_par)
        proposed_likelihood <- likelihood_function(tree, is_observed, proposed_par)
        proposed_posterior <- proposed_prior + proposed_likelihood    
        if(proposed_posterior > output_log[i-1, 4]){
            output_log[i, ] <- c(i, proposed_likelihood, proposed_prior,
                                  proposed_posterior, proposed_par)
            }else{
                mh_ratio <- exp(proposed_posterior - output_log[i-1, 4])
                if(mh_ratio > runif(1)){
                    output_log[i, ] <- c(i, proposed_likelihood, proposed_prior,
                              proposed_posterior, proposed_par)
                    }else{
                        output_log[i, ] <- output_log[i-1, ]
                        output_log[i, 1] <- i
                        }
                }
        if(i%%10 == 0){
	   out <- output_log[i, ]
	   names(out) <- NULL
           print(out)
        }
    }
    
    return(output_log)
}


###################
tr <- read.tree('tree_states.tree')
states <- read.table('tree_states.csv', head = T, sep = ',')

is_observed <- states[, 2]
names(is_observed) <- c(1:199)

m1 <- run_mcmc(tree = tr, is_observed = is_observed, start_par = c(1, 1), n_steps = 10000)

pdf('out.pdf')
par(mfrow = c(2, 2))
plot(m1[200:nrow(m1), 5], type = 'l', col = rgb(1, 0, 0, 0.5))
hist(m1[200:nrow(m1), 5], col = rgb(1, 0, 0, 0.5))
plot(m1[, 6], type = 'l', col = rgb(1, 0, 0, 0.5))
hist(m1[200:nrow(m1), 6], col = rgb(1, 0, 0, 0.5))
dev.off()