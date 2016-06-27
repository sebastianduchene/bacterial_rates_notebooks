library(NELSI)
options(width = 150)

likelihood_function <- function(tree, is_observed, par, basic_variance){

    get_variances <- function(path_times, path_is, alpha, beta, gamma, basic_variance){
      # Note that times should start from zero and increase
        single_variances <- vector()
        for(i in 2:length(path_times)){
           var_time <- (path_times[i] - path_times[i-1]) * alpha
           var_is <- abs(path_is[i-1] * beta)
	   var_base <- basic_variance * gamma
           if(var_is > var_time){
              single_variances[i-1] <- var_base
           }else{
           single_variances[i-1] <- (var_time - var_is) + var_base
       }
        }
        return(single_variances)
    }

    alpha <- par[1]
    beta <- par[2]
    gamma <- par[3]
    times <- allnode.times(tree)
    times <- abs(times - max(times))
    likelihoods <- vector()
    nodes_count <- vector()
    # Each iteration is a path
    for(tip in 1:length(tree$tip.label) ){
        path <- get.ancestor.nodes.branches(tree, tip)$ancestor.nodes
        path_times <- times[names(times) %in% path]
        path_times <- c(path_times[2:length(path_times)], path_times[1])
        path_is <- is_observed[names(is_observed) %in% path]
        path_is <- c(path_is[2:length(path_is)], path_is[1])
        ## And only for nodes for which the likelihood has not been calculated.
#        path_times <- path_times[!(names(path_times) %in% nodes_count)]
#        path_is <- path_is[!(names(path_is) %in% nodes_count)]
#        nodes_count <- c(nodes_count, names(path_times))
        variances <- get_variances(path_times, path_is, alpha, beta, gamma, basic_variance)
        single_likelihoods <- log(dnorm(path_is[-1], 0, variances))
        likelihoods[tip] <- sum(single_likelihoods)
    }
    #repeat for final tip/path
#    path <- get.ancestor.nodes.branches(tree, length(tree$tip.label))$ancestor.nodes
 #   path_times <- times[names(times) %in% path]
  #  path_times <- c(path_times[2:length(path_times)], path_times[1])
#    path_is <- is_observed[names(is_observed) %in% path]
#    path_is <- c(path_is[2:length(path_is)], path_is[1])
#    variances <- get_variances(path_times, path_is, alpha, beta, gamma, basic_variance)
#    single_likelihoods <- log(dnorm(path_is[-1], 0, variances))
#    likelihoods[length(tree$tip.label)] <- sum(single_likelihoods)

    sum_likelihoods <- sum(likelihoods)
#    if(sum_likelihoods == Inf | sum_likelihoods == -Inf) stop('Error. Likelihood overflow or underflow')
    return(sum_likelihoods)
}

prior_function <- function(par){
    if(any(par < 0)) stop('Error computing prior: alpha and beta cannot be <0')
    prior_alpha <- log(dunif(par[1], 0, 1000))
    prior_beta <- log(dunif(par[2], 0, 1000))
    prior_gamma <- log(dnorm(par[3], 0, 0.1))
    return(prior_alpha+prior_beta+prior_gamma)
}


run_mcmc <- function(tree, is_observed, start_par, n_steps, basic_variance){
    start_likelihood <- likelihood_function(tree, is_observed, start_par, basic_variance)
    start_prior <- prior_function(start_par)
    start_posterior <- start_likelihood + start_prior

    output_log <- matrix(NA, n_steps, 7)
    colnames(output_log) <- c('step', 'likelihood', 'prior', 'posterior', 'alpha', 'beta', 'gamma')
    output_log[1, ] <- c(1, start_likelihood, start_prior, start_posterior, start_par)

    proposal_function <- function(par){
        alpha_proposed <- abs(par[1] + rnorm(1, 0, 0.1))
        beta_propsed <- abs(par[2] + rnorm(1, 0, 0.1))
	gamma_proposed <- abs(par[3] + rnorm(1, 0, 0.1))
        return(c(alpha_proposed, beta_propsed, gamma_proposed))
    }

    for(i in 2:n_steps){
        proposed_par <- proposal_function(output_log[i-1, c(5:7)])
        proposed_prior <- prior_function(proposed_par)
        proposed_likelihood <- likelihood_function(tree, is_observed, proposed_par, basic_variance)
        attempts <- 0
        if(proposed_likelihood == -Inf){
            print('WARNING: Likelihood underflow')
            while(attempts < 100 && proposed_likelihood == -Inf){
                proposed_par <- proposal_function(output_log[i-1, c(5:7)])
                proposed_likelihood <- likelihood_function(tree, is_observed, proposed_par, basic_variance)
               attempts <- attempts +1
            }
            if(attempts == 100){
                stop('SORRY. Likelihood is still -Inf after 100 attempts')
            }
            print(paste('Success underflow solved after', attempts, 'attempts at step', i))
           proposed_prior <- prior_function(proposed_par)
        }
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
#tr <- read.tree('tree_states_alpha.tree')
#states <- read.table('tree_states_alpha.csv', head = T, sep = ',')

#is_observed <- states[, 2]
#	    is_observed <- rnorm(length(is_observed), 0, 0.1)
#names(is_observed) <- c(1:length(is_observed))



#m1 <- run_mcmc(tree = tr, is_observed = is_observed, start_par = c(1, 1, 1), n_steps = 500, basic_variance = .1)
#pdf('out_alpha.pdf')
#par(mfrow = c(3, 2))
#plot(m1[100:nrow(m1), 5])#, type = 'l', col = rgb(1, 0, 0, 0.5))
#hist(m1[100:nrow(m1), 5])#, col = rgb(1, 0, 0, 0.5))
#plot(m1[100:nrow(m1), 6])#, type = 'l', col = rgb(1, 0, 0, 0.5))
#hist(m1[100:nrow(m1), 6])#, col = rgb(1, 0, 0, 0.5))
#plot(m1[100:nrow(m1), 7])#, type = 'l', col = rgb(1, 0, 0, 0.5))
#hist(m1[100:nrow(m1), 7])#, col = rgb(1, 0, 0, 0.5))
#dev.off()
