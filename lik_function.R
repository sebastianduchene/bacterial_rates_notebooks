
library(NELSI)
options(width = 120)

likelihood_function <- function(tree, is_observed, par, basic_variance){
    get_variances <- function(path_times, path_is, alpha, beta, gamma, basic_variance){
      # Note that times should start from zero and increase
        single_variances <- vector()
        for(i in 2:length(path_times)){
            var_time <- path_is[i-1]* (alpha^(path_times[i] - path_times[i-1]))
            var_is <- abs((path_is[i-1] / max(path_is)) * beta)
	          var_base <- basic_variance * gamma
           single_variances[i-1] <- var_time + var_is + var_base
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
        variances <- get_variances(path_times, path_is, alpha, beta, gamma, basic_variance)
        single_likelihoods <- log(dnorm(path_is[-1], 0, variances))
        likelihoods[tip] <- sum(single_likelihoods)

    }
    sum_likelihoods <- sum(likelihoods)
    if(sum_likelihoods == Inf | sum_likelihoods == -Inf) stop('Error. Likelihood overflow or underflow')
    return(sum_likelihoods)
}


# Example of using likelihood function (uncomment to run)
########################################
#tree <- rtree(100)
#nodetimes <- allnode.times(tr)
#is_observed <- abs(rnorm(length(nodetimes), 0, 3)) # randomly evolving IS numbers through time. Note that the numebrs are nod discrete
#names(is_observed) <- 1:length(is_observed)

#par <- c(3, 2, 2)
#basic_variance <- 0.1
#likelihood_function(tree, is_observed, par, basic_variance)
