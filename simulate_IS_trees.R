

library(NELSI)

simulate_IS <- function(tree, root_state, par, basic_variance){
    alpha <- par[1] # time variance weight 
    beta <- par[2] # IS variance weight
    gamma <- par[3] # base variance weight. Usually set to 1 for a stochastic process
    node_ages <- allnode.times(tree)
    node_ages <- abs(node_ages - max(node_ages))
    root_node <- unique(tree$edge[!(tr$edge[, 1] %in% tree$edge[, 2]), 1])
    node_states <- cbind(as.numeric(names(node_ages)), rep(NA, length(node_ages)))
    node_states[node_states[, 1] == root_node, 2] <- root_state
    all_tips <- 1:length(tree$tip.label)

    tips_traces <- list()# This list stores some useful information for plotting IS along each path in the tree
    for(tip in all_tips){
        tip_temp <- get.ancestor.nodes.branches(tree, tip)
        nodes_i <- rev(tip_temp$ancestor.nodes)
        tips_traces[[tip]] <- list()
        tips_traces[[tip]][[1]] <- vector()#differences in node ages
        tips_traces[[tip]][[2]] <- vector()#IS number
        tips_traces[[tip]][[3]] <- vector()#differences in IS number
        for(i in 2:length(nodes_i)){
            if(is.na(node_states[node_states[, 1] == nodes_i[i], 2])){
                diff_times <- abs(node_ages[names(node_ages) == nodes_i[i]] - node_ages[names(node_ages) == nodes_i[i-1]])
                nu_time <- diff_times * alpha
                nu_state <- abs(node_states[node_states[, 1] == nodes_i[i-1], 2]) * beta
		nu_base <- basic_variance * par[3]
                nu_total <- nu_time + nu_state + nu_base
                new_state <- node_states[node_states[, 1] == nodes_i[i-1], 2] + (rnorm(1, 0, nu_total))
                node_states[node_states[, 1] == nodes_i[i], 2] <- new_state
                tips_traces[[tip]][[1]] <- c(tips_traces[[tip]][[1]], diff_times)
                tips_traces[[tip]][[2]] <- c(tips_traces[[tip]][[2]], new_state)
                tips_traces[[tip]][[3]] <- c(tips_traces[[tip]][[3]],
                                new_state - node_states[node_states[, 1] == nodes_i[i-1], 2])
            }
            nodes_in_path <- node_states[, 1] %in% nodes_i
            tips_traces[[tip]][[4]] <- cbind(node_ages[nodes_in_path], node_states[nodes_in_path, 2])
        }
    }

    return(node_states)
}



tr <- rtree(200)

sim_1 <- simulate_IS(tr, root_state = 200, par = c(1, 100, 1), 0.1)
sim_2 <- simulate_IS(tr, root_state = 200, par = c(1, 0.1, 1), 0.1)

par(mfrow = c(2, 2))
plot(sim_1[, 1], sim_1[, 2])
plot(tr, show.tip.label = F)
nodelabels(frame = 'circle')

plot(sim_2[, 1], sim_2[, 2])
plot(tr, show.tip.label = F)
