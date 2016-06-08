library(NELSI)
source('mcmc_IS_autocor.R')

tr <- read.tree('tree_05_10.tree')
states <- read.table('par_05_10.csv', head = T, sep = ',')

is_observed <- states[, 2]
names(is_observed) <- states[, 1]

m1 <- run_mcmc(tree = tr, is_observed = is_observed, start_par = c(0.5, 5), n_steps = 50000)

pdf('out_par_05_10.pdf')
par(mfrow = c(2, 2))
plot(m1[200:nrow(m1), 5], type = 'l', col = rgb(1, 0, 0, 0.5))
hist(m1[200:nrow(m1), 5], col = rgb(1, 0, 0, 0.5))
plot(m1[, 6], type = 'l', col = rgb(1, 0, 0, 0.5))
hist(m1[200:nrow(m1), 6], col = rgb(1, 0, 0, 0.5))
dev.off()
