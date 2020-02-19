# Set seed
set.seed(14)

# Require package
library(pROC)

# Read data
dat <- read.csv("wu.csv", header=FALSE)
dat$group <- c(rep("resistant", 21), rep("responder", 50))
dat <- dat[, c(2, 3)]
names(dat)[1] <- "delta_hamd"

# Impute missing data point
# One is missing in the reponder group. We can impute it because we know the mean, which is equal to point "Bar45".
mean(dat$delta_hamd[dat$group == "responder"])
# This is very close to the plotted value of 0.799. Imputation not possible on the basis of guessing an overplotted point. The last point is probably missing from the graph.
# We attempt nonetheless to impute one more at the mean
dat_imp <- dat
dat_imp <- rbind(dat_imp, c(8, "responder"))
dat_imp$delta_hamd <- as.numeric(dat_imp$delta_hamd)

# Estimate ROC
roc <- roc(dat$group, dat$delta_hamd)
plot(roc)
ci(roc)
roc

# Estimate ROC with imputation
roc_imp <- roc(dat_imp$group, dat_imp$delta_hamd)
plot(roc_imp)
ci(roc_imp)
roc_imp

# Note: the following code is deprecated - permutations not necessary
# Bootstrap samples
n_iterations <- 1000
rocs <- vector("list", n_iterations)
roc_tests <- vector("list", n_iterations)
sensitivities <- vector("list", n_iterations)
specificities <- vector("list", n_iterations)

for(i in 1:n_iterations){
# Randomise order
dat$group_random <- c(rep("resistant", 21), rep("responder", 50))
rows <- sample(nrow(dat))
dat$group_random <- dat$group_random[rows]

roc_boot <- roc(dat$group_random, dat$delta_hamd)
rocs[[i]] <- roc_boot
sensitivities[[i]] <- roc_boot$sensitivities
specificities[[i]] <- roc_boot$specificities

roc_tests[[i]] <- roc.test(roc, roc_boot)
}

plot(rocs[[1]], col = "pink", lwd = 0.5)
for(i in 2:n_iterations){
  lines(rocs[[i]], col = "pink", lwd = 0.5)
}
lines(roc)  
legend("bottomright", lwd = c(2, 0.5, 2), col = c("black", "pink", "red"), 
       legend = c("original", "bootstrapped", "average of bootstrapped"))

sens <- matrix(unlist(sensitivities), nrow = 69)
spec <- matrix(unlist(specificities), nrow = 69)
sens_avg <- rowMeans(sens)
spec_avg <- rowMeans(spec)
lines(sens_avg ~ spec_avg, col = "red", lwd = 2)

p <- vector()
for(i in 1:n_iterations){
  pval <- roc_tests[[i]]$p.value
  p <- c(p, pval)
}
sum(p < 0.05)

auroc <- vector()
for(i in 1:n_iterations){
  auroc_estimate <- rocs[[i]]$auc
  auroc <- c(auroc, auroc_estimate)
}
auroc_mean <- mean(auroc)
error <- qnorm(0.975)*auroc_mean/sqrt(length(auroc))
left <- auroc_mean - error
right <- auroc_mean + error

