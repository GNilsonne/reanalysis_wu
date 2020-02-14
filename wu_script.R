# Set seed
set.seed(14)

# Require package
library(pROC)

# Read data
dat <- read.csv("C:/Users/gusta/Box Sync/Gustavs_arbete/Pek/Commentary_Wu/wu.csv", header=FALSE)
dat$group <- c(rep("resistant", 21), rep("responder", 50))
dat <- dat[, c(2, 3)]
names(dat)[1] <- "delta_hamd"

# Estimate ROC
roc <- roc(dat$group, dat$delta_hamd)
plot(roc)
ci(roc)

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

