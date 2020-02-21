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

# Estimate ROC
roc <- roc(dat$group, dat$delta_hamd)
plot(roc)
ci(roc)
roc