rm(list=ls())

default_path <- "D:/Files/Courses/0 SL/HDSL_code" 
setwd(default_path)

source("HDSL.R")
source("simulation_setting.R")

####################################################
## generate simulation data

seed <- 2026
scenario <- 1
p <- 50
noise_sd = 0.3

dat <- simulatedata(seed, scenario = scenario, p = p, noise_sd=noise_sd)
X <- dat$X
y <- dat$y
n <- nrow(X)
p <- ncol(X)
U <- dat$U
V <- dat$V
B <- dat$B


####################################################
## estimation and evaluation

## HDSL 
rini <- 12
lamU <- 0.2
lamV <- 0.2
fmmM <- NA
alpha <- 0.01
method <- 'HDSL'
fit <- algorithm2(X = X, y = y, method = method, r = rini, 
                  lamU = lamU, lamV = lamV, fmmM = fmmM, alpha = alpha)
fit <- modelBIC(fit)
evaluate_result <- evaluate(fit = fit, U = U, V = V, B = B)
evaluate_print <- c(signif(unlist(evaluate_result), 4))
print(evaluate_print)


## HDSLV
rini <- 12
lamU <- 0.2
lamV <- NA
fmmM <- NA
alpha <- 0.01
method <- "HDSLV"
fit <- algorithm2(X = X, y = y, method = method, r = rini, 
                  lamU = lamU, lamV = lamV, fmmM = fmmM, alpha = alpha)
fit <- modelBIC(fit)
evaluate_result <- evaluate(fit = fit, U = U, V = V, B = B)
evaluate_print <- c(signif(unlist(evaluate_result), 4))
print(evaluate_print)


## FMMenet
rini <- NA
lamU <- NA
lamV <- NA
fmmM <- 3
alpha <- NA
method <- "FMMenet"
fit <- algorithm2(X = X, y = y, method = method, r = rini, 
                  lamU = lamU, lamV = lamV, fmmM = fmmM, alpha = alpha)
fit <- modelBIC(fit)
evaluate_result <- evaluate(fit = fit, U = U, V = V, B = B)
evaluate_print <- c(signif(unlist(evaluate_result), 4))
print(evaluate_print)



