library(sparsepca) # spca
library(ncvreg) # mcp regression
library(Rcpp)  
library(lsa)
library(mclust)
library(data.table)
library(flexmix)
library(MASS)
library(nnet)
library(class)
library(glmnet)
sourceCpp("new_admm.cpp")

############################################################
#####  Estimation  ######################################### 
############################################################

norm21 <- function(V) {
  V <- as.matrix(V) # in case V is a vector
  rr <- apply(V, 2, function(x){!all(x == 0)})
  V[, rr] <- apply(as.matrix(V[, rr]), 2, function(x) {x / sqrt(sum(x^2))})  
  return(V)
}

myplot <- function(Uhat, U, title = NULL, jit = 1e-4, ...) {
  plot(jitter(Uhat[, 1], amount = jit), 
       jitter(Uhat[, 2], amount = jit), 
       main = title, col = kmeans(U, 3)$cluster + 1, xlim = c(-6, 6), ylim = c(-6, 6), 
       ...)
  points(U[, 1], U[, 2], pch = 1, cex = 3, lwd = 2)
}

softthre <- function(z, lam) {
  if(z > lam) {return(z - lam)
  } else if (z >= -lam) {return(0)
  } else {return(z + lam)}
}

fMCP <- function(z, lam, gam = 3) {
  if(abs(z) <= gam * lam) {
    return(softthre(z, lam) / (1 - 1 / gam))
  } else {return(z)}
}

ncvupdate <- function(y, X, b, lam) {
  n <- length(y)
  residual <- y - X %*% b
  bnew <- b
  for(j in 1:length(b)) {
    zj <- t(X[, j]) %*% residual / n + b[j]
    bnew[j] <- fMCP(zj, lam)
  }
  return(bnew)
}

Vinitialize <- function(Vinit, p, r, X, V = NA, alpha = 0.1) {
  if(Vinit == "SPCA") Vhat <- spca(X, k = r, alpha = alpha, verbose = FALSE)$loadings  
  for(j in 1:ncol(Vhat)) {if(sum(Vhat[, j]) < 0) {Vhat[, j] <- -Vhat[, j]}}
  return(Vhat)
}


Uinitialize <- function(Uinit, n, r, y, X, Vhat) {
  if(Uinit == "FMME") {
    fmmM <- 20
    p <- ncol(X)
    ddd <- data.frame(y = y, X = X %*% Vhat)
    whichnot0 <- !apply(ddd[, -1], 2, function(x) {all(x == 0)})
    ddd <- ddd[, c(TRUE, whichnot0)]
    if(ncol(as.matrix(ddd)) != 1) {
      tempr <- ncol(ddd) - 1
      m1 <- flexmix(y ~ ., data = ddd, k = fmmM)
      Uhatfill <- as.matrix(t(parameters(m1)[, clusters(m1)])[, c(2:(min(r, tempr)+1))])
      Uhat <- matrix(0, n, r)
      Uhat[, whichnot0] <- Uhatfill
    } else {
      Uhat <- matrix(0, n, r)
      Vhat <- matrix(0, p, r)
      Bhat <- matrix(0, n, p)
    }
    Uhat <- Uhat + matrix(rnorm(n * r, sd = 0), ncol = r) # Make sd = 0 but add M (cluster size)
  }
  if(Uinit == "FMMElasso") {
    fmmM <- 20
    p <- ncol(X)
    ddd <- data.frame(y = y, X = X %*% Vhat)
    whichnot0 <- !apply(ddd[, -1], 2, function(x) {all(x == 0)})
    ddd <- ddd[, c(TRUE, whichnot0)]
    if(ncol(as.matrix(ddd)) != 1) {
      tempr <- ncol(ddd) - 1
      temp <- stepFlexmix(y ~ ., data = ddd, k = fmmM, nrep = 3, drop = TRUE,
                          model = FLXMRglmnet(family = "gaussian", alpha = 1, intercept = FALSE),
                          control = list(iter.max = 300))
      Uhatfill <- t(parameters(temp)[, clusters(temp)])[, -c(1, nrow(parameters(temp)))]
      Uhat <- matrix(0, n, r)
      Uhat[, whichnot0] <- Uhatfill
    } else {
      Uhat <- matrix(0, n, r)
      Vhat <- matrix(0, p, r)
      Bhat <- matrix(0, n, p)
    }
    Uhat <- Uhat + matrix(rnorm(n * r, sd = 0), ncol = r) # Make sd = 0 but add M (cluster size)
  }
  return(Uhat)
}

HDSL <- function(X, y, method, r, lamU, lamV = NA, alpha = NA, kap = 0) {
  tol <- 5e-5
  n <- nrow(X)
  p <- ncol(X)
  active_rank <- rep(TRUE, r)
  Vhat <- Vinitialize(Vinit = "SPCA", p = p, r = r, X = X, V = V, alpha = alpha)
  Vhat <- norm21(Vhat) ### normalize V and then calculate U
  if(method %in% c("HDSL")){
    Uhat <- Uinitialize(Uinit = "FMMElasso", n = n, r = r, y = y, X = X, Vhat = Vhat)
    Vhat_init <- Vhat
    Uhat_init <- Uhat
    active_rank <- !apply(Uhat, 2, function(x){all(x == 0)})
    Vhat[, !active_rank] <- 0 # MAKE ZERO
  }
  if(method %in% c("HDSLV")){
    Uhat <- Uinitialize(Uinit = "FMME", n = n, r = r, y = y, X = X, Vhat = Vhat)
    Vhat_init <- Vhat
    Uhat_init <- Uhat
  }
  Bhat <- Uhat %*% t(Vhat)
  iter <- 1
  admm_iter <- 50
  maxit <- 30
  loss <- lossfunc(y, X, Uhat, Vhat, lamU, lamV)
  ######################## PART I: V fixed, just ADMM for U ####################
  if(method %in% c("HDSLV")) {
    zeros <- matrix(0, n, r)
    Eta1 <- subG_ADMM_new(Yin = y, Xin = rep(1, n), Zin = X %*% Vhat,
                          Eta1in = Uhat, ridgein = 0, # if ridgein = 1, Eta1in ignored
                          Utin = zeros, kappain = kap * lamU, # if kappain = 0, Utin ignored, kappa is the ridge term
                          lambdain = lamU, admm_iterin = admm_iter, tolin = tol)$Eta1
    Uhat <- as.matrix(drop(round(Eta1, 4)))
    Bhat <- Uhat %*% t(Vhat)
  }
  ######################## PART II: V and U both updated further ###############
  if(method %in% c("HDSL")) {
    converge <- FALSE
    Bhatold <- Bhat
    Uupdate_count <- 0
    Vupdate_count <- 0
    while((iter <= maxit) & !converge){
      ########### update Vhat by coordesc ##########
      Z1 <- c()
      for(i in 1:n) {Z1 <- rbind(Z1, c(t(outer(X[i, ], Uhat[i, active_rank]))))}
      Vcheck <- ncvupdate(y = y, X = Z1, b = c(t(Vhat[, active_rank])), lam = lamV)
      Vhat_copy <- matrix(0, nrow(Vhat), ncol(Vhat))
      Vhat_copy[, active_rank] <- norm21(t(matrix(Vcheck, nrow = sum(active_rank)))) # do not norm21 V in iter
      loss_Vmove <- lossfunc(y, X, Uhat, Vhat_copy, lamU, lamV)
      ########### update Uhat #############
      zeros <- matrix(0, n, sum(active_rank))
      Eta1 <- subG_ADMM_new(Yin = y, Xin = rep(1, n), Zin = X %*% Vhat[, active_rank],
                            Eta1in = Uhat[, active_rank], ridgein = 0, # if ridgein = 1, Eta1in ignored
                            Utin = zeros, kappain = kap * lamU, # if kappain = 0, Utin ignored
                            lambdain = lamU, admm_iterin = admm_iter, tolin = tol)$Eta1
      Uhat_copy <- matrix(0, nrow(Uhat), ncol(Uhat))
      Uhat_copy[, active_rank] <- as.matrix(drop(round(Eta1, 4)))
      loss_Umove <- lossfunc(y, X, Uhat_copy, Vhat, lamU, lamV)
      ########### maximum improvement #########
      if(iter == 1 | loss_Vmove <= loss_Umove) {
        Vhat <- Vhat_copy
        Vupdate_count <- Vupdate_count + 1
        loss <- loss_Vmove
        active_rank <- active_rank & !apply(Vhat, 2, function(x){all(x == 0)}) # UPDATE RANK
      } else {
        Uhat <- Uhat_copy
        Uupdate_count <- Uupdate_count + 1
        loss <- loss_Umove
        active_rank <- active_rank & !apply(Vhat, 2, function(x){all(x == 0)}) # UPDATE RANK
      }
      ########### check convergence #######
      Bhat <- Uhat %*% t(Vhat)
      Bdiff <- sqrt(sum((Bhat - Bhatold)^2)) / sqrt(sum(Bhatold^2))
      if(is.nan(Bdiff)) {Bdiff <- 0}
      if(iter %% 5 == 0) {print(paste0("************************ Iter = ", iter, 
                                       "; Bdiff = ", round(Bdiff, 8),
                                       "; loss = ", round(loss, 3),
                                       "; Umove = ", Uupdate_count,
                                       "; Vmove = ", Vupdate_count))}
      converge <- Bdiff < tol
      Bhatold <- Bhat
      iter <- iter + 1
    }
  }
  Bhat <- Uhat %*% t(Vhat)
  list(Uhat = Uhat, Vhat = Vhat, Bhat = Bhat, iter = max(1, iter - 1), Vhat_init = Vhat_init, Uhat_init = Uhat_init)
}


FMM <- function(X, y, method, r, fmmM, alpha){
  n <- nrow(X)
  p <- ncol(X)
  Uhat <- NA; Vhat <- NA; Bhat <- NA; Vhat_init <- NA
  if(method %in% c("FMM")) {
    Vhat <- norm21(Vinitialize(Vinit = "SPCA", p = p, r = r, X = X, V = V, alpha = alpha))
    Vhat_init <- Vhat
    ddd <- data.frame(y = y, X = X %*% Vhat)
    ddd <- ddd[, !apply(ddd, 2, function(x) {all(x == 0)})]
    if(ncol(as.matrix(ddd)) != 1) {
      tempr <- ncol(ddd) - 1
      m1 <- flexmix(y ~ ., data = ddd, k = fmmM)
      Uhat <- as.matrix(t(parameters(m1)[, clusters(m1)])[, c(2:(min(r, tempr)+1))])
      Uhat <- as.matrix(cbind(Uhat, matrix(0, n, r - tempr)))
      Bhat <- Uhat %*% t(Vhat)  
    } else {
      Uhat <- matrix(0, n, r)
      Vhat <- matrix(0, p, r)
      Bhat <- matrix(0, n, p)
    }
  }
  if(method %in% c("FMMlasso")) {
    temp <- stepFlexmix(y ~ X, k = fmmM, nrep = 5, drop = TRUE,
                        model = FLXMRglmnet(family = "gaussian", alpha = 1, 
                                            intercept = FALSE),
                        control = list(iter.max = 300))
    Bhat <- t(parameters(temp)[, clusters(temp)])[, -c(1, nrow(parameters(temp)))]
  }
  if(method %in% c("FMMenet")) {
    temp <- stepFlexmix(y ~ X, k = fmmM, nrep = 5, drop = TRUE,
                        model = FLXMRglmnet(family = "gaussian", alpha = 0.5, 
                                            intercept = FALSE),
                        control = list(iter.max = 300))
    Bhat <- t(parameters(temp)[, clusters(temp)])[, -c(1, nrow(parameters(temp)))]
  }
  list(Uhat = Uhat, Vhat = Vhat, Bhat = Bhat, iter = NA, Vhat_init = Vhat_init, Uhat_init = NA)
}


# algorithm2
# subG_ADMM_new(): MCP fusion, non-sparse, Xin is Z1 in paper eq (2)
algorithm2 <- function(X, y, method, r, lamU = NA, lamV = NA, fmmM = NA, alpha = NA, kap = 0) {
  n <- nrow(X)
  p <- ncol(X)
  if(method %in% c("HDSL", "HDSLV")) {
    ifit <- HDSL(X = X, y = y, method = method, r = r, lamU = lamU, lamV = lamV, alpha = alpha, kap = kap)
  }
  if(method %in% c("FMMlasso", "FMMenet")) {ifit <- FMM(X = X, y = y, method = method, r = r, fmmM = fmmM, alpha = alpha)}
  Uhat = ifit$Uhat
  Vhat = ifit$Vhat
  Bhat = ifit$Bhat 
  iter = ifit$iter
  Uhat_init = ifit$Uhat_init
  Vhat_init = ifit$Vhat_init
  yhat <- rowSums(X * Bhat)
  RSS <- sum((y - yhat)^2)
  neg2lkhd <- n * log(RSS / n)
  Rsq <- 1 - RSS / sum((y - mean(y))^2)
  list(method = method, X = X, y = y, n = n, p = p,
       r = r, lamU = lamU, lamV = lamV, fmmM = fmmM, alpha = alpha, kap = kap, 
       Uhat = Uhat, Vhat = Vhat, Bhat = Bhat, yhat = yhat, iter = iter, Uhat_init = Uhat_init, Vhat_init = Vhat_init,
       RSS = RSS, neg2lkhd = neg2lkhd, Rsq = Rsq)
}



############################################################
#####  Evaluation  ######################################### 
############################################################

modelBIC <- function(fit) {
  method <- fit$method
  n <- fit$n
  p <- fit$p
  r <- fit$r
  Vhat <- fit$Vhat
  Uhat <- fit$Uhat
  Bhat <- fit$Bhat
  neg2lkhd <- fit$neg2lkhd
  if(method %in% c("HDSL", "HDSLV")) {
    rV <- sum(apply(Vhat, 2, function(x){!all(x == 0)}))
    DFVeach <- apply(Vhat, 2, function(x) {sum(unique(x) != 0)})
    DFV <- sum(DFVeach)
    Khat <- nrow(unique(round(Uhat, 1)))
    DFU <- Khat * rV
    DF <- DFU + DFV
    BIC <- neg2lkhd + 0.1 * log(log(n * rV + p * rV)) * log(n) * DFU + 1 * log(log(n * rV + p * rV)) * log(n) * DFV
  }
  if(method %in% c("FMMlasso", "FMMenet")) {
    DFBeach <- apply(unique(Bhat), 1, function(x) {sum(unique(x) != 0)})
    DF <-  sum(DFBeach)
    BIC <- neg2lkhd + log(n) * DF
  }
  fit$rV <- ifelse(exists("rV"), rV, NA)
  fit$Khat <- ifelse(exists("Khat"), Khat, NA)
  fit$DFU <- ifelse(exists("DFU"), DFU, NA)
  fit$DFVeach <- ifelse(exists("DFVeach"), DFVeach, NA)
  fit$DFV <- ifelse(exists("DFV"), DFV, NA)
  fit$DFBeach <- ifelse(exists("DFBeach"), DFBeach, NA)
  fit$DF <- DF
  fit$BIC <- BIC
  return(fit)
}


evaluate <- function(fit, U, V, B) {
  rtrue <- ncol(V)
  method <- fit$method
  X <- fit$X
  y <- fit$y
  r <- fit$r
  Bhat <- fit$Bhat
  Brmse <- sqrt(mean((Bhat - B)^2))
  Brmse_nz <- sqrt(mean((Bhat[B != 0] - B[B != 0])^2))
  B != 0
  if(method %in% c("HDSL", "HDSLV")) {
    Uhat <- as.matrix(fit$Uhat)
    Vhat <- as.matrix(fit$Vhat)
    if(rtrue < r) {
      V <- cbind(V, matrix(0, nrow(V), r - ncol(V)))
      U <- cbind(U, matrix(0, nrow(U), r - ncol(U)))
    }
    for(j in 1:min(rtrue, r)) {
      if(sum((Vhat[, j] - V[, j])^2) > sum((Vhat[, j] + V[, j])^2)) {
        Vhat[, j] <- -Vhat[, j]; Uhat[, j] <- -Uhat[, j]
      }
    }
    Vrmse <- sapply(1:min(rtrue, r),
                    function(i) {pmin(sqrt(mean((V[, i] - Vhat[, i])^2)), sqrt(mean((V[, i] + Vhat[, i])^2)))}
    )
    Vcos <- sapply(1:min(rtrue, r),
                   function(i) {pmax(cosine(V[, i], Vhat[, i]), cosine(V[, i], -Vhat[, i]))}
    )
    Vsens <- senspec(V[, 1:min(rtrue, r)], Vhat[, 1:min(rtrue, r)])$sensitivity
    Vspec <- senspec(V[, 1:min(rtrue, r)], Vhat[, 1:min(rtrue, r)])$specificity
    Urmse <- sapply(1:min(rtrue, r),
                    function(i) {pmin(sqrt(mean((U[, i] - Uhat[, i])^2)), sqrt(mean((U[, i] + Uhat[, i])^2)))}
    )
    Uari <- adjustedRandIndex(as.numeric(as.factor(apply(U, 1, paste0, collapse = "_"))),
                              as.numeric(as.factor(apply(round(Uhat, 1), 1, paste0, collapse = "_"))))
    Bari <- NA
    if(all(V == 0)) {Vcancor <- NA} else {Vcancor<- min(cancor(V, Vhat)$cor)}
    Vrank <- fit$rV
  }
  if(method %in% c("FMMlasso", "FMMenet")) {
    Bari <- adjustedRandIndex(as.numeric(as.factor(apply(round(B, 1), 1, paste0, collapse = "_"))),
                              as.numeric(as.factor(apply(round(Bhat, 1), 1, paste0, collapse = "_"))))
    Urmse <- rep(NA, min(2))
    Vcos <- rep(NA, min(2))
    Uari <- Bari
    Vsens <- NA 
    Vspec <- NA
    Vcancor <- NA
    Vrank <- NA
  }
  list(r = fit$r, lamU = fit$lamU, lamV = fit$lamV, fmmM = fit$fmmM, alpha = fit$alpha, iter = fit$iter, 
       Brmse = Brmse, 
       Brmse_nz = Brmse_nz,
       Bari = Bari,
       Uari = Uari,
       RSS = fit$RSS, neg2lkhd = fit$neg2lkhd, Rsq = fit$Rsq, 
       Khat = fit$Khat, DFU = fit$DFU, 
       DFV = fit$DFV, 
       DF = fit$DF, 
       BIC = fit$BIC,
       Vcancor = Vcancor,
       Vrank = Vrank
       )  
}

# sensitivity specificity
senspec <- function(Vtrue, Vpred) {
  d <- ncol(Vtrue)
  senslist <- c()
  speclist <- c()
  for(j in 1:d) {
    tb <- table(factor(Vtrue[, j] != 0, levels = c(TRUE, FALSE)), factor(Vpred[, j] != 0, levels = c(TRUE, FALSE)))
    sensitivity <- tb[1, 1] / sum(tb[1, ])
    specificity <- tb[2, 2] / sum(tb[2, ])
    senslist <- c(senslist, sensitivity)
    speclist <- c(speclist, specificity)
  }
  list(sensitivity = senslist, specificity = speclist)
}


lamrho <- function(lam, U) {
  tempx <- c(dist(U))
  gamma <- 3
  sum(as.numeric(abs(tempx) <= gamma * lam) * (lam * abs(tempx) - tempx^2 / (2 * gamma)) +
        as.numeric(abs(tempx) > gamma * lam) * (gamma * lam^2 / 2))
}

lossfunc <- function(y, X, U, V, lamU, lamV){sum((y - rowSums(X * (U %*% t(V))))^2) / 2 + 
  lamrho(lam = lamU, U = U) + lamV * sum(abs(V))}



############################################################
#####  Prediction  ######################################### 
############################################################

predict_HDSL <- function(X_test, X_train, fit, membership) {
  ## predict the outcome for HDSL or HDSLV
  ## membership: "softmax" or "KNN".
  Vhat <- fit$Vhat
  Uhat <- fit$Uhat
  Bhat <- fit$Bhat
  Ucenters <- unique(round(Uhat, 1))
  Khat <- nrow(Ucenters)
  cat(paste0("U has ",Khat, " Clusters"))
  
  ## get labels according to Uhat
  labels <- apply(Uhat, 1, function(row) {
    matches <- apply(Ucenters, 1, function(u_row) all(round(row, 1) == u_row))
    which(matches)[1]
  })
  labels <- factor(labels, levels = 1:Khat)
  
  if(membership %in% c('softmax')){
    if(Khat == 2) {
      family <- "binomial"
    } else {
      family <- "multinomial"
    }
    
    # tuning penalty for softmax model
    cv_model <- cv.glmnet(x = X_train[,], y = labels, 
                          intercept = FALSE,
                          family = family,
                          alpha = 0.5,  # enet
                          nfolds = 5,
                          type.measure = "class",
                          lambda =  10^seq(3, -4, length.out = 100))
    lambda_tuned <- cv_model$lambda.min
    cat(paste0('lambda_tuned = ', lambda_tuned))
    
    model <- glmnet(x = X_train[,], y = labels, 
                    intercept = FALSE,
                    family = family,
                    alpha = 0.5,   # enet
                    lambda = lambda_tuned)
    
    coef_min <- coef(model)
    labels_pred <- predict(model, newx = X_test, type = 'class', s = model$lambda)
    if(is.matrix(labels_pred)) {
      labels_pred <- as.vector(labels_pred)
    }
    labels_pred <- factor(labels_pred, levels = 1:Khat)
    U_test <- Ucenters[as.numeric(labels_pred), , drop = FALSE]
    y_pred <- rowSums(X_test * (U_test %*% t(Vhat)))
    
  } else if(membership %in% c('KNN')){
    Z_train <- as.matrix(X_train) %*% Vhat # use the reduced features rather than X
    Z_test <- as.matrix(X_test) %*% Vhat
    labels_pred <- knn(train = Z_train, test = Z_test, cl = labels, k = 10)
    U_test <- Ucenters[as.numeric(labels_pred), , drop = FALSE]
    y_pred <- rowSums(X_test * (U_test %*% t(Vhat)))
    labels_pred <- factor(labels_pred, levels = 1:Khat)
    
  } else{
    stop("membership must be either 'softmax' or 'KNN'")
  }
  return(y_pred)
}
  


predict_FMM <- function(X_test, X_train, y_train, method, r, fmmM_lst, alpha, aggregate = TRUE){
  ## estimate and predict using FMM
  ## method: "FMMlasso" or "FMMenet"

  X_train_df <- as.data.frame(X_train)
  X_test_df <- as.data.frame(X_test)
  colnames(X_train_df) <- paste0("X", 1:ncol(X_train))
  colnames(X_test_df) <- paste0("X", 1:ncol(X_test))
  train_data <- cbind(y = y_train, X_train_df)
  
  formula_str <- paste("y ~", paste(colnames(X_train_df), collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  if(method %in% c("FMMlasso")) {
    temp <- stepFlexmix(formula_obj, data = train_data, k = fmmM_lst, nrep = 2, drop = TRUE,
                        model = FLXMRglmnet(family = "gaussian", alpha = 1, 
                                            intercept = FALSE, adaptive = FALSE),
                        control = list(iter.max = 300))
  }
  if(method %in% c("FMMenet")) {
    temp <- stepFlexmix(formula_obj, data = train_data, k = fmmM_lst, nrep = 1, drop = TRUE,
                        model = FLXMRglmnet(family = "gaussian", alpha = 0.5, cluster=labels_train,
                                            intercept = FALSE, adaptive = FALSE),
                        control = list(iter.max = 300))
  }
 
  model <- getModel(temp, which = "ICL")  # select the best model
  y_pred <- predict(model, newdata = X_test_df, aggregate = aggregate)
  y_pred <- as.numeric(y_pred[[1]])
  return(y_pred)
}
