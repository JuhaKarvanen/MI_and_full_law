# File: simulation_examples_fig1.R
# Author: Juha Karvanen and Santtu Tikka
# Date: 2025-10-28
# Summary: R code to reproduce the simulation results for Examples 1a - 1f
# of the paper J. Karvanen, S. Tikka (2025), Multiple imputation and full law identifiability.

library("mice", quietly = TRUE, warn.conflicts = FALSE)
library("mgcv")
expit <- function(x) exp(x) / (1 + exp(x))
n <- 1000000
reporder <- c("statistic", "truth", "mi", "miri", "fmi", "cca", "aca", "eq")
set.seed(211020251)

# Apply all imputation methods, CCA and ACA and compile results
compute_results <- function(da, truth = c(0, 0, 1, sqrt(2), sqrt(2) / 2), formulas = FALSE) {
  res <- data.frame(
    statistic = c("$\\E(X)$", "$\\E(Y)$", "$\\sd(X)$", "$\\sd(Y)$", "$\\Cor(X,Y)$"),
    truth = truth,
    mi = rep(NA, 5),
    miri = rep(NA,5),
    fmi = rep(NA, 5),
    cca = rep(NA, 5),
    aca = rep(NA, 5),
    eq = rep(NA, 5)
  )
  rownames(res) <- c("x", "y", "sx", "sy", "r")

  # Standard multiple imputation (without response indicators, MI)
  im1 <- mice(da[, c("x", "y")], method = c("norm", "norm"), printFlag = FALSE)
  res[, "mi"] <- mi_results(im1)

  # Multiple imputation with response indicators (MIRI)
  if (formulas) {
    im2 <- mice(
      da,
      method = c("norm", "norm", "mean", "mean"),
      formulas = list(x = x ~ y * ry, y = y ~ x * rx, rx = rx ~ 1, ry = ry ~ 1),
      visitSequence = c("x", "y", "rx", "ry"),
      printFlag = FALSE
    )
  } else {
    predmatrix <- make.predictorMatrix(da)
    predmatrix["x", "rx"] <- 0
    predmatrix["y", "ry"] <- 0
    im2 <- mice(
      da,
      method = c("norm", "norm", "mean", "mean"),
      predictorMatrix = predmatrix,
      visitSequence = c("x", "y", "rx", "ry"),
      printFlag = FALSE
    )
  }
  res[, "miri"] <- mi_results(im2)

  # Factorizable multiple imputation (FMI)
  daf <- da
  daf$y[rx == 0 & ry == 1] <- NA
  im3 <- mice(
    daf[, c("x", "y")],
    method = c("sample", "norm"),
    visitSequence = c("x", "y"),
    maxit = 1,
    printFlag = FALSE
  )
  res[, "fmi"] <- mi_results(im3)

  # Complete case analysis
  dacc <- da[rx & ry, ]
  res$cca <- c(
    mean(dacc$x),
    mean(dacc$y),
    sd(dacc$x),
    sd(dacc$y),
    cor(dacc$x, dacc$y)
  )

  # Available case analysis
  res$aca <- c(
    mean(da$x, na.rm = TRUE),
    mean(da$y, na.rm = TRUE),
    sd(da$x, na.rm = TRUE),
    sd(da$y, na.rm = TRUE),
    cor(da$x, da$y, use = "pairwise.complete.obs")
  )

  res
}

# Results from Multiple Imputation
mi_results <- function(im) {
  meanx <- with(im, mean(x))
  meany <- with(im, mean(y))
  sdx <- with(im, sd(x))
  sdy <- with(im, sd(y))
  cor <- with(im, cor(x, y))
  c(
    mean(unlist(meanx$analyses)),
    mean(unlist(meany$analyses)),
    mean(unlist(sdx$analyses)),
    mean(unlist(sdy$analyses)),
    mean(unlist(cor$analyses))
  )
}

# Format table
format_table <- function(x, bias = FALSE) {
  if (bias) {
    for (col in c("mi", "miri", "fmi", "cca", "aca")) {
      x[[col]] <- paste0("$", formatC(x[[col]] - x$truth, digits = 2, format = "f"), "$")
    }
    if (any(is.na(x$eq))) {
      x$eq <- "NA"
    } else {
      x$eq <- paste0("$", formatC(x$eq - x$truth, digits = 2, format = "f"), "$")
    }
    x$truth <- paste0("$", formatC(x$truth, digits = 2, format = "f"), "$")
  } else {
    for (col in c("mi", "miri", "fmi", "cca", "aca")) {
      x[[col]] <- paste0("$", formatC(x[[col]], digits = 2, format = "f"), "$")
    }
    if (any(is.na(x$eq))) {
      x$eq <- "NA"
    } else {
      x$eq <- paste0("$", formatC(x$eq, digits = 2, format = "f"), "$")
    }
  }
  reporder <- c("statistic", "truth", "mi", "miri", "fmi", "cca", "aca", "eq")
  x[, reporder]
}

# Reporting
latextable <- function(td) {
  tddim <- dim(td)
  lt <- "\n"
  for (i in 1:tddim[1]) {
    for (j in 1:tddim[2]) {
      if (j == 1) lt <- paste(lt, td[i, j])
      else lt <- paste(lt, " & ", td[i, j])
    }
    if (i == tddim[1]) lt <- paste(lt, "\n")
    else lt <- paste(lt, "\\\\ \n")
  }
  cat(lt)
}

############################################################################
# Example 1a
############################################################################

x1 <- rnorm(n,0,1)
y1 <- x1 + rnorm(n, 0, 1)
rx <- as.numeric(0.7 > runif(n, 0, 1))
ry <- as.numeric(expit(x1) > runif(n, 0, 1))
x <- x1
x[!rx] <- NA
y <- y1
y[!ry] <- NA
da <- data.frame(x = x, y = y, rx = rx, ry = ry)

results_1a <- compute_results(da)

# Estimation with identifying formula
dax <- da[da$rx == 1, ]
dacc <- da[da$rx & da$ry, ]
yxmodel <- lm(y ~ x, data = dacc)
yxsigma <- summary(yxmodel)$sigma
ypred <- predict(yxmodel, newdata = dax)
ysim <- ypred + rnorm(length(ypred), 0, yxsigma)
results_1a["x", "eq"] <- mean(dax$x)
results_1a["sx", "eq"] <- sd(dax$x)
results_1a["y", "eq"] <- mean(ysim)
results_1a["sy", "eq"] <- sd(ysim)
results_1a["r", "eq"] <- cor(dax$x, ysim)

tab_1a <- format_table(results_1a)
tab_1ab <- format_table(results_1a, bias = TRUE)

cat("\nExample 1a:\n")
cat(reporder, "\n")
cat(latextable(tab_1a))
cat(reporder, "\n")
cat(latextable(tab_1ab))


############################################################################
# Example 1b
############################################################################

x1 <- rnorm(n, 0, 1)
y1 <- x1 + rnorm(n, 0, 1)
ry <- as.numeric(expit(x1) > runif(n, 0, 1))
rx <- as.numeric(0.4 + 0.3 * ry  > runif(n, 0, 1))
x <- x1
x[!rx] <- NA
y <- y1
y[!ry] <- NA
da <- data.frame(x = x, y = y, rx = rx, ry = ry)

results_1b <- compute_results(da)

# Estimation with identifying formula
dax <- da[da$rx == 1, ]
xmodel0 <- lm(x ~ 1, data = da[da$ry == 0, ])
xmodel1 <- lm(x ~ y, data = da[da$ry == 1, ])
xsigma0 <- summary(xmodel0)$sigma
xsigma1 <- summary(xmodel1)$sigma
xpred0 <- predict(xmodel0, newdata = da[da$ry == 0, ])
xpred1 <- predict(xmodel1, newdata = da[da$ry == 1, ])
xsim0 <-  xpred0 + rnorm(length(xpred0), 0, xsigma0)
xsim1 <-  xpred1 + rnorm(length(xpred1), 0, xsigma1)
xsim <- c(xsim0, xsim1)
yxmodel <- lm(y ~ x, data = dacc)
yxsigma <- summary(yxmodel)$sigma
ypred <- predict(yxmodel, newdata = data.frame(x = xsim))
ysim <- ypred + rnorm(length(ypred), 0, yxsigma)
results_1b["x", "eq"] <- mean(xsim)
results_1b["sx", "eq"] <- sd(xsim)
results_1b["y", "eq"] <- mean(ysim)
results_1b["sy", "eq"] <- sd(ysim)
results_1b["r", "eq"] <- cor(xsim, ysim)

# Reporting
tab_1b <- format_table(results_1b)
tab_1bb <- format_table(results_1b, bias = TRUE)

cat("\nExample 1b:\n")
cat(reporder, "\n")
cat(latextable(tab_1b))
cat(reporder, "\n")
cat(latextable(tab_1bb[, reporder]))


############################################################################
# Example 1c
############################################################################

x1 <- rnorm(n, 0, 1)
y1 <- x1 + rnorm(n, 0, 1)
rx <- as.numeric(expit(y1) > runif(n, 0, 1))
ry <- as.numeric(expit(x1) > runif(n, 0, 1))
x <- x1
x[!rx] <- NA
y <- y1
y[!ry] <- NA
da <- data.frame(x = x, y = y, rx = rx, ry = ry)

results_1c <- compute_results(da)

# Estimation with identifying formula
dax <- da[da$rx == 1, ]
day <- da[da$ry == 1, ]
dacc <- da[da$rx & da$ry, ]
yxmodel <- gam(y ~ s(x), data = dacc)
yxsigma <- sqrt(yxmodel$sig2)
ypred <- predict(yxmodel, newdata = dax)
ysim <- ypred + rnorm(length(ypred), 0, yxsigma)
ysim <- ifelse(dax$ry == 1, dax$y, ysim) # Comment for simulation only; uncomment for imputation
rxmodel <- glm(rx ~ y, family = binomial, data = day)
rxpred <- predict(rxmodel, type = "response", newdata = data.frame(y = ysim))
weights <- 1 / rxpred
sumweights <- sum(weights)
weightedcov <- cov.wt(data.frame(x = dax$x, y = ysim), wt = weights, cor = TRUE)
weightedsd <- sqrt(diag(weightedcov$cov))
results_1c["x", "eq"] <- sum(dax$x * weights) / sumweights
results_1c["sx", "eq"] <- weightedsd["x"]
results_1c["y", "eq"] <- sum(ysim * weights) / sumweights
results_1c["sy", "eq"] <- weightedsd["y"]
results_1c["r", "eq"] <- weightedcov$cor["x", "y"]

# Reporting
tab_1c <- format_table(results_1c)
tab_1cb <- format_table(results_1c, bias = TRUE)

cat("\nExample 1c:\n")
cat(reporder, "\n")
cat(latextable(tab_1c))
cat(reporder, "\n")
cat(latextable(tab_1cb))


############################################################################
# Example 1d
############################################################################

x1 <- rnorm(n, 0, 1)
y1 <- x1 + rnorm(n, 0, 1)
rx <- as.numeric(0.7 > runif(n, 0, 1))
ry <- as.numeric(expit((x1 + 0.4) * (2 * rx - 1)) > runif(n, 0, 1))
x <- x1
x[!rx] <- NA
y <- y1
y[!ry] <- NA
da <- data.frame(x = x, y = y, rx = rx, ry = ry)

results_1d <- compute_results(da)

# Estimation with identifying formula
dax <- da[da$rx == 1, ]
yxmodel <- lm(y ~ x, data = dacc)
yxsigma <- summary(yxmodel)$sigma
ypred <- predict(yxmodel, newdata = dax)
ysim <- ypred + rnorm(length(ypred), 0, yxsigma)
results_1d["x", "eq"] <- mean(dax$x)
results_1d["sx", "eq"] <- sd(dax$x)
results_1d["y", "eq"] <- mean(ysim)
results_1d["sy", "eq"] <- sd(ysim)
results_1d["r", "eq"] <- cor(dax$x, ysim)

# Reporting
tab_1d <- format_table(results_1d)
tab_1db <- format_table(results_1d, bias = TRUE)

cat("\nExample 1d:\n")
cat(reporder, "\n")
cat(latextable(tab_1d))
cat(reporder,"\n")
cat(latextable(tab_1db))


############################################################################
# Example 1e
############################################################################

x1 <- rnorm(n, 0, 1)
y1 <- x1 + rnorm(n, 0, 1)
rx <- as.numeric(expit(y1) > runif(n, 0, 1))
ry <- as.numeric(expit((x1 + 0.4) * (2 * rx - 1)) > runif(n, 0, 1))
x <- x1
x[!rx] <- NA
y <- y1
y[!ry] <- NA
da <- data.frame(x = x, y = y, rx = rx, ry = ry)

results_1e <- compute_results(da)

# Estimation with identifying formula
results_1e["x", "eq"] <- NA
results_1e["sx", "eq"] <- NA
results_1e["y", "eq"] <- NA
results_1e["sy", "eq"] <- NA
results_1e["r", "eq"] <- NA

# Reporting
tab_1e <- format_table(results_1e)
tab_1eb <- format_table(results_1e, bias = TRUE)

cat("\nExample 1e:\n")
cat(reporder,"\n")
cat(latextable(tab_1e))
cat(reporder,"\n")
cat(latextable(tab_1eb))


############################################################################
# Figure 1f
############################################################################

x1 <- rnorm(n, 0, 1)
y1 <- rnorm(n, 0, sqrt(2))
rx <- as.numeric(expit(y1) > runif(n, 0, 1))
ry <- as.numeric(expit((x1 + 0.4) * (2 * rx - 1)) > runif(n, 0, 1))
x <- x1
x[!rx] <- NA
y <- y1
y[!ry] <- NA
da <- data.frame(x = x, y = y, rx = rx, ry = ry)

results_1f <- compute_results(da, truth = c(0, 0, 1, sqrt(2), 0), formulas = TRUE)

# Estimation with identifying formula
dax <- da[da$rx == 1, ]
results_1f["x", "eq"] <- mean(dax$x)
results_1f["sx", "eq"] <- sd(dax$x)
yxmodel <- lm(y ~ rx, data = da)
yxsigma <- summary(yxmodel)$sigma
ypred <- predict(yxmodel, newdata = da)
ysim <- ypred + rnorm(length(ypred), 0, yxsigma)
results_1f["y", "eq"] <- mean(ysim)
results_1f["sy", "eq"] <- sd(ysim)
results_1f["r", "eq"] <- cor(dax$x, ysim[da$rx == 1])

# Reporting
tab_1f <- format_table(results_1f)
tab_1fb <- format_table(results_1f, bias = TRUE)

cat("\nExample 1f:\n")
cat(reporder, "\n")
cat(latextable(tab_1f))
cat(reporder, "\n")
cat(latextable(tab_1fb))


# All bias tables
cat(reporder, "\n")
cat(latextable(tab_1ab))
cat(latextable(tab_1bb))
cat(latextable(tab_1cb))
cat(latextable(tab_1db))
cat(latextable(tab_1eb))
cat(latextable(tab_1fb))

# save.image(file = "simulation_results_fig1_20251021.Rdata")
