# File: simulation_examples.R
# Author: Juha Karvanen and Santtu Tikka
# Date: 2025-10-28
# Summary: R code to reproduce the simulation results for Example 2a and 2b
# of the paper J. Karvanen, S. Tikka (2025), Multiple imputation and full law identifiability.

library("mice", quietly = TRUE, warn.conflicts = FALSE)
library("mgcv")

expit <- function(x) exp(x) / (1 + exp(x))
n <- 10000000
set.seed(211020252)

latextable <- function(td) {
  tddim <- dim(td)
  lt <- "\n"
  for (i in 1:tddim[1]) {
    for (j in 1:tddim[2]) {
      if (j == 1) lt <- paste(lt, td[i, j])
      else lt <- paste(lt," & ", td[i, j])
    }
    if (i == tddim[1]) lt <- paste(lt, "\n")
    else lt <- paste(lt, "\\\\ \n")
  }
  cat(lt)
}

repordera <- c("statistic", "truth", "mi", "miri", "fmi", "cca", "aca", "eq")
reporderb <- c("statistic", "truth", "mi", "miri", "fmi", "cca", "aca", "eqA", "eqB")


############################################################################
# Example 2a
############################################################################

results_2a <- data.frame(
  statistic = c(
    "$\\E(X)$", "$\\E(W)$", "$\\E(Z)$", "$\\E(Y)$",
    "$\\sd(X)$", "$\\sd(W)$", "$\\sd(Z)$", "$\\sd(Y)$",
    "$\\Cor(X,W)$", "$\\Cor(X,Z)$", "$\\Cor(X,Y)$",
    "$\\Cor(W,Z)$", "$\\Cor(W,Y)$", "$\\Cor(Z,Y)$"
  ),
  truth = c(
    0, 0, 0, 0, 1, 1, 1, 1, sqrt(2)/2, 0.5, 3/4, sqrt(2)/2, sqrt(2)/2, 3/4
  ),
  fmi = rep(NA, 14),
  mi = rep(NA, 14),
  cca = rep(NA, 14),
  aca = rep(NA, 14),
  eq = rep(NA, 14)
)
rownames(results_2a) <- c(
  "x", "w", "z", "y", "sx", "sw", "sz", "sy",
  "rxw", "rxz", "rxy", "rwz", "rwy", "rzy"
)

x1 <- rnorm(n, 0, 1)
w1 <- (x1 + rnorm(n, 0, 1)) / sqrt(2)
z1 <- (w1 + rnorm(n, 0, 1)) / sqrt(2)
y1 <- (x1 + z1 + rnorm(n, 0, 1)) / sqrt(4)
da1 <- data.frame(z = z1, w = w1, x = x1,  y = y1)
const1 <- 0
const2 <- 0
coef <- 1
rx <- as.numeric(coef * expit(w1 + const2) + const1 > runif(n, 0, 1))
rz <- as.numeric(0.7  > runif(n, 0, 1))
rw <- as.numeric(coef * expit((z1 + 0.4) * (2 * rz - 1) + const2) + const1 > runif(n, 0, 1))
ry <- as.numeric(coef * expit(x1 + const2) + const1 > runif(n, 0, 1))
x <- x1
x[!rx] <- NA
z <- z1
z[!rz] <- NA
w <- w1
w[!rw] <- NA
y <- y1
y[!ry] <- NA
da <- data.frame(x = x, w = w, z = z, y = y, rx = rx, rw = rw, rz = rz, ry = ry)

# Standard multiple imputation (without response indicators)
im_2a <- mice(
  da[, c("z", "w", "x", "y")],
  method = c("norm", "norm", "norm", "norm"),
  printFlag = FALSE
)
meanx <- with(im_2a, mean(x))
meanw <- with(im_2a, mean(w))
meanz <- with(im_2a, mean(z))
meany <- with(im_2a, mean(y))
sdx <- with(im_2a, sd(x))
sdw <- with(im_2a, sd(w))
sdz <- with(im_2a, sd(z))
sdy <- with(im_2a, sd(y))
corxw <- with(im_2a, cor(x, w))
corxz <- with(im_2a, cor(x, z))
corxy <- with(im_2a, cor(x, y))
corwz <- with(im_2a, cor(w, z))
corwy <- with(im_2a, cor(w, y))
corzy <- with(im_2a, cor(z, y))

results_2a["x", "mi"] <- mean(unlist(meanx$analyses))
results_2a["w", "mi"] <- mean(unlist(meanw$analyses))
results_2a["z", "mi"] <- mean(unlist(meanz$analyses))
results_2a["y", "mi"] <- mean(unlist(meany$analyses))
results_2a["sx", "mi"] <- mean(unlist(sdx$analyses))
results_2a["sw", "mi"] <- mean(unlist(sdw$analyses))
results_2a["sz", "mi"] <- mean(unlist(sdz$analyses))
results_2a["sy", "mi"] <- mean(unlist(sdy$analyses))
results_2a["rxw", "mi"] <- mean(unlist(corxw$analyses))
results_2a["rxz", "mi"] <- mean(unlist(corxz$analyses))
results_2a["rxy", "mi"] <- mean(unlist(corxy$analyses))
results_2a["rwz", "mi"] <- mean(unlist(corwz$analyses))
results_2a["rwy", "mi"] <- mean(unlist(corwy$analyses))
results_2a["rzy", "mi"] <- mean(unlist(corzy$analyses))

# Multiple imputation with response indicators (MIRI)
im_2a <- mice(
  da,
  method = c("norm", "norm", "norm", "norm", "mean", "mean", "mean", "mean"),
  formulas = list(
    x = x ~ (z + w + y) * rz * rw * ry, w = w ~ (z + x + y) * rz * rx * ry,
    z = z ~ (x + w + y) * rx * rw * ry, y = y ~ (z + w + x) * rz * rw *rx,
    rz = rz ~ 1, rw = rw ~ 1, rx = rx ~ 1, ry = ry ~ 1
  ),
  visitSequence = c("x", "w", "z", "y", "rz", "rw", "rx", "ry"),
  printFlag = FALSE
)
meanx <- with(im_2a, mean(x))
meanw <- with(im_2a, mean(w))
meanz <- with(im_2a, mean(z))
meany <- with(im_2a, mean(y))
sdx <- with(im_2a, sd(x))
sdw <- with(im_2a, sd(w))
sdz <- with(im_2a, sd(z))
sdy <- with(im_2a, sd(y))
corxw <- with(im_2a, cor(x, w))
corxz <- with(im_2a, cor(x, z))
corxy <- with(im_2a, cor(x, y))
corwz <- with(im_2a, cor(w, z))
corwy <- with(im_2a, cor(w, y))
corzy <- with(im_2a, cor(z, y))

results_2a["x", "miri"] <- mean(unlist(meanx$analyses))
results_2a["w", "miri"] <- mean(unlist(meanw$analyses))
results_2a["z", "miri"] <- mean(unlist(meanz$analyses))
results_2a["y", "miri"] <- mean(unlist(meany$analyses))
results_2a["sx", "miri"] <- mean(unlist(sdx$analyses))
results_2a["sw", "miri"] <- mean(unlist(sdw$analyses))
results_2a["sz", "miri"] <- mean(unlist(sdz$analyses))
results_2a["sy", "miri"] <- mean(unlist(sdy$analyses))
results_2a["rxw", "miri"] <- mean(unlist(corxw$analyses))
results_2a["rxz", "miri"] <- mean(unlist(corxz$analyses))
results_2a["rxy", "miri"] <- mean(unlist(corxy$analyses))
results_2a["rwz", "miri"] <- mean(unlist(corwz$analyses))
results_2a["rwy", "miri"] <- mean(unlist(corwy$analyses))
results_2a["rzy", "miri"] <- mean(unlist(corzy$analyses))

# Factorizable multiple imputation (FMI)
daf <- da
daf$w[daf$rz == 0] <- NA
daf$x[daf$rz * daf$rw == 0] <- NA
daf$y[daf$rz * daf$rw * daf$rx == 0] <- NA
daf <- daf[, c("z", "w", "x", "y")]

im_2a <- mice(
  daf,
  method = c("sample", "norm", "norm", "norm"),
  formulas = list(z = z ~ 1, w = w ~ z, x = x ~ w, y = y ~ x + z),
  visitSequence = c("z", "w", "x", "y"),
  maxit = 1,
  printFlag = FALSE
)
meanx <- with(im_2a, mean(x))
meanw <- with(im_2a, mean(w))
meanz <- with(im_2a, mean(z))
meany <- with(im_2a, mean(y))
sdx <- with(im_2a, sd(x))
sdw <- with(im_2a, sd(w))
sdz <- with(im_2a, sd(z))
sdy <- with(im_2a, sd(y))
corxw <- with(im_2a, cor(x, w))
corxz <- with(im_2a, cor(x, z))
corxy <- with(im_2a, cor(x, y))
corwz <- with(im_2a, cor(w, z))
corwy <- with(im_2a, cor(w, y))
corzy <- with(im_2a, cor(z, y))

results_2a["x", "fmi"] <- mean(unlist(meanx$analyses))
results_2a["w", "fmi"] <- mean(unlist(meanw$analyses))
results_2a["z", "fmi"] <- mean(unlist(meanz$analyses))
results_2a["y", "fmi"] <- mean(unlist(meany$analyses))
results_2a["sx", "fmi"] <- mean(unlist(sdx$analyses))
results_2a["sw", "fmi"] <- mean(unlist(sdw$analyses))
results_2a["sz", "fmi"] <- mean(unlist(sdz$analyses))
results_2a["sy", "fmi"] <- mean(unlist(sdy$analyses))
results_2a["rxw", "fmi"] <- mean(unlist(corxw$analyses))
results_2a["rxz", "fmi"] <- mean(unlist(corxz$analyses))
results_2a["rxy", "fmi"] <- mean(unlist(corxy$analyses))
results_2a["rwz", "fmi"] <- mean(unlist(corwz$analyses))
results_2a["rwy", "fmi"] <- mean(unlist(corwy$analyses))
results_2a["rzy", "fmi"] <- mean(unlist(corzy$analyses))

# Complete case analysis
dacc <- da[rx & rw & rz & ry, ]
results_2a["x", "cca"] <- mean(dacc$x)
results_2a["w", "cca"] <- mean(dacc$w)
results_2a["z", "cca"] <- mean(dacc$z)
results_2a["y", "cca"] <- mean(dacc$y)
results_2a["sx", "cca"] <- sd(dacc$x)
results_2a["sw", "cca"] <- sd(dacc$w)
results_2a["sz", "cca"] <- sd(dacc$z)
results_2a["sy", "cca"] <- sd(dacc$y)
results_2a["rxw", "cca"] <- cor(dacc$x, dacc$w)
results_2a["rxz", "cca"] <- cor(dacc$x, dacc$z)
results_2a["rxy", "cca"] <- cor(dacc$x, dacc$y)
results_2a["rwz", "cca"] <- cor(dacc$w, dacc$z)
results_2a["rwy", "cca"] <- cor(dacc$w, dacc$y)
results_2a["rzy", "cca"] <- cor(dacc$z, dacc$y)

# Available case analysis
results_2a["x", "aca"] <- mean(da$x, na.rm = TRUE)
results_2a["w", "aca"] <- mean(da$w, na.rm = TRUE)
results_2a["z", "aca"] <- mean(da$z, na.rm = TRUE)
results_2a["y", "aca"] <- mean(da$y, na.rm = TRUE)
results_2a["sx", "aca"] <- sd(da$x, na.rm = TRUE)
results_2a["sw", "aca"] <- sd(da$w, na.rm = TRUE)
results_2a["sz", "aca"] <- sd(da$z, na.rm = TRUE)
results_2a["sy", "aca"] <- sd(da$y, na.rm = TRUE)
results_2a["rxw", "aca"] <- cor(da$x, da$w, use = "pairwise.complete.obs")
results_2a["rxz", "aca"] <- cor(da$x, da$z, use = "pairwise.complete.obs")
results_2a["rxy", "aca"] <- cor(da$x, da$y, use = "pairwise.complete.obs")
results_2a["rwz", "aca"] <- cor(da$w, da$z, use = "pairwise.complete.obs")
results_2a["rwy", "aca"] <- cor(da$w, da$y, use = "pairwise.complete.obs")
results_2a["rzy", "aca"] <- cor(da$z, da$y, use = "pairwise.complete.obs")

# Reporting

tab_2a <- results_2a

tab_2a$truth <- paste0("$", formatC(tab_2a$truth, digits = 2, format = "f"), "$")
tab_2a$fmi <- paste0("$", formatC(tab_2a$fmi, digits = 2, format = "f"), "$")
tab_2a$mi <- paste0("$", formatC(tab_2a$mi, digits = 2, format = "f"), "$")
tab_2a$miri <- paste0("$", formatC(tab_2a$miri, digits = 2, format = "f"), "$")
tab_2a$cca <- paste0("$", formatC(tab_2a$cca, digits = 2, format = "f"), "$")
tab_2a$aca <- paste0("$", formatC(tab_2a$aca, digits = 2, format = "f"), "$")

tab_2ab <- results_2a
tab_2ab$fmi <- formatC(tab_2ab$fmi - tab_2ab$truth, digits = 2, format = "f")
tab_2ab$mi <- formatC(tab_2ab$mi- tab_2ab$truth, digits = 2, format = "f")
tab_2ab$miri <- formatC(tab_2ab$miri - tab_2ab$truth, digits = 2, format = "f")
tab_2ab$cca <- formatC(tab_2ab$cca - tab_2ab$truth, digits = 2, format = "f")
tab_2ab$aca <- formatC(tab_2ab$aca - tab_2ab$truth, digits = 2, format = "f")
tab_2ab$truth <- formatC(tab_2ab$truth, digits = 2, format = "f")

tab_2ab$truth <- paste0("$", formatC(tab_2ab$truth, digits = 2, format = "f"), "$")
tab_2ab$midecomp <- paste0("$", formatC(tab_2ab$midecomp, digits = 2, format = "f"), "$")
tab_2ab$mi <- paste0("$", formatC(tab_2ab$mi, digits = 2, format = "f"), "$")
tab_2ab$cca <- paste0("$", formatC(tab_2ab$cca, digits = 2, format = "f"), "$")
tab_2ab$aca <- paste0("$", formatC(tab_2ab$aca, digits = 2, format = "f"), "$")


cat("\nExample 2a:\n")
cat(repordera,"\n")
cat(latextable(tab_2a[,repordera]))
cat(repordera,"\n")
cat(latextable(tab_2ab[,repordera]))


############################################################################
# Example 2b
############################################################################

results_2b <- data.frame(
  statistic = c(
    "$\\E(X)$","$\\E(W)$","$\\E(Z)$","$\\E(Y)$",
    "$\\sd(X)$","$\\sd(W)$","$\\sd(Z)$","$\\sd(Y)$",
    "$\\Cor(X,W)$","$\\Cor(X,Z)$","$\\Cor(X,Y)$",
    "$\\Cor(W,Z)$","$\\Cor(W,Y)$","$\\Cor(Z,Y)$"
  ),
  truth = c(
    0, 0, 0, 0, 1, 1, 1, 1, sqrt(2)/2, 0.5, 3/4, sqrt(2)/2, sqrt(2)/2, 3/4
  ),
  mi = rep(NA, 14),
  miri = rep(NA, 14),
  fmi = rep(NA, 14),
  cca = rep(NA, 14),
  aca = rep(NA, 14),
  eqA = rep(NA, 14),
  eqB = rep(NA, 14)
)
rownames(results_2b) <- c(
  "x", "w", "z", "y", "sx", "sw", "sz", "sy",
  "rxw", "rxz", "rxy", "rwz", "rwy", "rzy"
)

x1 <- rnorm(n, 0, 1)
w1 <- (x1 + rnorm(n, 0, 1)) / sqrt(2)
z1 <- (w1 + rnorm(n, 0, 1)) / sqrt(2)
y1 <- (x1 + z1 + rnorm(n, 0, 1)) / sqrt(4)
da1 <- data.frame(z = z1, w = w1, x = x1, y = y1)
const1 <- 0
const2 <- 0
coef <- 1
rx <- as.numeric(coef * expit(w1 + const2) + const1  > runif(n, 0, 1))
rz <- as.numeric(coef * expit(x1 + const2) + const1  > runif(n, 0, 1))
rw <- as.numeric(coef * expit((z1 + 0.4) * (2 * rz - 1) + const2) + const1 > runif(n, 0, 1))
ry <- as.numeric(coef * expit(x1 + const2) + const1 > runif(n,0,1))
x <- x1
x[!rx] <- NA
z <- z1
z[!rz] <- NA
w <- w1
w[!rw] <- NA
y <- y1
y[!ry] <- NA
da <- data.frame(x = x, w = w, z = z, y = y,  rx = rx, rw = rw, rz = rz, ry = ry)

# # Estimation with identifying formula
# # p(X,W,Z,Y) = \frac{P(X, W, Z, Y, R_X = 1, R_W = 1, R_Z = 1, R_Y = 1)}{p(R_X = 1|W,R_W = 1)p(R_W = 1|Z,R_Z = 1)p(R_Z = 1,R_Y = 1|X,R_X = 1)}.
dacc <- da[rx & rw & rz & ry, ]

ncc <- nrow(dacc)
modelrx <- gam(rx ~ s(w), family = binomial(), data = da[da$rw == 1,])
modelrw <- gam(rw ~ s(z), family = binomial(), data = da[da$rz == 1,])
modelrzry <- gam(I(rz*ry) ~ s(x), family = binomial(), data = da[da$rx == 1,])

rxpred <- predict(modelrx, newdata = dacc, type = "response")
rwpred <- predict(modelrw, newdata = dacc, type = "response")
rzrypred <- predict(modelrzry, newdata = dacc, type = "response")

rweights <- 1 / (rxpred * rwpred * rzrypred)
rawsumrweights <- sum(rweights)
rweights <- ncc * rweights / rawsumrweights
sumrweights <- sum(rweights)

results_2b["x", "eqA"] <- sum(dacc$x * rweights) / sumrweights
results_2b["w", "eqA"] <- sum(dacc$w * rweights) / sumrweights
results_2b["z", "eqA"] <- sum(dacc$z * rweights) / sumrweights
results_2b["y", "eqA"] <- sum(dacc$y * rweights) / sumrweights
weightedcov <- cov.wt(dacc[, c("x","w","z","y")], wt = as.vector(rweights), cor = TRUE)
weightedsd <- sqrt(diag(weightedcov$cov))

results_2b["sx", "eqA"] <- weightedsd["x"]
results_2b["sw", "eqA"] <- weightedsd["w"]
results_2b["sz", "eqA"] <- weightedsd["z"]
results_2b["sy", "eqA"] <- weightedsd["y"]
results_2b["rxw", "eqA"] <- weightedcov$cor["x", "w"]
results_2b["rxz", "eqA"] <- weightedcov$cor["x", "z"]
results_2b["rxy", "eqA"] <- weightedcov$cor["x", "y"]
results_2b["rwz", "eqA"] <- weightedcov$cor["w", "z"]
results_2b["rwy", "eqA"] <- weightedcov$cor["w", "y"]
results_2b["rzy", "eqA"] <- weightedcov$cor["z", "y"]


# Estimation with identifying formula
# p(X,W,Z,Y) = \frac{p(R_X=1)p(X|R_X = 1)p(Z|X,R_X R_Z = 1)p(W,Y|X,Z,R_X R_W R_Z R_Y = 1)}{p(R_X = 1|W,R_W = 1)}
dax <- da[da$rx == 1,]
daw <- da[da$rw == 1,]
daxz <- da[da$rx == 1 & da$rz == 1,]
dacc <- da[da$rx == 1 & da$rw == 1 & da$rz == 1 & da$ry == 1,]

wmodel <- lm(w ~ x + z, data = dacc)
zmodel <- lm(z ~ x, data = daxz)
ymodel <- lm(y ~ x + z, data = dacc)
rxmodel <- glm(rx ~ w, family = binomial(), data = daw)

wsigma <- summary(wmodel)$sigma
zsigma <- summary(zmodel)$sigma
ysigma <- summary(ymodel)$sigma

ccind <- dax$rw == 1 & dax$rz == 1 & dax$ry == 1

zpred <- predict(zmodel, newdata = data.frame(x = dax$x))
zsim <- zpred + rnorm(length(zpred),0,zsigma)
zsim <- ifelse(dax$rz == 1, dax$z, zsim) # Comment for simulation only; uncomment for imputation

wpred <- predict(wmodel, newdata = data.frame(x = dax$x, z = zsim))
wsim <- wpred + rnorm(length(wpred),0,wsigma)
wsim <- ifelse(dax$rw == 1 & dax$rz == 1 & dax$ry == 1, dax$w, wsim) # Comment for simulation only; uncomment for imputation

ypred <- predict(ymodel, newdata = data.frame(x = dax$x, z = zsim))
ysim <- ypred + rnorm(length(ypred),0,ysigma)
ysim <- ifelse(dax$rw == 1 & dax$rz == 1 & dax$ry == 1, dax$y, ysim) # Comment for simulation only; uncomment for imputation

dasim <- data.frame(x = dax$x, w = wsim, z = zsim, y = ysim)

rxpred <- predict(rxmodel, newdata = data.frame(w = dasim$w), type = "response")
weights <- 1/rxpred
sumweights <- sum(weights)

results_2b["x", "eqB"] <- sum(dasim$x * weights) / sumweights
results_2b["w", "eqB"] <- sum(dasim$w * weights) / sumweights
results_2b["z", "eqB"] <- sum(dasim$z * weights) / sumweights
results_2b["y", "eqB"] <- sum(dasim$y * weights) / sumweights
weightedcov <- cov.wt(dasim, wt = weights, cor = TRUE)
weightedsd <- sqrt(diag(weightedcov$cov))

results_2b["sx", "eqB"] <- weightedsd["x"]
results_2b["sw", "eqB"] <- weightedsd["w"]
results_2b["sz", "eqB"] <- weightedsd["z"]
results_2b["sy", "eqB"] <- weightedsd["y"]
results_2b["rxw", "eqB"] <- weightedcov$cor["x", "w"]
results_2b["rxz", "eqB"] <- weightedcov$cor["x", "z"]
results_2b["rxy", "eqB"] <- weightedcov$cor["x", "y"]
results_2b["rwz", "eqB"] <- weightedcov$cor["w", "z"]
results_2b["rwy", "eqB"] <- weightedcov$cor["w", "y"]
results_2b["rzy", "eqB"] <- weightedcov$cor["z", "y"]


# Standard multiple imputation (without response indicators)
im_2b <- mice(
  da[, c("z", "w", "x", "y")],
  method = c("norm", "norm", "norm", "norm"),
  printFlag = FALSE
)
meanx <- with(im_2b, mean(x))
meanw <- with(im_2b, mean(w))
meanz <- with(im_2b, mean(z))
meany <- with(im_2b, mean(y))
sdx <- with(im_2b, sd(x))
sdw <- with(im_2b, sd(w))
sdz <- with(im_2b, sd(z))
sdy <- with(im_2b, sd(y))
corxw <- with(im_2b, cor(x, w))
corxz <- with(im_2b, cor(x, z))
corxy <- with(im_2b, cor(x, y))
corwz <- with(im_2b, cor(w, z))
corwy <- with(im_2b, cor(w, y))
corzy <- with(im_2b, cor(z, y))

results_2b["x","mi"] <- mean(unlist(meanx$analyses))
results_2b["w","mi"] <- mean(unlist(meanw$analyses))
results_2b["z","mi"] <- mean(unlist(meanz$analyses))
results_2b["y","mi"] <- mean(unlist(meany$analyses))
results_2b["sx","mi"] <- mean(unlist(sdx$analyses))
results_2b["sw","mi"] <- mean(unlist(sdw$analyses))
results_2b["sz","mi"] <- mean(unlist(sdz$analyses))
results_2b["sy","mi"] <- mean(unlist(sdy$analyses))
results_2b["rxw","mi"] <- mean(unlist(corxw$analyses))
results_2b["rxz","mi"] <- mean(unlist(corxz$analyses))
results_2b["rxy","mi"] <- mean(unlist(corxy$analyses))
results_2b["rwz","mi"] <- mean(unlist(corwz$analyses))
results_2b["rwy","mi"] <- mean(unlist(corwy$analyses))
results_2b["rzy","mi"] <- mean(unlist(corzy$analyses))



# Multiple imputation with response indicators (MIRI)
im_2b <- mice(
  da,
  method = c("norm", "norm", "norm", "norm", "mean", "mean", "mean", "mean"),
  formulas = list(
    x = x ~ (z + w + y) * rz * rw * ry, w = w ~ (z + x + y) * rz * rx * ry,
    z = z ~ (x + w + y) * rx * rw * ry, y = y ~ (z + w + x) * rz * rw * rx,
    rz = rz ~ 1, rw = rw ~ 1, rx = rx ~ 1, ry = ry ~ 1
  ),
  visitSequence = c("x", "w", "z", "y", "rz", "rw", "rx", "ry"),
  printFlag = FALSE
)
meanx <- with(im_2b, mean(x))
meanw <- with(im_2b, mean(w))
meanz <- with(im_2b, mean(z))
meany <- with(im_2b, mean(y))
sdx <- with(im_2b, sd(x))
sdw <- with(im_2b, sd(w))
sdz <- with(im_2b, sd(z))
sdy <- with(im_2b, sd(y))
corxw <- with(im_2b, cor(x, w))
corxz <- with(im_2b, cor(x, z))
corxy <- with(im_2b, cor(x, y))
corwz <- with(im_2b, cor(w, z))
corwy <- with(im_2b, cor(w, y))
corzy <- with(im_2b, cor(z, y))

results_2b["x", "miri"] <- mean(unlist(meanx$analyses))
results_2b["w", "miri"] <- mean(unlist(meanw$analyses))
results_2b["z", "miri"] <- mean(unlist(meanz$analyses))
results_2b["y", "miri"] <- mean(unlist(meany$analyses))
results_2b["sx", "miri"] <- mean(unlist(sdx$analyses))
results_2b["sw", "miri"] <- mean(unlist(sdw$analyses))
results_2b["sz", "miri"] <- mean(unlist(sdz$analyses))
results_2b["sy", "miri"] <- mean(unlist(sdy$analyses))
results_2b["rxw", "miri"] <- mean(unlist(corxw$analyses))
results_2b["rxz", "miri"] <- mean(unlist(corxz$analyses))
results_2b["rxy", "miri"] <- mean(unlist(corxy$analyses))
results_2b["rwz", "miri"] <- mean(unlist(corwz$analyses))
results_2b["rwy", "miri"] <- mean(unlist(corwy$analyses))
results_2b["rzy", "miri"] <- mean(unlist(corzy$analyses))


# Factorizable multiple imputation (FMI)
daf <- da
daf$w[daf$rz == 0] <- NA
daf$x[daf$rz * daf$rw == 0] <- NA
daf$y[daf$rz * daf$rw * daf$rx == 0] <- NA
daf <- daf[,c("z","w","x","y")]

im_2b <- mice(
  daf,
  method = c("sample", "norm", "norm", "norm"),
  formulas = list(z = z ~ 1, w = w ~ z, x = x ~ w, y = y ~ x + z),
  visitSequence = c("z", "w", "x", "y"),
  maxit = 1,
  printFlag = FALSE
)
meanx <- with(im_2b, mean(x))
meanw <- with(im_2b, mean(w))
meanz <- with(im_2b, mean(z))
meany <- with(im_2b, mean(y))
sdx <- with(im_2b, sd(x))
sdw <- with(im_2b, sd(w))
sdz <- with(im_2b, sd(z))
sdy <- with(im_2b, sd(y))
corxw <- with(im_2b, cor(x, w))
corxz <- with(im_2b, cor(x, z))
corxy <- with(im_2b, cor(x, y))
corwz <- with(im_2b, cor(w, z))
corwy <- with(im_2b, cor(w, y))
corzy <- with(im_2b, cor(z, y))

results_2b["x","fmi"] <- mean(unlist(meanx$analyses))
results_2b["w","fmi"] <- mean(unlist(meanw$analyses))
results_2b["z","fmi"] <- mean(unlist(meanz$analyses))
results_2b["y","fmi"] <- mean(unlist(meany$analyses))
results_2b["sx","fmi"] <- mean(unlist(sdx$analyses))
results_2b["sw","fmi"] <- mean(unlist(sdw$analyses))
results_2b["sz","fmi"] <- mean(unlist(sdz$analyses))
results_2b["sy","fmi"] <- mean(unlist(sdy$analyses))
results_2b["rxw","fmi"] <- mean(unlist(corxw$analyses))
results_2b["rxz","fmi"] <- mean(unlist(corxz$analyses))
results_2b["rxy","fmi"] <- mean(unlist(corxy$analyses))
results_2b["rwz","fmi"] <- mean(unlist(corwz$analyses))
results_2b["rwy","fmi"] <- mean(unlist(corwy$analyses))
results_2b["rzy","fmi"] <- mean(unlist(corzy$analyses))

# Complete case analysis (CCA)
dacc <- da[rx & rw & rz & ry, ]
results_2b["x", "cca"] <- mean(dacc$x)
results_2b["w", "cca"] <- mean(dacc$w)
results_2b["z", "cca"] <- mean(dacc$z)
results_2b["y", "cca"] <- mean(dacc$y)
results_2b["sx", "cca"] <- sd(dacc$x)
results_2b["sw", "cca"] <- sd(dacc$w)
results_2b["sz", "cca"] <- sd(dacc$z)
results_2b["sy", "cca"] <- sd(dacc$y)
results_2b["rxw", "cca"] <- cor(dacc$x, dacc$w)
results_2b["rxz", "cca"] <- cor(dacc$x, dacc$z)
results_2b["rxy", "cca"] <- cor(dacc$x, dacc$y)
results_2b["rwz", "cca"] <- cor(dacc$w, dacc$z)
results_2b["rwy", "cca"] <- cor(dacc$w, dacc$y)
results_2b["rzy", "cca"] <- cor(dacc$z, dacc$y)

# Available case analysis
results_2b["x", "aca"] <- mean(da$x, na.rm = TRUE)
results_2b["w", "aca"] <- mean(da$w, na.rm = TRUE)
results_2b["z", "aca"] <- mean(da$z, na.rm = TRUE)
results_2b["y", "aca"] <- mean(da$y, na.rm = TRUE)
results_2b["sx", "aca"] <- sd(da$x, na.rm = TRUE)
results_2b["sw", "aca"] <- sd(da$w, na.rm = TRUE)
results_2b["sz", "aca"] <- sd(da$z, na.rm = TRUE)
results_2b["sy", "aca"] <- sd(da$y, na.rm = TRUE)
results_2b["rxw", "aca"] <- cor(da$x, da$w, use = "pairwise.complete.obs")
results_2b["rxz", "aca"] <- cor(da$x, da$z, use = "pairwise.complete.obs")
results_2b["rxy", "aca"] <- cor(da$x, da$y, use = "pairwise.complete.obs")
results_2b["rwz", "aca"] <- cor(da$w, da$z, use = "pairwise.complete.obs")
results_2b["rwy", "aca"] <- cor(da$w, da$y, use = "pairwise.complete.obs")
results_2b["rzy", "aca"] <- cor(da$z, da$y, use = "pairwise.complete.obs")

# Reporting
cat("\nExample 2b:\n")
tab_2b <- results_2b
tab_2b$truth <- formatC(tab_2b$truth, digits = 2, format = "f")
tab_2b$fmi <- formatC(tab_2b$fmi, digits = 2, format = "f")
tab_2b$mi <- formatC(tab_2b$mi, digits = 2, format = "f")
tab_2b$miri <- formatC(tab_2b$miri, digits = 2, format = "f")
tab_2b$cca <- formatC(tab_2b$cca, digits = 2, format = "f")
tab_2b$aca <- formatC(tab_2b$aca, digits = 2, format = "f")
tab_2b$eqA <- formatC(tab_2b$eqA, digits = 2, format = "f")
tab_2b$eqB <- formatC(tab_2b$eqB, digits = 2, format = "f")

tab_2b$truth <- paste0("$", formatC(tab_2b$truth, digits = 2, format = "f"), "$")
tab_2b$fmi <- paste0("$", formatC(tab_2b$fmi, digits = 2, format = "f"), "$")
tab_2b$mi <- paste0("$", formatC(tab_2b$mi, digits = 2, format = "f"), "$")
tab_2b$cca <- paste0("$", formatC(tab_2b$cca, digits = 2, format = "f"), "$")
tab_2b$aca <- paste0("$", formatC(tab_2b$aca, digits = 2, format = "f"), "$")
tab_2b$eqA <- paste0("$", formatC(tab_2b$eqA, digits = 2, format = "f"), "$")
tab_2b$eqB <- paste0("$", formatC(tab_2b$eqB, digits = 2, format = "f"), "$")

tab_2bb <- results_2b
tab_2bb$fmi <- formatC(tab_2bb$fmi - tab_2bb$truth, digits = 2, format = "f")
tab_2bb$mi <- formatC(tab_2bb$mi- tab_2bb$truth, digits = 2, format = "f")
tab_2bb$miri <- formatC(tab_2bb$miri - tab_2bb$truth, digits = 2, format = "f")
tab_2bb$cca <- formatC(tab_2bb$cca - tab_2bb$truth, digits = 2, format = "f")
tab_2bb$aca <- formatC(tab_2bb$aca - tab_2bb$truth, digits = 2, format = "f")
tab_2bb$eqA <- formatC(tab_2bb$eqA - tab_2bb$truth, digits = 2, format = "f")
tab_2bb$eqB <- formatC(tab_2bb$eqB - tab_2bb$truth, digits = 2, format = "f")
tab_2bb$truth <- formatC(tab_2bb$truth, digits = 2, format = "f")

tab_2bb$truth <- paste0("$", formatC(tab_2bb$truth, digits = 2, format = "f"), "$")
tab_2bb$fmi <- paste0("$", formatC(tab_2bb$fmi, digits = 2, format = "f"), "$")
tab_2bb$mi <- paste0("$", formatC(tab_2bb$mi, digits = 2, format = "f"), "$")
tab_2bb$cca <- paste0("$", formatC(tab_2bb$cca, digits = 2, format = "f"), "$")
tab_2bb$aca <- paste0("$", formatC(tab_2bb$aca, digits = 2, format = "f"), "$")
tab_2bb$eqA <- paste0("$", formatC(tab_2bb$eqA, digits = 2, format = "f"), "$")
tab_2bb$eqB <- paste0("$", formatC(tab_2bb$eqB, digits = 2, format = "f"), "$")

cat(reporderb,"\n")
cat(latextable(tab_2b[, reporderb]))
cat(reporderb,"\n")
cat(latextable(tab_2bb[, reporderb]))

#save.image(file = "simulation_results_fig2_20251021.Rdata")
