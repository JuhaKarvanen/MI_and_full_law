# File: examples_mice.R
# Author: Juha Karvanen
# Date: 2024-10-24
# Summary: R code to reproduce the simulation results for Examples 1, 2, 3 and 4 
# of the paperJ. Karvanen, S. Tikka (2024), Multiple imputation and full law identifiability.

library(mice, quietly = TRUE, warn.conflicts = FALSE)
expit <- function(x) exp(x)/(1+exp(x))
n <- 1000000


results_1 <- data.frame(
  statistic = c("$\\E(X)$","$\\E(Y)$","$\\sd(X)$","$\\sd(Y)$","$\\Cor(X,Y)$"),
  truth = c(0,0,1,sqrt(2),sqrt(2)/2), mi = rep(NA,5), midecomp = rep(NA,5),  
            cca = rep(NA,5), aca = rep(NA,5), dosearch = rep(NA,5))
rownames(results_1) <- c("x","y","sx","sy","r")
results_2 <- results_1
results_3 <- results_1


############################################################################
# Example 1
############################################################################

set.seed(101020241)
x1 <- rnorm(n,0,1)
y1 <- x1 + rnorm(n,0,1) # X->Y
rx <- as.numeric( 0.7 > runif(n,0,1) ) 
ry <- as.numeric( expit(x1) > runif(n,0,1)) # X->RY
x <- x1
x[!rx] <- NA
y <- y1
y[!ry] <- NA
da <- data.frame(x = x, y = y, rx = rx, ry = ry)

# Standard multiple imputation
predmatrix <- make.predictorMatrix(da)
predmatrix["x","rx"] <- 0
predmatrix["y","ry"] <- 0
im1 <- mice(da, method = c("norm","norm","mean","mean"),
            predictorMatrix = predmatrix,
            visitSequence = c("x","y","rx","ry"), maxit = 10,
            printFlag = FALSE)
meanx <- with(im1,mean(x))
meany <- with(im1,mean(y))
sdx <- with(im1,sd(x))
sdy <- with(im1,sd(y))
cor <- with(im1,cor(x,y))

results_1["x","mi"] <- mean(unlist(meanx$analyses))
results_1["y","mi"] <- mean(unlist(meany$analyses))
results_1["sx","mi"] <- mean(unlist(sdx$analyses))
results_1["sy","mi"] <- mean(unlist(sdy$analyses))
results_1["r","mi"] <- mean(unlist(cor$analyses))

# Decomposable multiple imputation
dadecomp <- da
dadecomp$y[rx == 0 & ry == 1] <- NA 
im1decomp <- mice(dadecomp[,c("x","y")], method = c("sample","norm"),
                visitSequence = c("x","y"), maxit = 10,
                printFlag = FALSE)
meanx <- with(im1decomp,mean(x))
meany <- with(im1decomp,mean(y))
sdx <- with(im1decomp,sd(x))
sdy <- with(im1decomp,sd(y))
cor <- with(im1decomp,cor(x,y))

results_1["x","midecomp"] <- mean(unlist(meanx$analyses))
results_1["y","midecomp"] <- mean(unlist(meany$analyses))
results_1["sx","midecomp"] <- mean(unlist(sdx$analyses))
results_1["sy","midecomp"] <- mean(unlist(sdy$analyses))
results_1["r","midecomp"] <- mean(unlist(cor$analyses))

# Complete case analysis
dacc <- da[rx & ry,]
results_1["x","cca"] <- mean(dacc$x)
results_1["y","cca"] <- mean(dacc$y)
results_1["sx","cca"] <- sd(dacc$x)
results_1["sy","cca"] <- sd(dacc$y)
results_1["r","cca"] <- cor(dacc$x,dacc$y)

# Available case analysis
results_1["x","aca"] <- mean(da$x, na.rm = TRUE)
results_1["y","aca"] <- mean(da$y, na.rm = TRUE)
results_1["sx","aca"] <- sd(da$x, na.rm = TRUE)
results_1["sy","aca"] <- sd(da$y, na.rm = TRUE)
results_1["r","aca"] <- cor(da$x,da$y, use = "pairwise.complete.obs")

# Estimation with identifying formula
dax <- da[da$rx == 1,]
results_1["x","dosearch"] <- mean(dax$x)
results_1["sx","dosearch"] <- sd(dax$x)
yxmodel <- lm(y ~ x, data = dacc)
yxsigma <- summary(yxmodel)$sigma
ypred <- predict(yxmodel, newdata = dax)
ysim <- ypred + rnorm(length(ypred),0,yxsigma)
results_1["y","dosearch"] <- mean(ysim)
results_1["sy","dosearch"] <- sd(ysim)
results_1["r","dosearch"] <- cor(dax$x,ysim)

# Reporting
latextable<-function(td) {
  tddim<-dim(td);
  lt<-"\n";
  for (i in 1:tddim[1]) {
    for (j in 1:tddim[2]) {
      if (j==1) lt <- paste(lt, td[i,j])
      else lt<-paste(lt," & ",td[i,j])
    }
    if (i==tddim[1]) lt <- paste(lt, "\n")
    else lt<-paste(lt, "\\\\ \n")
  }
  cat(lt)
} 

tab_1 <- results_1
# tab_1$truth <- formatC(tab_1$truth, digits = 2, format = "f")
# tab_1$mi <- formatC(tab_1$mi, digits = 2, format = "f")
# tab_1$midecomp <- formatC(tab_1$midecomp, digits = 2, format = "f")
# tab_1$cca <- formatC(tab_1$cca, digits = 2, format = "f")
# tab_1$aca <- formatC(tab_1$aca, digits = 2, format = "f")
# tab_1$dosearch <- formatC(tab_1$dosearch, digits = 2, format = "f")
tab_1$truth <- paste0("$", formatC(tab_1$truth, digits = 2, format = "f"), "$")
tab_1$mi <- paste0("$", formatC(tab_1$mi, digits = 2, format = "f"), "$")
tab_1$midecomp <- paste0("$", formatC(tab_1$midecomp, digits = 2, format = "f"), "$")
tab_1$cca <- paste0("$", formatC(tab_1$cca, digits = 2, format = "f"), "$")
tab_1$aca <- paste0("$", formatC(tab_1$aca, digits = 2, format = "f"), "$")
tab_1$dosearch <- paste0("$", formatC(tab_1$dosearch, digits = 2, format = "f"), "$")

cat("\nExample 1:\n")
cat(latextable(tab_1[,-4]))



############################################################################
# Example 2
############################################################################

set.seed(101020242)
x1 <- rnorm(n,0,1)
y1 <- x1 + rnorm(n,0,1) # X->Y
rx <- as.numeric( 0.7 > runif(n,0,1) ) 
ry <- as.numeric( expit(x1 + 2*rx - 1) > runif(n,0,1)) #RX->RY, X->RY
x <- x1
x[!rx] <- NA
y <- y1
y[!ry] <- NA
da <- data.frame(x = x, y = y, rx = rx, ry = ry)

# Standard multiple imputation
predmatrix <- make.predictorMatrix(da)
predmatrix["x","rx"] <- 0
predmatrix["y","ry"] <- 0
im2 <- mice(da, method = c("norm","norm","mean","mean"),
            predictorMatrix = predmatrix,
            visitSequence = c("x","y","rx","ry"), maxit = 10,
            printFlag = FALSE)
meanx <- with(im2,mean(x))
meany <- with(im2,mean(y))
sdx <- with(im2,sd(x))
sdy <- with(im2,sd(y))
cor <- with(im2,cor(x,y))

results_2["x","mi"] <- mean(unlist(meanx$analyses))
results_2["y","mi"] <- mean(unlist(meany$analyses))
results_2["sx","mi"] <- mean(unlist(sdx$analyses))
results_2["sy","mi"] <- mean(unlist(sdy$analyses))
results_2["r","mi"] <- mean(unlist(cor$analyses))

# Decomposable multiple imputation
dadecomp <- da
dadecomp$y[rx == 0 & ry == 1] <- NA 
im2decomp <- mice(dadecomp[,c("x","y")], method = c("sample","norm"),
                visitSequence = c("x","y"), maxit = 10,
                printFlag = FALSE)
meanx <- with(im2decomp,mean(x))
meany <- with(im2decomp,mean(y))
sdx <- with(im2decomp,sd(x))
sdy <- with(im2decomp,sd(y))
cor <- with(im2decomp,cor(x,y))

results_2["x","midecomp"] <- mean(unlist(meanx$analyses))
results_2["y","midecomp"] <- mean(unlist(meany$analyses))
results_2["sx","midecomp"] <- mean(unlist(sdx$analyses))
results_2["sy","midecomp"] <- mean(unlist(sdy$analyses))
results_2["r","midecomp"] <- mean(unlist(cor$analyses))

# Complete case analysis
dacc <- da[ da$rx & da$ry,]
results_2["x","cca"] <- mean(dacc$x)
results_2["y","cca"] <- mean(dacc$y)
results_2["sx","cca"] <- sd(dacc$x)
results_2["sy","cca"] <- sd(dacc$y)
results_2["r","cca"] <- cor(dacc$x,dacc$y)

# Available case analysis
results_2["x","aca"] <- mean(da$x, na.rm = TRUE)
results_2["y","aca"] <- mean(da$y, na.rm = TRUE)
results_2["sx","aca"] <- sd(da$x, na.rm = TRUE)
results_2["sy","aca"] <- sd(da$y, na.rm = TRUE)
results_2["r","aca"] <- cor(da$x,da$y, use = "pairwise.complete.obs")

# Estimation with identifying formula
dax <- da[da$rx == 1,]
results_2["x","dosearch"] <- mean(dax$x)
results_2["sx","dosearch"] <- sd(dax$x)
yxmodel <- lm(y ~ x, data = dacc)
yxsigma <- summary(yxmodel)$sigma
ypred <- predict(yxmodel, newdata = dax)
ysim <- ypred + rnorm(length(ypred),0,yxsigma)
results_2["y","dosearch"] <- mean(ysim)
results_2["sy","dosearch"] <- sd(ysim)
results_2["r","dosearch"] <- cor(dax$x,ysim)

# Reporting
tab_2 <- results_2
# tab_2$truth <- formatC(tab_2$truth, digits = 2, format = "f")
# tab_2$mi <- formatC(tab_2$mi, digits = 2, format = "f")
# tab_2$midecomp <- formatC(tab_2$midecomp, digits = 2, format = "f")
# tab_2$cca <- formatC(tab_2$cca, digits = 2, format = "f")
# tab_2$aca <- formatC(tab_2$aca, digits = 2, format = "f")
# tab_2$dosearch <- formatC(tab_2$dosearch, digits = 2, format = "f")
tab_2$truth <- paste0("$", formatC(tab_2$truth, digits = 2, format = "f"), "$")
tab_2$mi <- paste0("$", formatC(tab_2$mi, digits = 2, format = "f"), "$")
tab_2$midecomp <- paste0("$", formatC(tab_2$midecomp, digits = 2, format = "f"), "$")
tab_2$cca <- paste0("$", formatC(tab_2$cca, digits = 2, format = "f"), "$")
tab_2$aca <- paste0("$", formatC(tab_2$aca, digits = 2, format = "f"), "$")
tab_2$dosearch <- paste0("$", formatC(tab_2$dosearch, digits = 2, format = "f"), "$")

cat("\nExample 2:\n")
cat(latextable(tab_2[,-4]))
cat(latextable(tab_2[,1:4]))





############################################################################
# Example 3
############################################################################

set.seed(101020243)
x1 <- rnorm(n,0,1)
y1 <- x1 + rnorm(n,0,1) # X->Y
rx <- as.numeric( 0.7 > runif(n,0,1) ) 
ry <- as.numeric( expit(x1 * (2*rx - 1)) > runif(n,0,1)) #RX->RY, X->RY
x <- x1
x[!rx] <- NA
y <- y1
y[!ry] <- NA
da <- data.frame(x = x, y = y, rx = rx, ry = ry)

# Standard multiple imputation
predmatrix <- make.predictorMatrix(da)
predmatrix["x","rx"] <- 0
predmatrix["y","ry"] <- 0
im3 <- mice(da, method = c("norm","norm","mean","mean"),
            predictorMatrix = predmatrix,
            visitSequence = c("x","y","rx","ry"), maxit = 10,
            printFlag = FALSE)
meanx <- with(im3,mean(x))
meany <- with(im3,mean(y))
sdx <- with(im3,sd(x))
sdy <- with(im3,sd(y))
cor <- with(im3,cor(x,y))

results_3["x","mi"] <- mean(unlist(meanx$analyses))
results_3["y","mi"] <- mean(unlist(meany$analyses))
results_3["sx","mi"] <- mean(unlist(sdx$analyses))
results_3["sy","mi"] <- mean(unlist(sdy$analyses))
results_3["r","mi"] <- mean(unlist(cor$analyses))

# Decomposable multiple imputation
dadecomp <- da
dadecomp$y[rx == 0 & ry == 1] <- NA 
im3decomp <- mice(dadecomp[,c("x","y")], method = c("sample","norm"),
                visitSequence = c("x","y"), maxit = 10,
                printFlag = FALSE)
meanx <- with(im3decomp,mean(x))
meany <- with(im3decomp,mean(y))
sdx <- with(im3decomp,sd(x))
sdy <- with(im3decomp,sd(y))
cor <- with(im3decomp,cor(x,y))

results_3["x","midecomp"] <- mean(unlist(meanx$analyses))
results_3["y","midecomp"] <- mean(unlist(meany$analyses))
results_3["sx","midecomp"] <- mean(unlist(sdx$analyses))
results_3["sy","midecomp"] <- mean(unlist(sdy$analyses))
results_3["r","midecomp"] <- mean(unlist(cor$analyses))


# Complete case analysis
dacc <- da[ da$rx & da$ry,]
results_3["x","cca"] <- mean(dacc$x)
results_3["y","cca"] <- mean(dacc$y)
results_3["sx","cca"] <- sd(dacc$x)
results_3["sy","cca"] <- sd(dacc$y)
results_3["r","cca"] <- cor(dacc$x,dacc$y)

# Available case analysis
results_3["x","aca"] <- mean(da$x, na.rm = TRUE)
results_3["y","aca"] <- mean(da$y, na.rm = TRUE)
results_3["sx","aca"] <- sd(da$x, na.rm = TRUE)
results_3["sy","aca"] <- sd(da$y, na.rm = TRUE)
results_3["r","aca"] <- cor(da$x,da$y, use = "pairwise.complete.obs")

# Estimation with identifying formula
dax <- da[da$rx == 1,]
results_3["x","dosearch"] <- mean(dax$x)
results_3["sx","dosearch"] <- sd(dax$x)
yxmodel <- lm(y ~ x, data = dacc)
yxsigma <- summary(yxmodel)$sigma
ypred <- predict(yxmodel, newdata = dax)
ysim <- ypred + rnorm(length(ypred),0,yxsigma)
results_3["y","dosearch"] <- mean(ysim)
results_3["sy","dosearch"] <- sd(ysim)
results_3["r","dosearch"] <- cor(dax$x,ysim)


# Reporting
tab_3 <- results_3
# tab_3$truth <- formatC(tab_3$truth, digits = 2, format = "f")
# tab_3$mi <- formatC(tab_3$mi, digits = 2, format = "f")
# tab_3$midecomp <- formatC(tab_3$midecomp, digits = 2, format = "f")
# tab_3$cca <- formatC(tab_3$cca, digits = 2, format = "f")
# tab_3$aca <- formatC(tab_3$aca, digits = 2, format = "f")
# tab_3$dosearch <- formatC(tab_3$dosearch, digits = 2, format = "f")
tab_3$truth <- paste0("$", formatC(tab_3$truth, digits = 2, format = "f"), "$")
tab_3$mi <- paste0("$", formatC(tab_3$mi, digits = 2, format = "f"), "$")
tab_3$midecomp <- paste0("$", formatC(tab_3$midecomp, digits = 2, format = "f"), "$")
tab_3$cca <- paste0("$", formatC(tab_3$cca, digits = 2, format = "f"), "$")
tab_3$aca <- paste0("$", formatC(tab_3$aca, digits = 2, format = "f"), "$")
tab_3$dosearch <- paste0("$", formatC(tab_3$dosearch, digits = 2, format = "f"), "$")



cat("\nExample 3:\n")
cat(latextable(tab_3[,-4]))
cat(latextable(tab_3[,1:4]))


############################################################################
# Example 4
############################################################################

set.seed(270120251)
results_4 <- data.frame(
  statistic = c("$\\E(X)$","$\\E(W)$","$\\E(Z)$","$\\E(Y)$",
                "$\\sd(X)$","$\\sd(W)$","$\\sd(Z)$","$\\sd(Y)$",
                "$\\Cor(X,W)$","$\\Cor(X,Z)$","$\\Cor(X,Y)$",  
                "$\\Cor(W,Z)$","$\\Cor(W,Y)$","$\\Cor(Z,Y)$"
  ),
  truth = c(0,0,0,0,1,1,1,1,sqrt(2)/2,0.5,3/4,sqrt(2)/2,sqrt(2)/2,3/4), 
  midecomp = rep(NA,14), mi = rep(NA,14),   cca = rep(NA,14), aca = rep(NA,14))
rownames(results_4) <- c("x","w","z","y","sx","sw","sz","sy","rxw","rxz","rxy","rwz","rwy","rzy")


x1 <- rnorm(n,0,1)
w1 <- (x1 + rnorm(n,0,1))/sqrt(2)
z1 <- (w1 + rnorm(n,0,1))/sqrt(2)
y1 <- (x1 + z1 + rnorm(n,0,1))/sqrt(4) 
da1 <- data.frame(z = z1, w = w1, x = x1,  y = y1)
rx <- as.numeric( expit(w1) > runif(n,0,1))
rz <- as.numeric( 0.7 > runif(n,0,1) ) 
rw <- as.numeric( expit(z1+2*rz-1) > runif(n,0,1))
ry <- as.numeric( expit(x1) > runif(n,0,1)) 
x <- x1
x[!rx] <- NA
z <- z1
z[!rz] <- NA
w <- w1
w[!rw] <- NA
y <- y1
y[!ry] <- NA
da <- data.frame(x = x, w = w, z = z, y = y,  rx = rx, rw = rw, rz = rz, ry = ry)

# Standard multiple imputation
im4 <- mice(da[,c("z","w","x","y")], method = c("norm","norm","norm","norm"),
            printFlag = FALSE)
meanx <- with(im4,mean(x))
meanw <- with(im4,mean(w))
meanz <- with(im4,mean(z))
meany <- with(im4,mean(y))
sdx <- with(im4,sd(x))
sdw <- with(im4,sd(w))
sdz <- with(im4,sd(z))
sdy <- with(im4,sd(y))
corxw <- with(im4,cor(x,w))
corxz <- with(im4,cor(x,z))
corxy <- with(im4,cor(x,y))
corwz <- with(im4,cor(w,z))
corwy <- with(im4,cor(w,y))
corzy <- with(im4,cor(z,y))

results_4["x","mi"] <- mean(unlist(meanx$analyses))
results_4["w","mi"] <- mean(unlist(meanw$analyses))
results_4["z","mi"] <- mean(unlist(meanz$analyses))
results_4["y","mi"] <- mean(unlist(meany$analyses))
results_4["sx","mi"] <- mean(unlist(sdx$analyses))
results_4["sw","mi"] <- mean(unlist(sdw$analyses))
results_4["sz","mi"] <- mean(unlist(sdz$analyses))
results_4["sy","mi"] <- mean(unlist(sdy$analyses))
results_4["rxw","mi"] <- mean(unlist(corxw$analyses))
results_4["rxz","mi"] <- mean(unlist(corxz$analyses))
results_4["rxy","mi"] <- mean(unlist(corxy$analyses))
results_4["rwz","mi"] <- mean(unlist(corwz$analyses))
results_4["rwy","mi"] <- mean(unlist(corwy$analyses))
results_4["rzy","mi"] <- mean(unlist(corzy$analyses))


# Decomposable multiple imputation
dadecomp <- da
dadecomp$w[dadecomp$rz == 0] <- NA
dadecomp$x[dadecomp$rz * dadecomp$rw == 0] <- NA
dadecomp$y[dadecomp$rz * dadecomp$rw * dadecomp$rx == 0] <- NA
dadecomp <- dadecomp[,c("z","w","x","y")]

im4decomp <- mice(dadecomp, method = c("sample","norm","norm","norm"),
            formulas = list(z = z ~ 1, w = w ~ z, x = x ~ w, y = y ~ x + z),
            visitSequence = c("z","w","x","y"), maxit = 5,
            printFlag = FALSE)
meanx <- with(im4decomp ,mean(x))
meanw <- with(im4decomp ,mean(w))
meanz <- with(im4decomp ,mean(z))
meany <- with(im4decomp ,mean(y))
sdx <- with(im4decomp ,sd(x))
sdw <- with(im4decomp ,sd(w))
sdz <- with(im4decomp ,sd(z))
sdy <- with(im4decomp,sd(y))
corxw <- with(im4decomp,cor(x,w))
corxz <- with(im4decomp,cor(x,z))
corxy <- with(im4decomp,cor(x,y))
corwz <- with(im4decomp,cor(w,z))
corwy <- with(im4decomp,cor(w,y))
corzy <- with(im4decomp,cor(z,y))

results_4["x","midecomp"] <- mean(unlist(meanx$analyses))
results_4["w","midecomp"] <- mean(unlist(meanw$analyses))
results_4["z","midecomp"] <- mean(unlist(meanz$analyses))
results_4["y","midecomp"] <- mean(unlist(meany$analyses))
results_4["sx","midecomp"] <- mean(unlist(sdx$analyses))
results_4["sw","midecomp"] <- mean(unlist(sdw$analyses))
results_4["sz","midecomp"] <- mean(unlist(sdz$analyses))
results_4["sy","midecomp"] <- mean(unlist(sdy$analyses))
results_4["rxw","midecomp"] <- mean(unlist(corxw$analyses))
results_4["rxz","midecomp"] <- mean(unlist(corxz$analyses))
results_4["rxy","midecomp"] <- mean(unlist(corxy$analyses))
results_4["rwz","midecomp"] <- mean(unlist(corwz$analyses))
results_4["rwy","midecomp"] <- mean(unlist(corwy$analyses))
results_4["rzy","midecomp"] <- mean(unlist(corzy$analyses))


# Complete case analysis
dacc <- da[rx & rw & rz & ry,]
results_4["x","cca"] <- mean(dacc$x)
results_4["w","cca"] <- mean(dacc$w)
results_4["z","cca"] <- mean(dacc$z)
results_4["y","cca"] <- mean(dacc$y)
results_4["sx","cca"] <- sd(dacc$x)
results_4["sw","cca"] <- sd(dacc$w)
results_4["sz","cca"] <- sd(dacc$z)
results_4["sy","cca"] <- sd(dacc$y)
results_4["rxw","cca"] <- cor(dacc$x,dacc$w)
results_4["rxz","cca"] <- cor(dacc$x,dacc$z)
results_4["rxy","cca"] <- cor(dacc$x,dacc$y)
results_4["rwz","cca"] <- cor(dacc$w,dacc$z)
results_4["rwy","cca"] <- cor(dacc$w,dacc$y)
results_4["rzy","cca"] <- cor(dacc$z,dacc$y)

# Available case analysis
results_4["x","aca"] <- mean(da$x, na.rm = TRUE)
results_4["w","aca"] <- mean(da$w, na.rm = TRUE)
results_4["z","aca"] <- mean(da$z, na.rm = TRUE)
results_4["y","aca"] <- mean(da$y, na.rm = TRUE)
results_4["sx","aca"] <- sd(da$x, na.rm = TRUE)
results_4["sw","aca"] <- sd(da$w, na.rm = TRUE)
results_4["sz","aca"] <- sd(da$z, na.rm = TRUE)
results_4["sy","aca"] <- sd(da$y, na.rm = TRUE)
results_4["rxw","aca"] <- cor(da$x,da$w, use = "pairwise.complete.obs")
results_4["rxz","aca"] <- cor(da$x,da$z, use = "pairwise.complete.obs")
results_4["rxy","aca"] <- cor(da$x,da$y, use = "pairwise.complete.obs")
results_4["rwz","aca"] <- cor(da$w,da$z, use = "pairwise.complete.obs")
results_4["rwy","aca"] <- cor(da$w,da$y, use = "pairwise.complete.obs")
results_4["rzy","aca"] <- cor(da$z,da$y, use = "pairwise.complete.obs")

# Reporting
cat("\nExample 4:\n")
tab_4 <- results_4
tab_4$truth <- formatC(tab_4$truth, digits = 2, format = "f")
tab_4$midecomp <- formatC(tab_4$midecomp, digits = 2, format = "f")
tab_4$mi <- formatC(tab_4$mi, digits = 2, format = "f")
tab_4$cca <- formatC(tab_4$cca, digits = 2, format = "f")
tab_4$aca <- formatC(tab_4$aca, digits = 2, format = "f")

tab_4$truth <- paste0("$", formatC(tab_4$truth, digits = 2, format = "f"), "$")
tab_4$midecomp <- paste0("$", formatC(tab_4$midecomp, digits = 2, format = "f"), "$")
tab_4$mi <- paste0("$", formatC(tab_4$mi, digits = 2, format = "f"), "$")
tab_4$cca <- paste0("$", formatC(tab_4$cca, digits = 2, format = "f"), "$")
tab_4$aca <- paste0("$", formatC(tab_4$aca, digits = 2, format = "f"), "$")

cat(latextable(tab_4))

#save.image(file = "simulation_results1234.Rdata")