# File: examples_mice.R
# Author: Juha Karvanen
# Date: 2024-10-24
# Summary: Simulation experiments of the paper
# J. Karvanen, S. Tikka (2024), Multiple imputation and full law identifiability.


library(mice, quietly = TRUE, warn.conflicts = FALSE)
set.seed(101020241)
expit <- function(x) exp(x)/(1+exp(x))
n <- 1000000

# Example (a)
results_a <- data.frame(#statistic = c("$\\bar{x}$","$\\bar{y}$","$\\sigma_x$","$\\sigma_y$","$\\rho_{x,y}$"),
                        statistic = c("$\\E(X)$","$\\E(Y)$","$\\sd(X)$","$\\sd(Y)$","$Cor{X,Y}$"),
                        truth = c(0,0,1,sqrt(2),sqrt(2)/2), mi = rep(NA,5), cca = rep(NA,5), aca = rep(NA,5), dosearch = rep(NA,5))
rownames(results_a) <- c("x","y","sx","sy","r")
results_b1 <- results_a
results_b2 <- results_a



x1 <- rnorm(n,0,1)
y1 <- x1 + rnorm(n,0,1) # X->Y
rx <- as.numeric( 0.7 > runif(n,0,1) ) 
ry <- as.numeric( expit(x1) > runif(n,0,1)) # X->RY
x <- x1
x[!rx] <- NA
y <- y1
y[!ry] <- NA
da <- data.frame(x = x, y = y, rx = rx, ry = ry)
md.pattern(da, plot = FALSE)

predmatrix <- make.predictorMatrix(da)
predmatrix["x","rx"] <- 0
predmatrix["y","ry"] <- 0
ima <- mice(da, method = c("norm","norm","mean","mean"),
            predictorMatrix = predmatrix,
            visitSequence = c("x","y","rx","ry"), maxit = 10,
            printFlag = FALSE)
a_meanx <- with(ima,mean(x))
a_meany <- with(ima,mean(y))
a_sdx <- with(ima,sd(x))
a_sdy <- with(ima,sd(y))
a_cor <- with(ima,cor(x,y))

results_a["x","mi"] <- mean(unlist(a_meanx$analyses))
results_a["y","mi"] <- mean(unlist(a_meany$analyses))
results_a["sx","mi"] <- mean(unlist(a_sdx$analyses))
results_a["sy","mi"] <- mean(unlist(a_sdy$analyses))
results_a["r","mi"] <- mean(unlist(a_cor$analyses))


dacc <- da[rx & ry,]
results_a["x","cca"] <- mean(dacc$x)
results_a["y","cca"] <- mean(dacc$y)
results_a["sx","cca"] <- sd(dacc$x)
results_a["sy","cca"] <- sd(dacc$y)
results_a["r","cca"] <- cor(dacc$x,dacc$y)

results_a["x","aca"] <- mean(da$x, na.rm = TRUE)
results_a["y","aca"] <- mean(da$y, na.rm = TRUE)
results_a["sx","aca"] <- sd(da$x, na.rm = TRUE)
results_a["sy","aca"] <- sd(da$y, na.rm = TRUE)
results_a["r","aca"] <- cor(da$x,da$y, use = "pairwise.complete.obs")

dax <- da[da$rx == 1,]
results_a["x","dosearch"] <- mean(dax$x)
results_a["sx","dosearch"] <- sd(dax$x)
yxmodel <- lm(y ~ x, data = dacc)
yxsigma <- summary(yxmodel)$sigma
ypred <- predict(yxmodel, newdata = dax)
ysim <- ypred + rnorm(length(ypred),0,yxsigma)
results_a["y","dosearch"] <- mean(ysim)
results_a["sy","dosearch"] <- sd(ysim)
results_a["r","dosearch"] <- cor(dax$x,ysim)


# Example (b1)
x1 <- rnorm(n,0,1)
y1 <- x1 + rnorm(n,0,1) # X->Y
rx <- as.numeric( 0.7 > runif(n,0,1) ) 
 ry <- as.numeric( expit(x1 + 2*rx - 1) > runif(n,0,1)) #RX->RY, X->RY
x <- x1
x[!rx] <- NA
y <- y1
y[!ry] <- NA
db <- data.frame(x = x, y = y, rx = rx, ry = ry)
md.pattern(db, plot = FALSE)

predmatrix <- make.predictorMatrix(db)
predmatrix["x","rx"] <- 0
predmatrix["y","ry"] <- 0
imb <- mice(db, method = c("norm","norm","mean","mean"),
            predictorMatrix = predmatrix,
            visitSequence = c("x","y","rx","ry"), maxit = 10,
            printFlag = FALSE)
b_meanx <- with(imb,mean(x))
b_meany <- with(imb,mean(y))
b_sdx <- with(imb,sd(x))
b_sdy <- with(imb,sd(y))
b_cor <- with(imb,cor(x,y))

results_b1["x","mi"] <- mean(unlist(b_meanx$analyses))
results_b1["y","mi"] <- mean(unlist(b_meany$analyses))
results_b1["sx","mi"] <- mean(unlist(b_sdx$analyses))
results_b1["sy","mi"] <- mean(unlist(b_sdy$analyses))
results_b1["r","mi"] <- mean(unlist(b_cor$analyses))


dbcc <- db[ db$rx & db$ry,]
results_b1["x","cca"] <- mean(dbcc$x)
results_b1["y","cca"] <- mean(dbcc$y)
results_b1["sx","cca"] <- sd(dbcc$x)
results_b1["sy","cca"] <- sd(dbcc$y)
results_b1["r","cca"] <- cor(dbcc$x,dbcc$y)

results_b1["x","aca"] <- mean(db$x, na.rm = TRUE)
results_b1["y","aca"] <- mean(db$y, na.rm = TRUE)
results_b1["sx","aca"] <- sd(db$x, na.rm = TRUE)
results_b1["sy","aca"] <- sd(db$y, na.rm = TRUE)
results_b1["r","aca"] <- cor(db$x,db$y, use = "pairwise.complete.obs")

dbx <- db[db$rx == 1,]
results_b1["x","dosearch"] <- mean(dbx$x)
results_b1["sx","dosearch"] <- sd(dbx$x)
yxmodel <- lm(y ~ x, data = dbcc)
yxsigma <- summary(yxmodel)$sigma
ypred <- predict(yxmodel, newdata = dbx)
ysim <- ypred + rnorm(length(ypred),0,yxsigma)
results_b1["y","dosearch"] <- mean(ysim)
results_b1["sy","dosearch"] <- sd(ysim)
results_b1["r","dosearch"] <- cor(dbx$x,ysim)


# Example (b2)
x1 <- rnorm(n,0,1)
y1 <- x1 + rnorm(n,0,1) # X->Y
rx <- as.numeric( 0.7 > runif(n,0,1) ) 
ry <- as.numeric( expit(x1 * (2*rx - 1)) > runif(n,0,1)) #RX->RY, X->RY
x <- x1
x[!rx] <- NA
y <- y1
y[!ry] <- NA
db <- data.frame(x = x, y = y, rx = rx, ry = ry)
md.pattern(db, plot = FALSE)

predmatrix <- make.predictorMatrix(db)
predmatrix["x","rx"] <- 0
predmatrix["y","ry"] <- 0
imb <- mice(db, method = c("norm","norm","mean","mean"),
            predictorMatrix = predmatrix,
            visitSequence = c("x","y","rx","ry"), maxit = 10,
            printFlag = FALSE)
b_meanx <- with(imb,mean(x))
b_meany <- with(imb,mean(y))
b_sdx <- with(imb,sd(x))
b_sdy <- with(imb,sd(y))
b_cor <- with(imb,cor(x,y))

results_b2["x","mi"] <- mean(unlist(b_meanx$analyses))
results_b2["y","mi"] <- mean(unlist(b_meany$analyses))
results_b2["sx","mi"] <- mean(unlist(b_sdx$analyses))
results_b2["sy","mi"] <- mean(unlist(b_sdy$analyses))
results_b2["r","mi"] <- mean(unlist(b_cor$analyses))


dbcc <- db[ db$rx & db$ry,]
results_b2["x","cca"] <- mean(dbcc$x)
results_b2["y","cca"] <- mean(dbcc$y)
results_b2["sx","cca"] <- sd(dbcc$x)
results_b2["sy","cca"] <- sd(dbcc$y)
results_b2["r","cca"] <- cor(dbcc$x,dbcc$y)

results_b2["x","aca"] <- mean(db$x, na.rm = TRUE)
results_b2["y","aca"] <- mean(db$y, na.rm = TRUE)
results_b2["sx","aca"] <- sd(db$x, na.rm = TRUE)
results_b2["sy","aca"] <- sd(db$y, na.rm = TRUE)
results_b2["r","aca"] <- cor(db$x,db$y, use = "pairwise.complete.obs")

dbx <- db[db$rx == 1,]
results_b2["x","dosearch"] <- mean(dbx$x)
results_b2["sx","dosearch"] <- sd(dbx$x)
yxmodel <- lm(y ~ x, data = dbcc)
yxsigma <- summary(yxmodel)$sigma
ypred <- predict(yxmodel, newdata = dbx)
ysim <- ypred + rnorm(length(ypred),0,yxsigma)
results_b2["y","dosearch"] <- mean(ysim)
results_b2["sy","dosearch"] <- sd(ysim)
results_b2["r","dosearch"] <- cor(dbx$x,ysim)

#save.image(file = "simulation_results.Rdata")

# Reporting (creates Table 1)

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

tab_a <- results_a
tab_a$truth <- formatC(tab_a$truth, digits = 2, format = "f")
tab_a$mi <- formatC(tab_a$mi, digits = 2, format = "f")
tab_a$cca <- formatC(tab_a$cca, digits = 2, format = "f")
tab_a$aca <- formatC(tab_a$aca, digits = 2, format = "f")
tab_a$dosearch <- formatC(tab_a$dosearch, digits = 2, format = "f")

cat(latextable(tab_a))

tab_b <- results_b1
tab_b$truth <- formatC(tab_b$truth, digits = 2, format = "f")
tab_b$mi <- formatC(tab_b$mi, digits = 2, format = "f")
tab_b$cca <- formatC(tab_b$cca, digits = 2, format = "f")
tab_b$aca <- formatC(tab_b$aca, digits = 2, format = "f")
tab_b$dosearch <- formatC(tab_b$dosearch, digits = 2, format = "f")

cat(latextable(tab_b))


tab_b2 <- results_b2
tab_b2$truth <- formatC(tab_b2$truth, digits = 2, format = "f")
tab_b2$mi <- formatC(tab_b2$mi, digits = 2, format = "f")
tab_b2$cca <- formatC(tab_b2$cca, digits = 2, format = "f")
tab_b2$aca <- formatC(tab_b2$aca, digits = 2, format = "f")
tab_b2$dosearch <- formatC(tab_b2$dosearch, digits = 2, format = "f")

cat(latextable(tab_b2))
