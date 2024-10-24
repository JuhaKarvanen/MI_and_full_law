# File: examples_dosearch.R
# Author: Juha Karvanen
# Date: 2024-10-24
# Summary: Derives the identifying functionals for the examples of the paper
# J. Karvanen, S. Tikka (2024), Multiple imputation and full law identifiability.

library(dosearch)

data <- "
  p(X*,Y*,R_X,R_Y)
"
targetlaw <- "p(X,Y)"
fulllaw <- "p(X,Y,R_X,R_Y)"
imputeX <- "p(X | Y, R_X, R_Y)"
imputeY <- "p(Y | X, R_X, R_Y)"
imputeXY <- "p(X, Y | R_X, R_Y)"


graphA <- "
  X -> Y
  X -> R_Y
"
graphB <- "
  X -> Y
  X -> R_Y
  R_X -> R_Y
"

# Example (a)
dosearch(data, targetlaw, graphA, missing_data = "R_X : X, R_Y : Y")
dosearch(data, fulllaw, graphA, missing_data = "R_X : X, R_Y : Y")
dosearch(data, imputeX, graphA, missing_data = "R_X : X, R_Y : Y")
dosearch(data, imputeY, graphA, missing_data = "R_X : X, R_Y : Y")
dosearch(data, imputeXY, graphA, missing_data = "R_X : X, R_Y : Y")

# Example (b)
dosearch(data, targetlaw, graphB, missing_data = "R_X : X, R_Y : Y")
dosearch(data, fulllaw, graphB, missing_data = "R_X : X, R_Y : Y")


