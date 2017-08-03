install.packages("INLA", repos="https://www.math.ntnu.no/inla/R/stable")
install.packages("sp","MatrixModels")
install.packages("devtools")
library(INLA)

library(devtools)
install_github("gkeele/Diploffect.INLA")
library(Diploffect.INLA)
