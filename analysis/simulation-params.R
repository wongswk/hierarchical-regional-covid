# Values of parameters for simulation study

sampFrac_1 <- c(0.06,0.39,0.07,0.02,0.11)
sampFrac_2 <- c(0.1,0.42,0.07,0.15,0.22)
sampFrac_3 <- c(0.27,0.38,0.22,0.24,0.26)
sampFrac_4 <- c(0.49,0.61,0.39,0.65,0.55)+0.05
f1 <- rep(1,5)
f2 <- c(0.13,0.3,0.29,0.17,0.39)+0.03
f3 <- c(0.79,0.72,0.66,0.95,0.59)
f4 <- c(0.62,0.66,0.67,0.52,0.63)
f5 <- c(0.45,0.46,0.39,0.68,0.63)
f6 <- c(0.99,0.79,0.87,0.79,0.75)
f7 <- c(0.46,0.46,0.61,0.59,0.49)+0.03
theta <- c(8,8,5,3,11)
R0b <- 3
pop_size <- c(843375,1225195,297570,795116,1889225)

show(rbind(R0b,  f2, f3, f4, f5, f6, f7, sampFrac_1, sampFrac_2, sampFrac_3, sampFrac_4, theta))
allpars <- rbind(R0b,  f2, f3, f4, f5, f6, f7, sampFrac_1, sampFrac_2, sampFrac_3, sampFrac_4, theta)
library(xtable)
xtable(allpars[,c(2,5,4,1,3)])
