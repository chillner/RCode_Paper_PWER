###########################################################
#Single stage designs PWER-control (Sections 5.1 and 5.2)
###########################################################
library(multcomp)

##Correlation of the test statistics:
#one treatment:
corpi12_1 <- function(pi12){2*pi12/(1+pi12)}
#multiple treatments:
corpi12_2 <- function(pi12){
  3*pi12/(2*(1+2*pi12))  
}

#critical values for PWER- and FWER-control
mtcritpwer <- function(alpha, pi12, corr){
  Sigma <- matrix(c(1, corr, corr, 1), nr = 2, byrow = T)
  res <- function(crit){(1-pi12)*(1-pnorm(crit)) + pi12*(1-pmvnorm(upper=rep(crit,2), corr = Sigma)[1])-alpha}
  return(uniroot(res, interval = c(1,10))$root)
}
mtcritfwer <- function(alpha, pi12, corr){
  Sigma <- matrix(c(1, corr, corr, 1), nr = 2, byrow = T)
  res <- qmvnorm(1-alpha, tail = "lower.tail", sigma = Sigma)$quantile
  return(res)
}

#plot sample size increase for PWER/FWER vs. unadjusted case:
beta <- 0.2; alpha = 0.025
#grid for pi12:
pi12 <- seq(0, 1, 0.01)
factor_fwer <- factor_pwer <- numeric(length(pi12))

#Calculate sample size factors for each value in pi12
for(i in 1:length(pi12)){
  factor_fwer[i] <- (qnorm(1-beta) + mtcritfwer(alpha=alpha, pi12=pi12[i], corr = corpi12_2(pi12[i])))^2/(qnorm(1-beta)+qnorm(1-alpha))^2
  factor_pwer[i] <- (qnorm(1-beta) + mtcritpwer(alpha=alpha, pi12=pi12[i], corr = corpi12_2(pi12[i])))^2/(qnorm(1-beta)+qnorm(1-alpha))^2
}
#Plots
plot(pi12, factor_fwer, xlab = expression(pi[12]), ylim = c(1, 1.25), type = "l",
	ylab = "Factor of sample size increase for PWER/FWER", col = "red")
lines(pi12, factor_pwer, col = "blue")
text(x=0.1, y=1.215, labels = "FWER", col = "red")
text(x=0.1, y=1.045, labels = "PWER", col = "blue")

#plug in corpi12_1(pi12[i]) in the for-loop above for corr in mtcritpwer/fwer to get the one treatment case (T1 = T2) in Section 5.2
#plug in corpi12_2(pi12[i]) in the for-loop above for corr in mtcritpwer/fwer to get the two treatment case (T1 != T2) in Section 5.2
#plug in 0 in the for-loop above for corr in mtcritpwer/fwer to get the independent case in Section 5.1
