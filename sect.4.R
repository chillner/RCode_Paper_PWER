#################################################
##upper boundaries conditional FWER (Section 4)
#################################################
library(multcomp)
#m = 3 intersecting sets
#We test Hi: theta_i <= 0 and reject it if Z_i >= c; i in {1,...,m}
#crit: critical value
#corr: m x m correlation matrix
#poptype: 1 = not nested; 2 = nested
#p = (p{1}, p{2}, p{3}, p{1,2}, p{1,3}, p{2,3}, p{1,2,3}) for poptype == 1
#p = (p{1}, p{1,2}, p{1,2,3}) for poptype == 2
#J0 = {1,2,3} for both types

#correlations:
#not nested:
#Corr(Ti,Tj) = sum_{J: i,j in J}(pi_J)/sqrt(pi_i*pi_j)
#nested (special case of above case)
#Corr(Ti,Tj) = sqrt(pi[j]/pi[i]) with j>i where pi[1] = p{1}+p{1,2}+p{1,2,3} = 1; pi[2] = p{1,2}+p{1,2,3}; p[3] = p{1,2,3}

corrmat <- function(p, poptype = 1){
  V <- diag(3)
  if(poptype == 1){
    p_i <- c(p[1]+p[4]+p[5]+p[7], p[2]+p[4]+p[6]+p[7], p[3]+p[5]+p[6]+p[7])
    V[1,2] <- V[2,1] <- (p[4]+p[7])/sqrt(p_i[1]*p_i[2])
    V[1,3] <- V[3,1] <- (p[5]+p[7])/sqrt(p_i[1]*p_i[3])
    V[2,3] <- V[3,2] <- (p[6]+p[7])/sqrt(p_i[2]*p_i[3])
  } else if(poptype == 2){
    p_i <- sapply(1:3, function(i) sum(p[3:i]))
    V[1,2] <- V[2,1] <- sqrt(p_i[2]/p_i[1]) 
    V[1,3] <- V[3,1] <- sqrt(p_i[3]/p_i[1])
    V[2,3] <- V[3,2] <- sqrt(p_i[3]/p_i[2])
  }
  return(V)
}
#PWER
pwer <- function(crit, p, poptype = 1){
  V <- corrmat(p, poptype)
  if(poptype == 1){
    (p[1]+p[2]+p[3])*pnorm(crit, lower.tail = F)+
     p[4]*(1-pmvnorm(upper = rep(crit,2), corr = V[1:2,1:2])[1])+ 
     p[5]*(1-pmvnorm(upper = rep(crit,2), corr = V[c(1,3),c(1,3)])[1])+ 
     p[6]*(1-pmvnorm(upper = rep(crit,2), corr = V[c(2,3),c(2,3)])[1])+
     p[7]*(1-pmvnorm(upper = rep(crit,3), corr = V)[1])
  } else if(poptype == 2){
    p[1]*pnorm(crit, lower.tail = F)+
    p[2]*(1-pmvnorm(upper = rep(crit,2), corr = V[1:2,1:2])[1])+
    p[3]*(1-pmvnorm(upper = rep(crit,3), corr = V)[1])
  }
}
#PWER without PJ0
cpwerJ0 <- function(crit, p, poptype = 1){
  V <- corrmat(p, poptype)
  if(poptype == 1){
    ((p[1]+p[2]+p[3])*pnorm(crit, lower.tail = F)+
    p[4]*(1-pmvnorm(upper = rep(crit,2), corr = V[1:2,1:2])[1])+ 
    p[5]*(1-pmvnorm(upper = rep(crit,2), corr = V[c(1,3),c(1,3)])[1])+ 
    p[6]*(1-pmvnorm(upper = rep(crit,2), corr = V[c(2,3),c(2,3)])[1]))/(1-p[7]) 
  } else if(poptype == 2){
    (p[1]*pnorm(crit, lower.tail = F)+
    p[2]*(1-pmvnorm(upper = rep(crit,2), corr = V[1:2,1:2])[1]))/(1-p[3])
  }
}
#FWER cond. on J0
fwerJ0 <- function(crit, p, poptype = 1, J0=1:3){
  V <- corrmat(p, poptype=poptype)
  1-pmvnorm(upper=rep(crit,length(J0)), corr = V[J0,J0])[1]
}
#critical value/quantile 
critval <- function(p, poptype = 1, J0empty = F, alpha = .025){
  if(!J0empty){
    f <- function(x) cpwerJ0(crit = x, p = p, poptype = poptype) - alpha
  } else{
    f <- function(x) pwer(crit = x, p = p, poptype = poptype) - alpha
  }
  uniroot(f, interval = c(0, 5))$root 
}
#function to find upper bound and left side of inequality (6)
findUpperBound <- function(p, poptype = 1, alpha = .025){
  #critical value actual PWER:
  critpwer <- critval(p=p, poptype = poptype, J0empty = T, alpha = alpha)
  #critical value complementary PWER:
  critpwer_no_J0 <- critval(p=p, poptype = poptype, J0empty = F, alpha = alpha)
  #left side of inequality: fwer_I
  x1 <- fwerJ0(crit = critpwer, p=p, poptype = poptype, J0 = 1:3)
  #upper bound: 
  x2 <- fwerJ0(crit = critpwer_no_J0, p=p, poptype = poptype, J0 = 1:3)
  return(c(maxFWER = x1, upperBound = x2))
}
#################
#examples:
#################
#not nested (poptype=1)
p <- c(2,1,1,1,1,.5,.5)/7
findUpperBound(p=p)


#nested (poptype=2)
p = c(0.6, 0.35, 0.05)
findUpperBound(p=p, poptype = 2)
