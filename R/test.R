
N = 7369+70180
R2xz = 0.21/100
k = 3


f.value <- 1 + N * R2xz / (1 - R2xz)
f.value2 <- (R2xz / (1-R2xz))*((N-1-k) /k)


print(f.value)
print(f.value2)







getBinaryPC <- function(rsq, pval, ratio, N, beta) {

  # Logarithms
  beta = log(beta)

  if (sign(beta) == -1) {
    power = pnorm(sqrt(N*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*beta-qnorm(1-pval/2, lower.tail = F), lower.tail = F)
  } else {
    power = pnorm(sqrt(N*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*beta-qnorm(1-pval/2))
  }
  return (power)
}


getContinuousPC <- function(rsq, sig, n, b1, power) {

  if (sign(b1) == -1) {

    res = pnorm(sqrt(n*rsq)*b1-qnorm(1-sig/2, lower.tail = F), lower.tail = F)
  } else {
    res = pnorm(sqrt(n*rsq)*b1-qnorm(1-sig/2))
  }
  return (res)
}

continuousPC(rsq=0.02, sig=0.05, n=1050, b1=0.34)




changeOR <- function(OR) {
  if(OR < 1) {
    OR = 1/OR
  }
  return(OR)
}


xaix1 = seq(from=0.1,to=1.0, by = 0.03)
xaix2 = unlist(lapply(xaix1, changeOR))
xaix3 = c(xaix1, xaix2)

a = NULL
for (i in xaix3) {

  res = (binaryPC(rsq=rsq, sig=sig, ratio=ratio, n=n, b1=i))
  a = c(a, res)
}

plot(xaix3, a, ylim= c(0.0, 1.0))



(binaryPC(rsq=0.21, sig=0.05, ratio=0.1764706, n=1050, b1=1.7)) *100


binaryPCOld <- function(rsq, sig, K, n, b1) {

  b_MR <- K * ( b1/ (1 + K * (b1 - 1)) -1)
  v_MR <- (K * (1-K) - b_MR^2) / (n*rsq)
  NCP <- b_MR^2 / v_MR

  # 2-sided test
  threschi <- qchisq(1 - sig, 1) # threshold chi(1) scale
  return(1 - pchisq(threschi, 1, NCP))
}



rsq = 0.30 # squared correlation
sig = 0.05 # significance level (alpha)
#ratio = 100/50 # ratio of cases:controls = 1:ratio
K = 0.333 # ratio of cases:controls = 1:ratio
n = 1000+50 # Sample size

b = NULL
for (i in xaix3) {

  res = binaryPCOld(rsq=rsq, sig=sig, K= K, n=n, b1=i)
  b = c(b, res)
}


plot(xaix3, a, ylim= c(0.0, 1.0), col='red', pch = 5, cex = 0.5)
par(new =TRUE)
plot(xaix3, b, ylim= c(0.0, 1.0), col='blue', pch = 15,cex = 0.5)



binaryPCOld(rsq=rsq, sig=sig, K= K, n=n, b1=1.315789  )





seq(from=0, to = 2.5, by = 0.01)








rsq = 0.02 # squared correlation
b1 = 0.2 # causal effect (log odds ratio per SD
b1 = log(1.2) # or log of OR per SD)
sig = 0.05 # significance level (alpha)
pow = 0.8 # power level (1-beta)
ratio = 1 # ratio of cases:controls = 1:ratio
cat("Sample size required for ", pow*100, "% power: ",
    (qnorm(1-sig/2)+qnorm(pow))^2/b1^2/rsq/(ratio/(1+ratio))/(1/(1+ratio)))
n = 40000 # Sample size
cat("Power of analysis with ", n, "participants: ",
    pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2)))














############## Description ##########################################
# Power calculation using the function
# power.calculator(alpha=0.05,beta,V_x,V_y,r_sq,nsample,sig_uv=NA)

#Input
#alpha:the specified significance level
#beta: the true causal effect
#V_x: the variance of the risk exposure (X)
#V_y: the variance of the outcome (Y)
#r_sq: the proportion of the total variation of exposure X that can
#     explained by IVs (Z)
#nsample: a 3-dimentional vector (n1,n2,n3), with n1 being the number of
#         subjects having measure on (Y,Z), n2 being the number of subjects
#         having measures on (Y,X,Z), and n3 being the number of subjects
#         having measure on (X,Z)
#sig_uv: the covariance of two error terms u and v (see definition in the
#        paper), can be null if n2=0

#Output
# the estimated power

# Example

# For an independent two-sample study (i.e., the outcome data and risk
# exposure data have no overlapping subjects), there is no need to provide
# value for sig_uv. Call the power calculation function as the following,
# power.calculator(alpha=0.05,beta,V_x,V_y,r_sq,c(n1,0,n3))

# For a one-sample study or a two-sample study with overlapping subjects, we
# must provide a value for sig_uv. For a one-sample study, call the power
# calculation function as the following,
# power.calculator(alpha=0.05,beta,V_x,V_y,r_sq,c(0,n2,0),sig_uv).
# For a two-sample study with overlapping subjects, use the following call,
# power.calculator(alpha=0.05,beta,V_x,V_y,r_sq,c(n1,n2,n3),sig_uv)

######### Source code for the R function ###################
power.calculator <- function(alpha=0.05,beta,V_x,V_y,r_sq,nsample,sig_uv=NA){
  n1<- nsample[1]; n2<- nsample[2]; n3<- nsample[3]
  if (n2!=0&is.na(sig_uv)){
    stop('the covariance of the two error terms much be given.')
  }
  if (n2==0){ sig_uv <- 0}

  V_v <- V_x*(1-r_sq)
  V_w <- V_y-V_x*r_sq*beta^2
  delta <- sig_uv/V_v
  nsigma <- (V_w/V_v+(n1+n2)*beta^2/(n2+n3) - 2*n2*(beta^2+beta*delta)/(n2+n3))   /   ((n1+n2)*r_sq/(1-r_sq))


  pnorm(-qnorm(1-alpha/2)+beta/sqrt(nsigma))+ pnorm(-qnorm(1-alpha/2)-beta/sqrt(nsigma))
}


power.calculator <- function(alpha=0.05,beta,V_x,V_y,r_sq,nsample,sig_uv=NA){
  n1<- nsample[1]; n2<- nsample[2]; n3<- nsample[3]
  if (n2!=0&is.na(sig_uv)){
    stop('the covariance of the two error terms much be given.')
  }
  if (n2==0){ sig_uv <- 0}

  V_v <- V_x*(1-r_sq)
  V_w <- V_y-V_x*r_sq*beta^2
  #delta <- sig_uv/V_v
  nsigma <- (V_w/V_v+(n1+n2)*beta^2/(n2+n3) )/((20+0)*0.2/(1-0.2))


  pnorm(-1.95 +beta/sqrt(nsigma))+ pnorm(-1.95-beta/sqrt(nsigma))
}

pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2))



pnorm( sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio))) * log(b1) - qnorm(1-sig/2)    )

pnorm( sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio))) * log(b1) - qnorm(1-sig/2)    ) +
  pnorm(-sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio))) * log(b1) - qnorm(1-sig/2))










