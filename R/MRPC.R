#' Mendelian Randomization Power Calculation
#'
#' @description
#' It calculate the Mendelian randomization power based on the odds ratio, the variance explained by genotype,
#' and number of case and control.
#'
#' The power to detect a causal effect can be calculated as :
#' \itemize{
#'  \item {For linear} : {\eqn{\Phi(\beta_1\rho_{GX} \sqrt{N} -z_{(1- {\alpha\over 2})})}}
#'  \item {For binary} : {\eqn{\Phi(\beta_1\rho_{GX} \sqrt{N P(Y=1)P(Y=0)} -z_{(1- {\alpha\over 2})})}}
#'}
#' where \eqn{\rho_{GX}} is the variance explained by the genotype (\eqn{\approx 2f(1-f)\beta^2}).
#' It can be calculated using effect size \eqn{\beta} and minor allele frequency.
#'
#'
#' @param OR odds ratio
#' @param rsq The variance explained by the additive effect on the genotype
#' @param N The total sample number
#' @param Ncase The number of case
#' @param Ncontrol The number of control
#' @param outcome either 'linear' or 'binary'
#' @param K the number of instrumental variable
#' @seealso \link{https://pubmed.ncbi.nlm.nih.gov/24608958/}
#' @export
#'
#' @examples
#' MRPC(OR=c(1.2,1.3,1.4,1.5,1.6), rsq=c(0.1, 0.2, 0.3), N=1000, Ncase=500, Ncontrol=500, outcome='binary', K=c(1,2,3))
#' MRPC(OR=c(1.2,1.3,1.65,1.87), rsq=c(0.111111111, 0.022), N=1000, Ncase=500, Ncontrol=500, outcome='binary', K=c(1,2))
MRPC <- function(x, ...) {
  UseMethod("MRPC")
}


#' @export
MRPC.default <- function(x, ...) {
  MRPowerCalculator(...)
}



#' @describeIn MRPC.default The \code{default} interface.
#' @importFrom purrr map
#' @export
MRPowerCalculator <- function(OR, rsq, N, pval=0.05, Ncase=NULL, Ncontrol=NULL, outcome='binary', K=1) {

  if (missing(OR)) {
    stop("The parameter 'OR' must be provided as numeric vector type")
  }
  if (missing(rsq)) {
    stop("The parameter 'rsq' must be provided as numeric vector type")
  }
  if (missing(N)) {
    stop("The parameter 'N' should be provided as numeric vector type")
  }
  if (!outcome %in% c('linear', 'binary')) {
    stop("The parameter 'outcome' must be either 'linear' or 'binary'")
  }
  if (length(K) == 1 && length(K) != length(rsq)) {
    K = rep(K[1], time=length(rsq))
    message("The number of IV(s) should be same size as rsq")
  }
  f.statistics <- ((N-K-1) / K) * (rsq/(1-rsq))

  power <- purrr::map(OR, MRFindPower, rsq, N, pval,Ncase = Ncase, Ncontrol=Ncontrol, outcome = outcome)
  #power <- paste0(format(round(power*100,2), nsamll=2), "%")
  power = do.call("cbind",power)
  colnames(power) = OR
  rownames(power) = rsq
  attr(power, 'method') <- 'MR'

  result <- list(power = power, outcome = outcome, f.statistics=f.statistics, K=K)
  class(result) <- 'GPC'

  return (result)
}



#' @describeIn MRPC
#' @importFrom stats pnorm qnorm
#' @export
MRFindPower <- function(OR, rsq, N, pval=0.05, Ncase=NULL, Ncontrol=NULL, outcome) {
  if (outcome == 'binary') {
    if (is.null(Ncase) | is.null(Ncontrol)) {
      stop("Please provide the number of Case and Control")
    }
    ratio <- Ncase/Ncontrol
    beta <- log(OR)
    lower.tail <- ifelse(sign(beta)==-1, FALSE, TRUE)
    power <- pnorm(sqrt(N*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*beta-qnorm(1-pval/2, lower.tail = lower.tail),
                   lower.tail = lower.tail)
  } else{
    beta <- OR
    lower.tail <- ifelse(sign(beta)==-1, FALSE, TRUE)
    power <- pnorm(sqrt(N*rsq)*beta-qnorm(1-pval/2, lower.tail = lower.tail), lower.tail = lower.tail)
  }
  return(power)
  #return(paste0(format(round(power*100,2),nsmall = 2), "%"))
}

