#' Mendelian Randomization Power Calculation
#'
#' @param OR odds ratio
#' @param rsq R square (The variance explanied by the )
#'
#' @return
#' @export
#'
#' @examples
#' MRPC(OR=c(1.2,1.3,1.4,1.5,1.6), rsq=c(0.1, 0.2, 0.3), N=1000, Ncase=500, Ncontrol=500, model='binary')
MRPC <- function(x, ...) {
  UseMethod("MRPC")
}


#' @export
MRPC <- function(x, ...) {
  MRPowerCalculator(...)
}



#' @describeIn GPC.default The \code{default} interface.
#' @importFrom purrr map
#' @export
MRPowerCalculator <- function(OR, rsq, N, pval=0.05, Ncase=NULL, Ncontrol=NULL, model='binary', K=NULL) {

  if (missing(OR)) {
    stop("The parameter 'OR' must be provided as numeric vector type")
  }
  if (missing(rsq)) {
    stop("The parameter 'rsq' must be provided as numeric vector type")
  }
  if (missing(N)) {
    stop("The parameter 'N' should be provided as numeric vector type")
  }
  if (!model %in% c('linear', 'binary')) {
    stop("The parameter 'model' must be either 'linear' or 'binary'")
  }
  f.statistics <- ((N-K-1) / K) * (rsq/(1-rsq))

  power <- purrr::map(OR, MRFindPower, rsq, N, pval,Ncase = Ncase, Ncontrol=Ncontrol, model = model)
  power = do.call("cbind",power)
  colnames(power) = OR
  rownames(power) = rsq
  attr(power, 'method') <- 'MR'

  result <- list(power = power, model = model, f.statistics=f.statistics)
  class(result) <- 'GPC'

  return (result)
}



#' @describeIn GPC
#' @importFrom stats qchisq pchisq
#' @export
MRFindPower <- function(OR, rsq, N, pval=0.05, Ncase=NULL, Ncontrol=NULL, model) {
  if (model == 'binary') {
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
  return(paste0(format(round(power*100,2),nsmall = 2), "%"))
}
