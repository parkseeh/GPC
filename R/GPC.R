#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
GPC <- function(x, ...) {
  UseMethod("GPC")
}


#GPC(OR = c(1.2,1.3,1.4,1.5,1.6), maf = c(0.1, 0.2, 0.3), N = 500, model ='linear')
#GWASPowerCalculator(OR = c(1.2,1.3,1.4,1.5,1.6), maf = c(0.1, 0.2, 0.3), N = 500, model ='linear')
#GPC.default(OR = c(1.2,1.3,1.4,1.5,1.6), maf = c(0.1, 0.2, 0.3), N = 500, model ='linear')


#' @export
GPC.default <- function(x, ...) {
  GWASPowerCalculator(x, ...)
}


#' Title
#'
#' @param beta
#' @param maf
#' @param N
#' @param pval
#' @param model
#'
#' @return
#' @export
#'
#' @examples
GWASPowerCalculator <- function(OR, maf, N, pval=5e-8, model='binary', Ncase,...) {

  # OR = c(1.2,1.3,1.4,1.5,1.6)
  # maf = c(0.1, 0.2, 0.3)
  # N = 500


  if (missing(OR)) {
    stop("The parameter 'OR' must be provided as numeric vector type")
  }
  if (missing(maf)) {
    stop("The parameter 'maf' must be provided as numeric vector type")
  }
  if (missing(N)) {
    stop("The parameter 'N' should be provided as numeric vector type")
  }
  if (!model %in% c('linear', 'binary')) {
    stop("The parameter 'model' must be either 'linear' or 'binary'")
  }


  power <- purrr::map(OR, findPower, maf, N, pval, model = model, Ncase = Ncase)
  power = do.call("cbind",power)
  colnames(power) = OR
  rownames(power) = maf

  #attr(power, "model") <- model

  result <- list(power = power,
                 model = model)


  #class(power) <- 'GPC'
  return (result)


}


#' Title
#'
#' @param OR
#' @param maf
#' @param N
#' @param pval
#' @param model
#' @param Ncase
#'
#' @return
#' @export
#'
#' @examples
findPower <- function(OR, maf, N, pval=5e-8, model, Ncase) {
  beta <- log(OR)
  if (model == 'linear') {
    sigma <- sqrt(1 - 2*maf*(1-maf)*beta^2)
    SE <- sigma / sqrt(2*maf*(1-maf)*N)

  } else if (model == 'binary') {
    if (missing(Ncase)) {
      stop("For binary power calculation, the number of cases must be provided on 'Ncase' parameter")
    }
    phi <- Ncase / N
    SE <- 1/sqrt(2*maf*(1-maf)*N*phi*(1-phi))
  } else {
    stop("The paramter 'model' should be either 'linear' or 'binary'")
  }

  NCP <- (beta/SE)^2
  q.thresh <- qchisq(p = pval, df = 1, ncp = 0, lower.tail = FALSE)
  power <- pchisq(q = q.thresh, df = 1, ncp = NCP, lower.tail = FALSE)


  return(power * 100)
}





