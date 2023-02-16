#' GWAS Power Calculator
#'
#' @description
#' It calculate the GWAS power based on the odds ratio, minor allele frequency,
#' number of subject and the p-value (default = 5e-8)
#'
#' @details
#' Statistical power is the probability that the test will reject the null hypothesis \eqn{H_0}
#' at the given significance threshold when the data follow a specific alternation hypothesis \eqn{H_1}.
#' In GWAS, \eqn{H_1} is specified by fixing the study design such as total sample size, case and controls counts
#' and parameters describing the variants such as MAF and effect size.
#'
#' By assuming that the sampling distribution of effect size estimate \eqn{\beta} is Normal,
#' we have \eqn{ \hat\beta \sim N(\beta, SE^2)}, where \eqn{\beta = 0}, the Wald's test statistic
#' \eqn{z = \hat\beta / SE} is distributed approximately as \eqn{z \sim N(0, 1)}, which is used to derived p-value.
#' The other way to get a p-value is from the chi-square distribution as \eqn{z^2 \sim \chi^2_1}.
#' The chi-square distribution is often used since we need to consider only the upper tail of the distribution
#' to compute p-value.
#'
#' However, we need to consider the variant has a non-zero effect such that \eqn{\beta \neq 0}.  Then
#' \eqn{z \sim N(\beta/SE, 1)} and \eqn{z^2 \sim \chi^2_1((\beta/SE)^2)}, which is called a chi-square distribution
#' with 1 degree of freedom and non-centrality parameter \eqn{NCP = (\beta / SE)^2}. When \eqn{\beta = 0},
#' we use the standard central chi-square distribution \eqn{\chi^2_1 = \chi^2_1(0)}.
#'
#' @section Ingredients of power:
#' \itemize{
#'  \item {Sample Sieze N} : {Increasing sample size increases power}
#'  \item {Effect Size \eqn{\beta}} : {Increasing the absolute value of effect size increases power}
#'  \item {Minor Allele Frequency \eqn{f}} : {Increasing MAF increases power}
#'  \item {Significant threshold \eqn{\alpha}} : {Increasing the threshold increases power}
#'  \item {Case-Control Proportion \eqn{\phi}} : {Moving \eqn{\phi} closer to 0.5 increases power}
#' }
#'
#'
#' @section For linear models:
#' For the linear model \eqn{y = \mu + x\beta + \epsilon}
#'
#' SE of \eqn{\hat\beta} is \eqn{ SE = {\sigma \over \sqrt {Var(x)n}}} = \eqn{ \sigma \over \sqrt {2f(1-f)n}} approximately,
#' where \eqn{f} is \eqn{maf}.  The  \eqn{Var(x)} is, under Hardy-Weinberg equilibrium,
#' approximately \eqn{2f(1-f)}, and \eqn{\sigma} is the standard deviation of the error term
#' \eqn{\epsilon : \sigma^2 = Var(y) - \beta^2Var(x)}. In a typical GWAS, the effects of variants
#' on the total phenotype variance are small, and then we can assume that the error variance
#' \eqn{\sigma^2 = Var(y)} approximately, which is 1 if the phenotype is processed by quantile normalization.
#'
#' Now we can write down the \eqn{NCP} of linear GWAS model as
#' \deqn{NCP = ({\beta \over SE})^2 = {2f(1-f)n\beta^2 \over \sigma^2}}
#'
#'
#' @section For binary models:
#' For binary case-control analyzed by logistics regression,
#' \eqn{SE = {1 \over \sqrt {Var(x)n\phi (1-\phi)}} = {1 \over \sqrt {2f(1-f)n\phi (1-\phi)}}}.
#'
#' Now we can write down the \eqn{NCP} of binary (logistics) GWAS model as
#' \deqn{NCP = ({\beta \over SE})^2 = 2f(1-f)n \phi (1-\phi) \beta^2}
#'
#'
#'
#'
#' @param OR A vector of odds ratio
#' @param maf A vector of minor allele frequency
#' @param N The total sample number
#' @param pval p-value; (the default value is 5e-8)
#' @param model The regression model used either linear or binary
#' @param Ncase The number of Case (applicable only if model is 'binary')
#' @param meta Put 'fix' if you want to calculate power for fixed effect of meta analysis
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(GPC)
#' GPC(OR = c(1.2,1.3,1.4,1.5,1.6), maf = c(0.1, 0.2, 0.3), N = 500, model ='linear')
#' GPC(OR = c(1.2,1.3,1.4,1.5,1.6), maf = c(0.1, 0.2, 0.3), N = 1500, model ='binary', Ncase = 500)
#' }
GPC <- function(x, ...) {
  UseMethod("GPC")
}



#' @export
GPC.default <- function(x, ...) {
  GWASPowerCalculator(...)
}



#' @describeIn GPC.default The \code{default} interface.
#' @importFrom purrr map
#' @export
GWASPowerCalculator <- function(OR, maf, N, pval=5e-8, model='binary', Ncase, meta=NULL) {

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

  power <- purrr::map(OR, findPower, maf, N, pval, model = model, Ncase = Ncase, meta = meta)
  power = do.call("cbind",power)
  colnames(power) = OR
  rownames(power) = maf
  attr(power, 'method') <- 'GWAS'

  result <- list(power = power, model = model)
  class(result) <- 'GPC'

  return (result)
}



#' @describeIn GPC
#' @importFrom stats qchisq pchisq
#' @export
findPower <- function(OR, maf, N, pval=5e-8, model, Ncase, meta = NULL) {

  beta <- log(OR)

  if (model == 'linear') {
    if (is.null(meta)) {
      if (length(N) > 1) {
        stop("Two GWAS detected. Please put only one GWAS")
      }
      N.eff <- N
    } else {
      N.eff <- sum(N)
    }
    sigma <- sqrt(1 - 2*maf*(1-maf)*beta^2)
    SE <- sigma / sqrt(2*maf*(1-maf)*N.eff)

  } else if (model == 'binary') {
    if (missing(Ncase)) {
      stop("For binary power calculation, the number of cases must be provided on 'Ncase' parameter")
    }
    phi <- Ncase / N
    if (is.null(meta)) {
      if (length(N) > 1 || length(Ncase) > 1) {
        stop("Two GWAS detected. Please put only one GWAS")
      }
      N.eff <- N*phi*(1-phi)
    } else {
      if (length(N) != length(Ncase)) {
        stop("The number of study and case are not equal")
      }
      N.eff <- sum(N*phi*(1-phi))
    }

    SE <- 1/sqrt(2*maf*(1-maf)*N.eff)
  } else {
    stop("The paramter 'model' should be either 'linear' or 'binary'")
  }

  NCP <- (beta/SE)^2
  q.thresh <- qchisq(p = pval, df = 1, ncp = 0, lower.tail = FALSE)
  power <-  pchisq(q = q.thresh, df = 1, ncp = NCP, lower.tail = FALSE)

  return(power)
  #return(paste0(format(round(power*100,2),nsmall = 2), "%"))
}



#' @param x An object of class 'GPC'
#' @export
print.GPC <- function(x) {
  obj <- x$power
  method <- attributes(obj)$method
  result <- lineCount(obj)
  line.length <- result$line.length
  head.line <- paste(rep("_", line.length+2), collapse = "")
  tail.line <- paste(rep("-", line.length+2), collapse = "")

  res <- result$res
  columnNames <- result$columnNames
  column.length <- result$column.length
  line.length <- result$line.length


  cat("\n")
  cat(centerprint(paste0(method, " Power Calculation"), width = line.length))
  cat("\n\n")
  cat(head.line, "\n")
  cat(centerprint(paste0("Odds Ratio"), width = line.length))
  cat('\n')
  cat(tail.line, '\n')

  for (i in 1:length(columnNames)) {
    if (i == 1) {
      #cat(paste0(columnNames[i]))
      cat((centerprint(columnNames[i], width = column.length[i])))
    } else {
      cat((centerprint(columnNames[i], width = column.length[i]+1)))
    }
  }
  cat('\n')
  cat(tail.line, "\n")

  for (i in 1:dim(res)[1]){
    cat(paste0(rownames(res)[i]), "|", sep="")
    #cat((centerprint(paste0(rownames(res)[i], "|"), width = column.length[1])))
    for(j in 1:length(columnNames)-1){

      cat(sapply(res[i,j], centerprint, width = column.length[j+1] + 1))
    }
    cat("\n")
  }
  cat(tail.line, '\n\n')

  if (method == 'MR') {
    cat("F-statistics :\n")
    for (i in 1:length(rownames(res))) {
      cat("For", rownames(res)[i], "--->", round(x$f.statistics[i],2), "on", x$K[i], 'number of IV(s)\n')
    }
  }

}



#' @param x A matrix object
#' @export
lineCount <- function(x) {
  method <- attributes(x)$method
  obj <- as.data.frame(x)
  rowN <- rownames(obj)

  reformat <- function(i) {
    sprintf("%.5f", as.numeric(i))
  }

  rowN <- unlist(lapply(rowN, reformat))
  obj <- apply(obj, 2,function(x) paste0(format(round(x*100,2),nsmall = 2), "%"))
  if (is.character(obj) && !is.matrix(obj) && length(obj) >= 1) {
    obj <- t(as.matrix(obj, byrow=F))
    rownames(obj) <- rowN
  } else if (is.matrix(obj)) {
    rownames(obj) <- rowN
  }


  methodName <- ifelse(method=='GWAS', 'MAF', "Rsq")
  columnNames <- c(methodName, colnames(obj))
  columnNum <- nchar(columnNames)
  column.nchar <- c(max(nchar(rownames(obj))),unname(apply(obj,2, function(x) max(nchar(x)))))


  column.length <- apply(rbind(columnNum, column.nchar), 2, max)
  line.length <- sum(column.length) + length(columnNames) - 1

  result <- list(res = obj,
                 columnNames = columnNames,
                 column.length = column.length,
                 line.length = line.length)

  return(result)
}


#' @param x A string
#' @param width A length of string
#' @export
centerprint <- function(x,...,width=10){
  mwidth <- max(nchar(x),width)
  sp <- (mwidth-nchar(x))/2
  front <- end <- ""
  front <- space(ceiling(sp))
  end <- space(floor(sp))
  x <- paste(front, x, end, sep="")
  return(x)
}


#' @param num an integer
#' @export
space <- function(num){
  ret <- c()
  if (num < 1) {
    return(ret)
  }
  for (i in 1:num) {
    ret <- paste(" ", ret, sep="")
  }
  return(ret)
}
