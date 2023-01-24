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
    GWASPowerCalculator(...)
}


#' Title
#'
#' @param maf
#' @param N
#' @param pval
#' @param model
#' @param OR
#' @param Ncase
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
GWASPowerCalculator <- function(OR, maf, N, pval=5e-8, model='binary', Ncase) {

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


    result <- list(power = power,
                   model = model)


    class(result) <- 'GPC'
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
    fmt <- sprintf("%s%df","%3.",3)

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
    power <-  pchisq(q = q.thresh, df = 1, ncp = NCP, lower.tail = FALSE)


    return(paste0(round(power*100,2), "%"))
}

lineCount <- function(x) {
    obj <- as.data.frame(x)
    columnNames <- c("MAF|", colnames(obj))
    columnNum <- nchar(columnNames)
    column.nchar <- c(max(nchar(rownames(obj))),unname(sapply(obj,function(x) max(nchar(x)))))

    column.length <- apply(rbind(columnNum, column.nchar), 2, max)
    line.length <- sum(column.length) + length(columnNames) - 1


    result <- list(res = obj,
                   columnNames = columnNames,
                   column.length = column.length,
                   line.length = line.length)

    return(result)

}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
print.GPC <- function(x) {
    obj <- x$power

    result <- lineCount(obj)
    line.length <- result$line.length
    head.line <- paste(rep("_", line.length+1), collapse = "")
    tail.line <- paste(rep("-", line.length+1), collapse = "")

    res <- result$res
    columnNames <- result$columnNames
    column.length <- result$column.length
    line.length <- result$line.length

    cat("\n")
    cat(centerprint(paste0("GWAS Power Calculation"), width = line.length))
    cat("\n\n")
    cat(head.line, "\n")
    cat(centerprint(paste0("Odds Ratio"), width = line.length))
    cat('\n')
    cat(tail.line, '\n')


    for (i in 1:length(columnNames)) {
        if (i == 1) {
            cat(paste0(columnNames[i]))
        } else {
            cat((centerprint(columnNames[i], width = column.length[i]+1.5)))
        }

    }
    cat('\n')
    cat(tail.line, "\n")


    for (i in 1:dim(res)[1]){
        cat(paste0(rownames(res)[i]), "|", sep="")
        for(j in 1:length(columnNames)){

            cat(sapply(res[i,j], centerprint, width = column.length[j] + 2))
        }
        cat("\n")
    }
    cat(tail.line, '\n\n')



}


#' Prints the string in the center within the width value
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



#findPower(OR = c(1.2), maf = c(0.1), N = 500, model ='linear')
#GPC(OR = c(1.2,1.3,1.4,1.5,1.6), maf = c(0.1, 0.2, 0.3), N = 500, model ='linear')
#GWASPowerCalculator(OR = c(1.2,1.3,1.4,1.5,1.6), maf = c(0.1, 0.2, 0.3), N = 500, model ='linear')
#GPC.default(OR = c(1.2,1.3,1.4,1.5,1.6), maf = c(0.1, 0.2, 0.3), N = 500, model ='linear')



