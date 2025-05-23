#' @useDynLib UniPro, .registration=TRUE
#' @importFrom utils modifyList head tail
NULL

#' Differential Evolution
#'
#' DE: A function to construct Designs using various criterions
#' algorithm
#' @param n number of rows(observations)
#' @param m number of columns(variables/factors)
#' @param s number of levels. Default = n
#' @param p the exponent for method maximinLHD. default is 15
#' @param NP initial population size. Default = 100
#' @param itermax maximum number of iterations. Default = 1000
#' @param pMut the probability of mutation
#' @param pCR Probability of CrossOver
#' @param pGBest Probability of using global best
#' @param replicates number of replications. Determines the output
#' @param seed seed to be set for reproducibility
#' @param trace Non-negative integer. If positive, tracing information on the progress of the optimization is produced.
#'  Higher values may produce more tracing information.
#' @param ncores number f cores to be used for parallelization/multi-threading.
#' @param method The criterion to be optimized. Only UniPro, maximinLHD, maxPro are implemented
#' @return A list with the phi value, the final design and total computational time taken:
#' \item{bestX}{The final UniPro Design. Only returned when replicates = 1}
#' \item{optValues}{The phi value of the final Design. A vector of length replicates}
#' \item{timeTaken}{Total computational time taken}
#' @examples
#' DE(10, 3) # uses UniPro by default
#' DE(30, 3, replicates=10, method = 'maximinLHD')
#' @usage DE(n, m, s = n, p = 15L, NP = 100L, itermax = 1500L, pMut = NULL, pCR = NULL,
#'           pGBest = NULL, replicates = 1L, seed = sample(1e7,1), ncores = NULL,
#'           trace = FALSE, method = c("UniPro", "maximinLHD"))
#'
#' @name DE
#'
NULL

#' Uniform Projection Designs
#'
#' UniPro: A function to construct a Uniform Projection Design using Differential Evolution
#' algorithm
#'
#' @usage UniPro(n, m, s = n, NP = 100L, itermax = 1500L, pMut = NULL, pCR = NULL,
#'        pGBest = NULL, replicates = 1L, seed = sample(1e7,1), ncores = NULL,
#'        trace = FALSE)
#'
#' @rdname DE
#' @examples
#' UniPro(10, 3)
#' UniPro(30, 3, replicates=10)
#'
#' @export
#'
#'


UniPro <-function(n, m, s = n, NP = 100L, itermax = 1500L, pMut = NULL,
                  pCR = NULL, pGBest = NULL, replicates = 1L,
                  seed = sample(1e7,1), ncores = NULL, trace = FALSE){
  DE(n = n, m = m, s = s, NP = NP, itermax = itermax, pMut = pMut, pCR = pCR,
     pGBest = pGBest, replicates = replicates, seed = seed, ncores = ncores,
     trace = trace, method = 'UniPro')
}

#' @export
DE <-function(n, m, s = n, p = 15L, NP = 100L, itermax = 1500L, pMut = NULL,
              pCR = NULL, pGBest = NULL, replicates = 1L,
              seed = sample(1e7,1), ncores = NULL, trace = FALSE,
              method = c("UniPro", "maximinLHD")){

  if(is.null(pGBest)) pGBest <- 0.95
  methods <- c("UniPro", "maximinLHD", "MaxPro")
  .method <- as.integer(charmatch(method[1], methods, 0))
  if(.method < 1 || .method > 3) stop("invalid method", method[1])
  if(is.null(ncores)) ncores <- as.integer(max(1, detectCores() - 2))
  if(ncores < 2) stop("ncores must be at least 2")
  args <- list(n = as.integer(n), m = as.integer(m), s = as.integer(s),
               NP = as.integer(NP), itermax = as.integer(itermax),
               pMut = as.double(pMut), pCR = as.double(pCR),
               pGBest = as.double(pGBest), replicates = as.integer(replicates),
               ncores = as.integer(ncores), method = .method,
               p = as.integer(p), trace = as.integer(trace))
  u <- 0
  if(is.null(pMut) || is.null(pCR)){
    set.seed(seed)
    args1 <- modifyList(args, list(trace = 0L, replicates = 1L, itermax = 100L))
    fn <- fn <- function(pMut, pCR){
      newargs <- modifyList(args1, list(pMut = pMut, pCR = pCR, seed=sample(1e7, 1)))
      do.call(.Call, c("DE", newargs))
    }
    pmut <- seq(0.1, 0.9, by = 0.2)
    v <- mapply(fn, pmut, 1 - pmut)
    idx <- which.min(unlist(v["measure", ]))
    args$pMut <- pmut[idx]
    args$pCR <- 1 - pmut[idx]
    u <- sum(unlist(v["timeTaken",]))
  }
  set.seed(seed)
  res <- do.call(.Call, c("DE", args, PACKAGE = 'UniPro'))
  res$timeTaken <- res$timeTaken + u
  structure(res, method = method[1])
}



# #' Maximum Projection Designs
# #'
# #' MaxPro: A function to construct a Maximum Projection Designs using Differential Evolution
# #' algorithm
# #'
# #' @usage MaxPro(n, m, NP = 100L, itermax = 1000L, pMut =0.2, pCR = 0.3,
# #'        pGBest = 0.9, replicates = 1L, seed = sample(1e7,1), ncores = NULL)
# #'
# #' @rdname DE
# #' @examples
# #' MaxPro(10, 3)
# #' MaxPro(30, 3, replicates=10)
# #'
# #'
# #' @name maxpro
# NULL
#
#
#
# MaxPro <-function(n, m, NP = 100, itermax = 1500, pMut = 0.2,
#                   pCR = 0.3, pGBest = 0.9, replicates = 1,
#                   seed = sample(1e7,1), ncores = NULL){
#   DE(n = n, m = m, NP = NP, itermax = itermax, pMut = pMut,
#      pCR = pCR, pGBest = pGBest, replicates = replicates,
#      seed = seed, ncores = ncores, method = "MaxPro")
# }
#
# #' Maximin-Distance Latin Hypercube Designs
# #'
# #' maximinLHD: A function to construct  Maximin-Distance Latin Hypercube Designs using Differential Evolution
# #' algorithm
# #'
# #' @usage maximinLHD(n, m, p = 15L, NP = 100L, itermax = 1000L, pMut = 0.2, pCR = 0.3,
# #'        pGBest = 0.9, replicates = 1L, seed = sample(1e7,1), ncores = NULL)
# #'
# #' @examples
# #' maximinLHD(10, 3)
# #' maximinLHD(30, 3, replicates=10)
# #'
# #' @name maximin
# NULL
#
#
#
# maximinLHD <-function(n, m, p = 15L, NP = 100, itermax = 1500, pMut = 0.2,
#                   pCR = 0.3, pGBest = 0.9, replicates = 1,
#                   seed = sample(1e7,1), ncores = NULL){
#   DE(n = n, m = m, p = p, NP = NP, itermax = itermax, pMut = pMut,
#      pCR = pCR, pGBest = pGBest, replicates = replicates,
#      seed = seed, ncores = ncores, method = "maximinLHD")
# }

printFun <- function(x){
  len <- length(x)
  y <- c(sprintf("%5.3f",head(x, 3)),
         if(len > 6) "...",
         sprintf("%5.3f",tail(tail(x, -3), 3)))
  cat(" measure", if(len > 1)c(" [1:",len,"]"),": [", toString(y), "]\n", sep="")
}


#' @exportS3Method base::print
print.DE <- function(x, ...){
  dm <- dim(x$Design)
  cat(" method:", attr(x, "method"))
  cat("    Size: ", dm[1], "x", dm[2], "\n", sep="")
  cat(sprintf(" Total Time Taken: %8.4f Secs\n", x$timeTaken))
  printFun(x$measure)
  invisible(x)
}


#' unipromeasure
#'
#' A function that computes the UniPro criterion value for a balanced design.
#' Equation 5.1 and 5.2 of the Uniform Projection Designs paper by Sun, Wang and Xu.
#' @usage unipromeasure(X, s = nrow(X))
#' @param X Design
#' @param s number of levels. Default = n
# #' @param method The criterion to be optimized. Only UniPro, maximinLHD are implemented
#' @return measure value.
#' @examples
#' unipromeasure(replicate(3, sample(10)))
#'
#' @export
#'
unipromeasure <- function(X, s = nrow(X)){
  measure(x = X, s = s, method = 'UniPro')
}


# NULL
#
#
measure <- function(x, s = nrow(x), p = 15L, method = c("UniPro", "maximinLHD")){
  m = charmatch(method[1], c("UniPro", "maximinLHD", "MaxPro"), 0L)
  if(!m) stop("not implemented for", method[1])
  .Call("phi", `mode<-`(x, "integer"), as.integer(s),
        as.integer(p), as.integer(m), PACKAGE = 'UniPro')
}

#
#
#
# #' maxpromeasure
# #'
# #'
# #' A function that computes the MaxPro criterion value for a given design.
# #' Equation 5 of the  Maximum projection designs for computer experiments
# #' paper BY V. ROSHAN JOSEPH,EVRENGUL
# #'
# #' @usage maxpromeasure(X)
# #' @examples
# #' maxpromeasure(replicate(3, sample(10)))
# #'
# #' @name maxpromeasure
# NULL
#
# maxpromeasure <- function(x){
#   measure(x = x, method = 'MaxPro')
# }
#
# #' maxminmeasure
# #'
# #'
# #' A function that computes the maximin criterion value for a given design.
# #' (Morris and Mitchell 1995)
# #'
# #' @usage maximinmeasure(X)
#
# #' @examples
# #' maximinmeasure(replicate(3, sample(10)))
# #'
# #' @name maximinmeasure
# NULL
#
# maximinmeasure <- function(x, p = 15){
#   measure(x = x, p = p, method = 'maximinLHD')
# }
#
#
detectCores <- function(){
  .Call("detectCores", PACKAGE = 'UniPro')
}
