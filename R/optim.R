
#' @export

DE <- function(fn, domain, pop_size = 50, maxiter=1000,
               weight=0.5, pCR=0.5, trace=0){
  lower <- domain$lower
  upper <- domain$upper
  p <- length(lower)
  x <- matrix(runif(pop_size*p, lower, upper), p)
  val <- apply(x, 2, fn)
  best_val <- Inf
  for(it in 1:maxiter){
    s <- replicate(3, sample(pop_size))
    x_mut <- x[,s[,1]] + weight*(x[,s[,2]] - x[,s[,3]])
    x_mut <- pmin(pmax(x_mut, lower), upper)
    id_x <- rmultinom(pop_size, 1, rep(1,p)) | matrix(runif(p*pop_size) < pCR, p)
    x_mut <- id_x * x_mut + (1-id_x)*x
    val1 <- apply(x_mut, 2, fn)
    idx <- val1 < val
    x[,idx] <- x_mut[,idx]
    val[idx] <- val1[idx]
    if(trace){
      i <-which.min(val)
      if(val[i] < best_val){
        best_val <- val[i]
        x_opt <- x[,i]
        cat(sprintf("%5i",it),", Best: f([", toString(sprintf("%8.4f", x_opt)), "]) = ", best_val,"\n", sep="")
      }
    }
  }
  if(!trace) i <-which.min(val)
  list(par = x[,i], value=val[i])
}
fn <- egoOptim::ackley
domain <- egoOptim::domain('ackley')(20)
pop_size <- 50
DE(fn, domain, 100,1000)


#' @export

ES <- function(fn, domain,  mu, lambda, maxiter=2000, type=2){
  lower <- domain$lower
  upper <- domain$upper
  p <- length(lower)
  parents <- matrix(runif(p*mu, lower, upper), p)
  t1 = 1/sqrt(2*p)
  t2 <- 1/sqrt(2*sqrt(p))
  sg <- 1
  ii <- seq(mu)
  best_val <- Inf
  for(i in 1:maxiter){
    sg <- switch(type,
            '1' = 0.15,
            '2' = sg * matrix(exp(t1*rnorm(1) + t2*rnorm(p*mu)), p))
    if(is.null(sg) & is.numeric(type)) sg <- type
    children <- matrix(rnorm(lambda*p, parents, sg), p)
    val <- apply(children, 2, fn)
    ord <- order(val)[seq(mu)]
    parents <- pmax(pmin(children[,ord], upper), lower)
    if(val[ord[1]] < best_val){
      best_val <- val[ord[1]]
      x_opt <- parents[,1]

      cat(sprintf("%5i",i),", Best: f([", toString(sprintf("%.4f", x_opt)), "]) = ", best_val,"\n", sep="")
    }
  }
  list(par = x_opt, values = best_val)
}

#' @export
ESplus <- function(fn, domain,  mu, lambda, maxiter=2000){
  lower <- domain$lower
  upper <- domain$upper
  p <- length(lower)
  parents <- matrix(runif(p*mu, lower, upper), p)
  val <- apply(parents, 2, fn)
  k <- sqrt(2*log(lambda))
  center <- parents[,which.min(val)]
  t1 = 1/sqrt(2*p)
  t2 <- 1/sqrt(2*sqrt(p))
  sg <-  1
  ii <- seq(mu)
  ord <- ii
  best_val <- Inf
  for(i in 1:maxiter){
    sg <- sg * matrix(exp(rnorm(p*mu, 0, t2)), p)
    children <- matrix(c(parents) + c(sg) * rnorm(lambda*p), p)
    val <- c(val[ord], apply(children, 2, fn))
    ord <- order(val)[seq(mu)]
    parents[,ord<=mu] <- parents[,ord[ord<=mu]]
    parents[,ord>mu] <- children[,ord[ord>mu] - mu]
    parents <- pmax(pmin(parents, upper), lower)
    if(val[ord[1]] < best_val){
      best_val <- val[ord[1]]
      x_opt <- parents[,1]
      cat(sprintf("%5i",i),", Best: f([", toString(sprintf("%.4f", x_opt)), "]) = ", best_val,"\n", sep="")
    }
  }
  list(par = x_opt, values = best_val)
}
# fn <- egoOptim::hart6
# p <- 2
# domain <- egoOptim::domain("hart6")
# mu <- 10
# lambda <- 50
# ES(fn, domain, mu,lambda, 10000)
# ESplus(fn, domain, mu,lambda, 10000)
