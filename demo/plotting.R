
## Install egoOptim

devtools::install_github('oonyambu/egoOptim')


library(tidyverse)
comp <- function (fun, lower, upper, ..., p = NULL, rho = 0.3, maximize = FALSE,
                  reps = 20L, file = NULL, control = list(), overwrite = FALSE)
{
  fun_name <- gsub("\\W", "", deparse1(substitute(fun)))
  time_file <- sub("\\.([^.]+$)", "_time.\\1", file)
  nsteps <- if (is.null(nsteps <- control$nsteps))
    5
  else nsteps
  if (is.null(control$budget))
    control$budget <- 50
  if (!is.null(file)) {
    if (file.exists(file) && overwrite) {
      file.remove(c(file, time_file))
    }
    cat("\n", strrep("\n#", 80), "\n", sep = "", file = time_file,
        append = TRUE)
    cat(strrep("#", 80), "\n# nsteps: ", nsteps, "\n# rho: ",
        rho, "\n# Lower Bound: (", toString(lower), ")",
        "\n# Upper Bound: (", toString(upper), ")", "\n# Budget: ",
        control$budget, "\n", file = file, append = TRUE,
        sep = "")
  }
  optimal <- control$trueglobal
  if (missing(lower)) {
    dom <- domain(fun)
    fun <- getFromNamespace(fun, "egoOptim")
    if (is.function(dom)) {
      if (is.null(p))
        stop("the dimension `p` must be provided for ",
             fun_name)
      else dom <- dom(p)
    }
    if (is.null(optimal))
      optimal <- dom$opt$f[1]
    lower <- if (!is.null(dom$lower))
      dom$lower
    else rep(0, p)
    upper <- if (!is.null(dom$upper))
      dom$upper
    else rep(1, p)
  }
  if (maximize & is.null(optimal))
    optimal <- -1
  control$trueglobal <- optimal
  res <- setNames(vector("list", 3), c("RSO", "EGO", "TREGO"))
  errors_list <- res
  RScontrol <- modifyList(control, list(basicEGO = FALSE))
  EGcontrol <- modifyList(control, list(basicEGO = TRUE))
  TRcontrol <- modifyList(control, list(basicEGO = TRUE, method = "TREGO"))
  time <- matrix(NA, nrow = reps, ncol = 3 + length(control$expansion_rate),
                 dimnames = list(NULL, unique(c(names(res), paste0("RSO",
                                                                   control$expansion_rate)))))
  parallel::mclapply(seq_len(reps), function(i) {
    X <- lhs::maximinLHS(5 * length(lower), length(lower))
    X1 <- mapply(rescale, data.frame(X), data.frame(rbind(lower,
                                                          upper)))
    cat("\n\nComputing f(x)...")
    y1 <- apply(X1, 1, function(x) (-1)^(maximize) * fun(x,
                                                         ...))
    cat("Done\nRSO ITERATION:", i, "\n")
    t1 <- proc.time()[["elapsed"]]
    res[["RSO"]][[i]] <<- optimize_fun(fun, lower, upper,
                                       ..., X = X1, y = y1, maximize = maximize, rho = rho,
                                       control = modifyList(RScontrol, list(expansion_rate = 0)))
    time[i, "RSO"] <<- proc.time()[["elapsed"]] - t1
    if (!is.null(file)) {
      cat("RSO", time[i, "RSO"], "\n", file = time_file,
          append = TRUE)
      cat("RSO", i, res[["RSO"]][[i]]$errors, "\n", file = file,
          append = TRUE)
    }
    cat("EGO ITERATION:", i, "\n")
    t2 <- proc.time()[["elapsed"]]
    res[["EGO"]][[i]] <<- optimize_fun(fun, lower, upper,
                                       ..., X = X1, y = y1, maximize = maximize, rho = rho,
                                       control = EGcontrol)
    time[i, "EGO"] <<- proc.time()[["elapsed"]] - t2
    if (!is.null(file)) {
      cat("EGO", time[i, "EGO"], "\n", file = time_file,
          append = TRUE)
      cat("EGO", i, res[["EGO"]][[i]]$errors, "\n", file = file,
          append = TRUE)
    }
    cat("TREGO ITERATION:", i, "\n")
    t3 <- proc.time()[["elapsed"]]
    res[["TREGO"]][[i]] <<- optimize_fun(fun, lower, upper,
                                         ..., X = X1, y = y1, maximize = maximize, rho = rho,
                                         control = TRcontrol)
    time[i, "TREGO"] <<- proc.time()[["elapsed"]] - t3
    if (!is.null(file)) {
      cat("TREGO", time[i, "TREGO"], "\n", file = time_file,
          append = TRUE)
      cat("TREGO", i, res[["TREGO"]][[i]]$errors, "\n",
          file = file, append = TRUE)
    }
    if (!is.null(control$expansion_rate) && any(control$expansion_rate >
                                                0)) {
      for (j in control$expansion_rate) {
        nms <- paste0("RSO", j)
        cat(nms, "ITERATION: ", i, "\n", sep = "")
        t4 <- proc.time()[["elapsed"]]
        res[[nms]][[i]] <<- optimize_fun(fun, lower,
                                         upper, ..., X = X1, y = y1, maximize = maximize,
                                         rho = rho, control = modifyList(RScontrol,
                                                                         list(expansion_rate = j)))
        time[i, nms] <<- proc.time()[["elapsed"]] -
          t4
        if (!is.null(file)) {
          cat(nms, time[i, nms], "\n", file = time_file,
              append = TRUE)
          cat(nms, i, res[[nms]][[i]]$errors, "\n",
              file = file, append = TRUE)
        }
      }
    }
  }, mc.cores = if (.Platform$OS.type == "windows")
    1
  else parallel::detectCores())
  len <- (control$budget - if (is.null(n <- control$n))
    10
    else n)/nsteps
  lapply(res, function(x) sapply(x, getElement, "errors"))
  res

}
environment(comp) <- asNamespace("egoOptim")

lineplot <- function(dat, transformer = identity, nstep_max = NULL, dir='..'){
  dat %>%
    mutate(y = match.fun(transformer)(y))%>%
    dplyr::filter(nstep <= if(is.null(nstep_max)) max(nstep) else nstep_max) %>%
    ggplot(aes(nstep, y, color = method)) +
    stat_summary(geom = 'pointrange',fun.data = ~mean_se(.x, 2)) +
    stat_summary(geom = 'line', fun = mean, linewidth=1)
}

plot_data <- function(data_path, transformer = identity){
  read.table(data_path)|>
    pivot_longer(!V1:V2, names_to = 'nstep',
                 names_transform = ~(parse_number(.x)-2)*5,
                 values_to = 'y')|>
    rename(method=V1)|>
    lineplot(transformer = transformer)
}


# branin:
domBranin <- egoOptim::domain('branin')
branin_result <- comp(egoOptim::branin, domBranin$lower,
     domBranin$upper, file = 'branin.txt',
     overwrite = TRUE, control = list(trueglobal = domBranin$opt$f))

### If using the stored data
plot_data('branin.txt', log10)

## Otherwise you can directly use `plotComparison` function
egoOptim::plotComparison(branin_result)
egoOptim::plotComparison(branin_result, n = 30)
# In the paper I recall using the plotComparison function and cutting at nstep=30


#Camel 6
domcamel6 <- egoOptim::domain('camel6')
camel6_result <- comp(egoOptim::camel6, domcamel6$lower,
     domcamel6$upper, file = 'camel6.txt',
     overwrite = TRUE, control = list(trueglobal = domcamel6$opt$f))
egoOptim::plotComparison(camel6_result, n = 30)


#hartmann 4
domhart4 <- egoOptim::domain('hart4') # 4d
hart4_result <- comp(egoOptim::hart4, domhart4$lower,
     domhart4$upper, file = 'hart4.txt',
     overwrite = TRUE, control = list(trueglobal = domhart4$opt$f))
egoOptim::plotComparison(camel6_result)



# DE Optimization

bounds <- data.frame(n = 30, m = 3, NP = c(10L, 100L),
               itermax = c(500L, 1500L),
               pMut = c(0.05, 0.95),
               pCR = c(0.05, 0.95),
               pGBest = c(0.05, 0.95))


ftest <- function(X, size, times = 10){

  replicate(times,
  do.call(Meta4Design::UniPro,
          as.list(setNames(c(size, X), names(bounds))))$opt.value)%>%
    mean()
}

lower <- unlist(bounds[1,-(1:2)])
upper <- unlist(bounds[2,-(1:2)])

### You can change the size as needed -- budget set to 50
DE_res <- comp(ftest, lower, upper, size = c(30, 3),
               file = 'ftest.txt',overwrite = TRUE)
plot_data('ftest.txt')


