setMethod("initialize", "RankData",
          function(.Object, ...){
            arg = list(...)
            fields = names(arg)
            # init ranking
            if("ranking" %in% fields){
              .Object@ranking = arg[["ranking"]]
            } else if ("ordering" %in% fields){
              .Object@ranking = OrderingToRanking(arg[["ordering"]])
            } else {
              stop("Either ordering or ranking matrix should be given")
            }
            # init nobj
            if("nobj" %in% fields){
              .Object@nobj = arg[["nobj"]]
            } else {
              .Object@nobj = max(.Object@ranking)
            }
            # init count
            if("count" %in% fields){
              .Object@count = arg[["count"]]
            } else {
              .Object@count = rep(1,nrow(.Object@ranking))
            }
            # init nobs
            if("nobs" %in% fields){
              .Object@nobs = arg[["nobs"]]
              stopifnot(.Object@nobs==sum(.Object@count))
            } else {
              .Object@nobs = sum(.Object@count)
            }
            # init ndistinct
            .Object@ndistinct = nrow(.Object@ranking)
            # handle topq
            if ("topq" %in% fields){  # the field is given
              if (min(arg[["topq"]]) < .Object@nobj-1){  # true topq case
                if (any(arg[["topq"]]>=.Object@nobj)){
                  warning("topq value should range between 1 and nobs-1")
                  arg[["topq"]][arg[["topq"]]>=.Object@nobj] = .Object@nobj-1
                }
                max_rank = apply(.Object@ranking,1,max)
                max_rank = as.numeric(max_rank) - 1
                if (!setequal(max_rank,arg[["topq"]])){
                  warning("The supplied top-q vector is not valid")
                  .Object@topq = unique(max_rank)
                } else {
                  .Object@topq = arg[["topq"]]
                }
                .Object@q_ind = c(1,cumsum(as.numeric(table(max_rank)[as.character(.Object@topq)]))+1)
                .Object@subobs = numeric(length(arg[["topq"]]))
                for (i in 1:length(arg[["topq"]])){
                  .Object@subobs[i] = sum(.Object@count[ .Object@q_ind[i]: (.Object@q_ind[i+1]-1) ])
                }
              } else {
                .Object@q_ind = -1
                .Object@subobs = -1
              }
            } else {
              .Object@topq = .Object@nobj - 1
              .Object@q_ind = -1
              .Object@subobs = -1
            }
            .Object
          }
)




setGeneric("SingleClusterModel",
           def=function(dat,init,ctrl,modal_ranking){
             standardGeneric("SingleClusterModel")
             })

# single cluster model method for Weighted Kendall Distance
setMethod(
  "SingleClusterModel",
  signature = c("RankData","RankInit","RankControlWeightedKendall"),
  definition = function(dat,init,ctrl,modal_ranking) {

    param_len = max(dat@topq)
    paramCoeff = CWeightGivenPi(dat@ranking,modal_ranking)
    paramCoeff = matrix(paramCoeff,ncol = dat@ndistinct,byrow = TRUE)

    #The iteration for EM algorithm.
    itr <- 0
    cores <- parallel::detectCores()

    while(TRUE){
      if (itr == 0){
        #No previously calculated parameters. Use the supplied initial values.
        params <- init@param.init[[1]]
        alphaTau <- unlist(init@alpha.init)
      }
      else{
        #Exist previously calculated parameters. Use the previous results.
        params <- param.est
        alphaTau <- unlist(alpha.est)
      }
      product <- apply(params*paramCoeff, 2, sum)
      medium <- product + alphaTau
      listPDF <- NULL
      leng <- ncol(paramCoeff)
      #prob is used to store the probability of each pi(i).
      prob <- numeric(length = leng)

      #for (i in 1:leng){
      #  listPDF <- c (listPDF,
      #                local({
      #                  i <- i;
      #                  function (x) adjustedPDFMono(x = x, mediumMono = medium[i],
      #                                               alpha = alphaTau, param = params)}))
      #}
      f <- function(mediumMono){
        function(x){
          adjustedPDFMono(x = x, mediumMono = mediumMono, alpha = alphaTau, param = params)
        }
      }
      listPDF <- parallel::mclapply(medium, f, mc.cores = cores)
      #Now we have listPDF as a list of PDFs.

      #For each iteration, we need to calculate a new set of M and N, since both are
      #dependent on timesteps.
      funcM <- NULL
      for (i in 1:leng){
        funcM <- c(funcM, local({ i <- i;
          function(x) {listPDF[[i]](x) * (-x)}
        }))
      }
      #f2 <- function(listPDFMono){
      #  function(x){
      #    listPDFMono(x) * (-x)
      #  }
      #}
      #funcM <- parallel::mclapply(listPDF, f2,mc.cores = cores)
      #print("c")
      funcN <- NULL
      for (i in 1:leng){
        funcN <- c(funcN, local({ i <- i;
          function(x) {listPDF[[i]](x) * (log(x))}
        }))
      }
      #f3 <- function(listPDFMono){
      #  function(x){
      #    listPDFMono(x) * (log(x))
      #  }
      #}
      #browser()
      #funcN <- lapply(listPDF, f3)

      #M <- numeric(leng)
      #for (i in 1:leng){
      #  M[i] <- pracma::integral(funcM[[i]], xmin = 0, xmax = Inf, reltol = 3e-7)
      #}
      f1 <- function(f, xmin, xmax, reltol){
        pracma::integral(f, xmin = xmin, xmax = xmax, reltol = reltol)
      }

      M <- unlist(parallel::mclapply(funcM, f1, xmin = 0, xmax = Inf, reltol = 3e-7,
                              mc.cores=cores), use.names = FALSE)
      #NVector <- numeric(leng)
      #for (i in 1:leng){
      #  NVector[i] <- pracma::integral(funcN[[i]], xmin = 0, xmax = Inf, reltol = 3e-7)
      #}
      NVector <- unlist(parallel::mclapply(funcN, f1, xmin = 0, xmax = Inf, reltol = 3e-7,
                                    mc.cores = cores), use.names = FALSE)

      N <- NVector %*% dat@count
      param.coeff <- paramCoeff %*% (dat@count*M)
      param.coeff <- as.numeric(param.coeff)[1:param_len]
      #param.coeff records the number of cumulative adjacent swappings needed with
      #constant M incorporated.

      if (length(dat@topq) == 1 && dat@topq == dat@nobj-1) {
        #counter <- 0
        #reference <- numeric()

        obj <- function(params) {
          alpha <- params[length(params)]

          a <- params[-length(params)] %*% param.coeff + alpha*(M%*%dat@count) + dat@nobs * (alpha *
            log(alpha) - log(gamma(alpha))) + (alpha-1)*N -
            ElogC(params[-length(params)], leng) %*% dat@count
          #counter <<- counter + 1

          #if (counter == 1)
           # reference <<- a

          print(a)
          print(params)

          as.numeric(-1 * a)
        }

        ElogC <- function(param, leng){
          output <- numeric(leng)
          integrand <- NULL
          for (i in 1:leng){
            integrand <- c(integrand, local({
              i <- i;
              function(x){
                LogC(x*param)*listPDF[[i]](x)
              }
            }))
          }
          unlist(parallel::mclapply(integrand, f1,xmin = 0, xmax = Inf, reltol = 3e-7,
                             mc.cores = cores), use.names=FALSE)

        }
        #ElogC <- compiler::cmpfun(ElogC1)


        #obj <- compiler::cmpfun(obj1)
        #obj() calculates the opposite of complete log likelihood of l(pi1, pi2, ..., pin)

        #objAlf <- function (alf){
        #  a <- params %*% param.coeff + alf*(M%*%dat@count) + dat@nobs * (alf *                                                                                             log(alpha) - log(gamma(alpha))) + (alpha-1)*N -
        #    ElogCAlf %*% dat@count

        #  as.numeric(-1 * a)
        #}
        tt = t.gen(param_len)
        gradient <- function(params) {
          phis <- params[-length(params)]
          grad <- GHCExt(phis, tt, listPDF, f1, cores)
          #grad is a matrix with nobj rows and ndistinct columns.
          grad <- grad %*% dat@count
          gradt <- numeric(length(params))
          gradt[-length(params)] <- grad - param.coeff
          gradt[length(params)] <- (-1)*(M%*%dat@count) - dat@nobs*(log(params[length(params)])
                                    - digamma( params[length(params)]) + 1) - N
          gradt
        }
        #gradient <- compiler::cmpfun(gradient1)

        #grad <- function(alf){
        #  gradt <- (-1)*(M%*%dat@count) - dat@nobs*(log(alf)
        #            - digamma(alf) + 1) - N
        #}
        #reference <- (-1)*obj(c(params, alphaTau))
        #print("arrived2")
        if (itr < 2)
          fctr <- 6e9
        else
          fctr <- 5e8

        opt_res <- optim(
          par = c(params, alphaTau), fn = obj, gr = gradient, method = "L-BFGS-B",
          lower = c(rep(0, length(params)), 1e-8), upper = rep(Inf, length(params) + 1),
          control = list(factr = fctr)
        )


        #browser()
        #, control = ctrl@optimx_control
        #,lower = c(rep(0,param_len), 1e-8), upper = rep(Inf,param_len+1),
        print("values of paramsTau:")
        print(params)
        print("values of alphaTau:")
        print(alphaTau)
        #print(gradient(c(params, alphaTau)))

        param.est <- opt_res[[1]][1:param_len]
        alpha.est <- opt_res[[1]][param_len + 1]
        prob <- FindProb(dat, ctrl, modal_ranking, c(param.est, alpha.est))
        log_likelihood <- sum(log(prob) %*% dat@count)
        print("log likelihood:")
        print(log_likelihood)
      }
      else {
        obj = function(param) {
          norm_vec = numeric(length(dat@topq))
          for (i in 1:length(dat@topq)) {
            j = dat@topq[i]
            norm_vec[i] = LogC(c(param[1:j],rep(0,dat@nobj-1 - j)))
          }
          a = -1 * param %*% param.coeff - dat@subobs %*% norm_vec
          as.numeric(-1 * a)
        }
        opt_res = optimx::optimx(
          par = init@param.init[[init@clu]][1:param_len],fn = obj,lower = rep(0,param_len),upper =
            rep(Inf,param_len),method = "L-BFGS-B",control = ctrl@optimx_control
        )
        param.est = unlist(opt_res[1:param_len])
        log_likelihood = -1 * obj(param.est)
      }

      if (itr == 0){
        crit <- log_likelihood
      }
      if (any(abs(param.est) <= 1e-6)){
        print("Detected phi parameters close to 0.")
        break
      }

      if (itr > 0){
        #if (log_likelihood - reference < abs(reference)*2e-9){
        if (log_likelihood - old_likelihood < abs(old_likelihood)*3e-6){
          #The EM algorithm failed to ascend
          print(paste("failed to ascend at iteration ", itr))
          #Use the old parameter values and likelihoods and jump out of the loop
          #param.est <- params
          #alpha.est <- alphaTau
          #log_likelihood <- old_goodness
          break
        }

        #if (itr %% 25==0)
          #browser()

        #otherwise we accept the new model and carry on to the next iteration.
        #old_goodness <- log_likelihood
      }
      old_likelihood <- log_likelihood

      itr <- itr+1
      print(itr)
      print(c(opt_res[[1]], log_likelihood))

      if (itr > ctrl@EM_limit){
        message("Algorithm did not converge in ",ctrl@EM_limit,
                    " iterations in function SingleClusterModel for modal ranking",
                    modal_ranking)
        break
      }
    }
    #param.est = c(param.est,rep(0,dat@nobj - 1 - param_len))
    list(
      param.est = param.est,
      w.est = paramTow(param.est),
      alpha.est = alpha.est, log_likelihood = log_likelihood,
      criteria = crit
    )
    }
)

# single cluster model method for Kendall distance
setMethod("SingleClusterModel",
          signature = c("RankData","RankInit","RankControlKendall"),
          definition = function(dat,init,ctrl,modal_ranking){
            param.coeff <- FindV(dat@ranking,modal_ranking)
            param.coeff <- rowSums(param.coeff)%*%dat@count
            param.coeff <- as.numeric(param.coeff)

            obj <- function(param){
              param*param.coeff + dat@nobs*LogC_Component(rep(param,dat@nobj-1))
            }

            opt_res = stats::optimize(f=obj,interval =c(0,100))
            list(param.est=opt_res$minimum,log_likelihood=-1*opt_res$objective)
          }
)

# single cluster model method for Phi Component Model
setMethod("SingleClusterModel",
          signature = c("RankData","RankInit","RankControlPhiComponent"),
          definition = function(dat,init,ctrl,modal_ranking){
            param.coeff <- FindV(dat@ranking,modal_ranking)
            param.coeff <- t(param.coeff)%*%dat@count
            param.coeff <- as.numeric(param.coeff)
            param_len <- dat@nobj-1
            obj <- list()
            opt_res <- list()
            for (i in 1:param_len){
              obj[[i]] <- function(param){
                rhs <- exp(-param)/(1-exp(-param))-(dat@nobj+1-i)*exp(-(dat@nobj+1-i)*param)/(1-exp(-(dat@nobj+1-i)*param))
                ret <- param.coeff[i] - dat@nobs* rhs
                ret
              }
              opt_res[[i]] <- stats::uniroot(f=obj[[i]],interval=c(0,100),f.lower=-Inf)
            }
            param.est <- vapply(opt_res,function(x)x$root,numeric(1))
            log_likelihood <- -param.coeff%*%param.est - dat@nobs*LogC_Component(param.est)
            list(param.est=param.est,log_likelihood=log_likelihood)
          }
)

setMethod("SingleClusterModel",
          signature = c("RankData","RankInit","RankControlWtau"),
          definition = function(dat, init, ctrl, modal_ranking){
            # the same order as param.expand
            param_len <- dat@nobj
            param.coeff <- Wtau(dat@ranking, modal_ranking)
            param.coeff <- t(param.coeff)%*%dat@count
            param.coeff <- as.numeric(param.coeff)
            allperm <- AllPerms(dat@nobj)
            all_coeff <- Wtau(allperm, modal_ranking)
            obj <- function(param){
              param.expand <- outer(param,param)[upper.tri(diag(length(param)))]
              LogC <- log(sum(exp(-1*all_coeff %*% param.expand)))
              loglike <- param.coeff %*% param.expand + dat@nobs*LogC
              as.numeric(loglike)
            }
            opt_res = optimx::optimx(
              par = init@param.init[[init@clu]][1:param_len],fn = obj,
              method = "Nelder-Mead",control = ctrl@optimx_control
            )
            param.est = unlist(opt_res[1:param_len])
            log_likelihood = -1 * opt_res[[param_len + 1]]
            list(
              param.est = param.est,log_likelihood = log_likelihood
            )
          }
)

# single cluster model method for Kendall distance
setMethod("SingleClusterModel",
          signature = c("RankData","RankInit","RankControlSpearman"),
          definition = function(dat,init,ctrl,modal_ranking){
            param.coeff <- colSums((t(dat@ranking) - modal_ranking)^2)
            param.coeff <- as.numeric(param.coeff)%*%dat@count
            param.coeff <- 0.5*as.numeric(param.coeff)
            allperm <- AllPerms(dat@nobj)
            all_coeff <- 0.5*as.numeric(colSums((t(allperm) - modal_ranking)^2))

            obj <- function(param){
              LogC <- log(sum(exp(-1 * all_coeff * param)))
              param*param.coeff + dat@nobs*LogC
            }

            opt_res <- stats::optimize(f=obj,interval =c(0,100))
            list(param.est=opt_res$minimum,log_likelihood=-1*opt_res$objective)
          }
)

setMethod("SingleClusterModel",
          signature = c("RankData","RankInit","RankControlFootrule"),
          definition = function(dat,init,ctrl,modal_ranking){
            param.coeff <- apply(dat@ranking, 1, function(x){ sum(abs(x-modal_ranking))} )
            param.coeff <- as.numeric(param.coeff)%*%dat@count
            param.coeff <- as.numeric(param.coeff)
            allperm <- AllPerms(dat@nobj)
            all_coeff <- as.numeric(apply(allperm, 1, function(x){ sum(abs(x-modal_ranking))} ))

            obj <- function(param){
              LogC <- log(sum(exp(-1 * all_coeff * param)))
              param*param.coeff + dat@nobs*LogC
            }

            opt_res <- stats::optimize(f=obj,interval =c(0,100))
            list(param.est=opt_res$minimum,log_likelihood=-1*opt_res$objective)
          }
)

setMethod("SingleClusterModel",
          signature = c("RankData","RankInit","RankControlHamming"),
          definition = function(dat,init,ctrl,modal_ranking){
            param.coeff <- apply(dat@ranking, 1, function(x){sum(x != modal_ranking)} )
            param.coeff <- as.numeric(param.coeff)%*%dat@count
            param.coeff <- as.numeric(param.coeff)
            allperm <- AllPerms(dat@nobj)
            all_coeff <- as.numeric(apply(allperm, 1, function(x){sum(x != modal_ranking)} ))

            obj <- function(param){
              LogC <- log(sum(exp(-1 * all_coeff * param)))
              param*param.coeff + dat@nobs*LogC
            }

            opt_res <- stats::optimize(f=obj,interval =c(0,100))
            list(param.est=opt_res$minimum,log_likelihood=-1*opt_res$objective)
          }
)

setMethod("SingleClusterModel",
          signature = c("RankData","RankInit","RankControlCayley"),
          definition = function(dat,init,ctrl,modal_ranking){
            param.coeff <- FindCayley(dat@ranking, modal_ranking)%*%dat@count
            param.coeff <- as.numeric(param.coeff)
            allperm <- AllPerms(dat@nobj)
            all_coeff <- as.numeric(FindCayley(allperm, modal_ranking))

            obj <- function(param){
              LogC <- log(sum(exp(-1 * all_coeff * param)))
              param*param.coeff + dat@nobs*LogC
            }

            opt_res <- stats::optimize(f=obj,interval =c(0,100))
            list(param.est=opt_res$minimum,log_likelihood=-1*opt_res$objective)
          }
)

setGeneric("FindProb",
           def=function(dat,ctrl,modal_ranking,param){standardGeneric("FindProb")}
)


setMethod("FindProb",
          signature=c("RankData","RankControlWeightedKendall"),
          definition = function(dat, ctrl, modal_ranking, param){
            phis <- param[-length(param)]
            alpha <- param[length(param)]
            distance = phis %*% matrix(CWeightGivenPi(dat@ranking,modal_ranking),
                                        ncol = dat@ndistinct, byrow = TRUE)
            medium <- distance + alpha
            leng <- dat@ndistinct
            prob <- NULL
            cores <- parallel::detectCores()

            if (length(dat@topq) == 1 && dat@topq == dat@nobj-1) {

              f2 <- function(mediumMono){
                getMarginalPiMono(mediumMono = mediumMono, alpha = alpha, param = phis,
                                    reltol = 5e-7)
              }
              prob <- unlist(parallel::mclapply(medium, f2, mc.cores = cores),
                             use.names = FALSE)
            }
            else {
              #The following is left unchanged. Will be changed when incomplete rankings are
              #being considered.
              cond_prob = dat@subobs/dat@nobs
              prob = exp(-1*distance)
              for (i in 1:length(dat@topq)) {
                j = dat@topq[i]
                norm_c = exp(LogC(c(phis[1:j],rep(0,length(phis) - j))) - lgamma(dat@nobj-j+1))
                prob[dat@q_ind[i]:(dat@q_ind[i+1]-1)] = prob[dat@q_ind[i]:
                                                              (dat@q_ind[i+1]-1)]/norm_c*cond_prob[i]
              }

            }
            prob
          }
)

setMethod("FindProb",
          signature=c("RankData","RankControlKendall"),
          definition = function(dat,ctrl,modal_ranking,param){
            param = param[1]
            distance = FindV(dat@ranking,modal_ranking) %*% rep(param,dat@nobj-1)
            C = exp(LogC_Component(rep(param,dat@nobj-1)))
            prob = exp(-1*distance)/C
            prob
          }
)


setMethod("FindProb",
          signature=c("RankData","RankControlPhiComponent"),
          definition = function(dat,ctrl,modal_ranking,param){
            distance = FindV(dat@ranking,modal_ranking) %*% param
            C = exp(LogC_Component(param))
            prob = exp(-1*distance)/C
            prob
          }
)

setMethod("FindProb",
          signature=c("RankData","RankControlWtau"),
          definition<- function(dat,ctrl,modal_ranking,param){
            allperm <- AllPerms(dat@nobj)
            all_coeff <- Wtau(allperm, modal_ranking)
            param.expand <- outer(param,param)[upper.tri(diag(length(param)))]
            C <- sum(exp(-1*all_coeff %*% param.expand))
            distance<- Wtau(dat@ranking, modal_ranking) %*% param.expand
            prob<- exp(-1*distance)/C
            prob
          }
)

setMethod("FindProb",
          signature=c("RankData","RankControlSpearman"),
          definition<- function(dat,ctrl,modal_ranking,param){
            allperm <- AllPerms(dat@nobj)
            all_coeff <- 0.5*as.numeric(colSums((t(allperm) - modal_ranking)^2))
            C <- sum(exp(-1*all_coeff * param))
            param.coeff <- colSums((t(dat@ranking) - modal_ranking)^2)
            param.coeff <- 0.5*as.numeric(param.coeff)
            distance <- param.coeff * param
            prob <- exp(-1*distance)/C
            prob
          }
)


setMethod("FindProb",
          signature=c("RankData","RankControlFootrule"),
          definition<- function(dat,ctrl,modal_ranking,param){
            allperm <- AllPerms(dat@nobj)
            all_coeff <- as.numeric(apply(allperm, 1, function(x){ sum(abs(x-modal_ranking))} ))
            C <- sum(exp(-1*all_coeff * param))
            param.coeff <- apply(dat@ranking, 1, function(x){ sum(abs(x-modal_ranking))} )
            param.coeff <- as.numeric(param.coeff)
            distance <- param.coeff * param
            prob <- exp(-1*distance)/C
            prob
          }
)


setMethod("FindProb",
          signature=c("RankData","RankControlHamming"),
          definition<- function(dat,ctrl,modal_ranking,param){
            allperm <- AllPerms(dat@nobj)
            all_coeff <- as.numeric(apply(allperm, 1, function(x){sum(x != modal_ranking)} ))
            C <- sum(exp(-1*all_coeff * param))
            param.coeff <- apply(dat@ranking, 1, function(x){sum(x != modal_ranking)} )
            param.coeff <- as.numeric(param.coeff)
            distance <- param.coeff * param
            prob <- exp(-1*distance)/C
            prob
          }
)

setMethod("FindProb",
          signature=c("RankData","RankControlCayley"),
          definition<- function(dat,ctrl,modal_ranking,param){
            allperm <- AllPerms(dat@nobj)
            all_coeff <- as.numeric(FindCayley(allperm, modal_ranking))
            C <- sum(exp(-1*all_coeff * param))
            param.coeff <- FindCayley(dat@ranking, modal_ranking)
            distance <- param.coeff * param
            prob <- exp(-1*distance)/C
            prob
          }
)
