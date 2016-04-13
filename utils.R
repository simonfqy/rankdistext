.onUnload <- function (libpath) {
  library.dynam.unload("rankdistext", libpath)
}

paramTow = function(param.true){
  w.true = rev(cumsum(rev(param.true)))
  w.true
}

wToparam = function(w.true){
  param.true = numeric(length(w.true))
  param.true[1:(length(w.true)-1)] = -diff(w.true)
  param.true[length(param.true)] = w.true[length(w.true)]
  param.true
}

UpdateCount = function(dat,count){
  dat@count = count
  dat@nobs = sum(count)
  if (length(dat@topq)!=1 && min(dat@topq) < dat@nobj-1){
    dat@subobs = numeric(length(dat@topq))
    for (i in 1:length(dat@topq)){
      dat@subobs[i] = sum(dat@count[ dat@q_ind[i]: (dat@q_ind[i+1]-1) ])
    }
  }
  dat
}



# used in SearchPi0: make a optimization result into a model
AddInfo <- function(solveres,dat,pi0){
  solveres$nobs = dat@nobs
  solveres$nobj = dat@nobj
  solveres$pi0.ranking = pi0
  solveres
}

# neighbour for incomplete rankings
SearchPi0Ext <- function(dat,init,ctrl){
  n = dat@nobj
  curr_best_ranking = init@modal_ranking.init[[1]]
  #already observed rankings.
  obsKeys <- RanktoHash(dat@rankings)
  #a hash table with boolean values signifying whether an observed ranking has been checked.
  hash::.set(obsCheckedTab, keys = obsKeys, values=replicate(dat@ndistinct, FALSE))
  if (ctrl@SearchPi0_show_message){
    message("<<< initial ranking ",curr_best_ranking," >>>")
  }
  if (max(dat@topq) < n-1){
    curr_best_ranking[curr_best_ranking>max(dat@topq)+1]=max(dat@topq)+1
  }
  print(system.time(curr_solve <- SingleClusterModel(dat,init,ctrl,curr_best_ranking)))
  currkeys = RanktoHash(curr_best_ranking)
  if( hash::has.key(currkeys, obsCheckedTab)){
    #Mark the keys corresponding to the curr_best_ranking as TRUE, or, visited.
    hash::.set(obsCheckedTab, keys = currkeys, values = TRUE)
  }
  curr_model = AddInfo(curr_solve,dat,curr_best_ranking)
  curr_criteria <- curr_solve$criteria
  #Found only means some candidates are found, not necessarily the true answer.
  found <- FALSE
  candidates <- NULL
  if(all(curr_solve$param.est > 1e-6)){
    found <- TRUE
    candidates <- c(candidates, curr_model)
  }
  #In the extended model, the SearchPi0_FUN is actually the expected log likelihood
  #for complete variables, no longer mere log likelihood.
  FUN <- ctrl@SearchPi0_FUN
  curr_goodness = FUN(curr_model)
  hashtable = hash::hash()
  hash::.set(hashtable,keys = RanktoHash(curr_best_ranking),values=TRUE)
  SearchPi0_step = 0
  nonCand <- NULL
  while(!found){
    #The boolean found is false.
    SearchPi0_step <- SearchPi0_step+1
    if (SearchPi0_step > ctrl@SearchPi0_limit){
      if (ctrl@SearchPi0_show_message){
        message("Search Pi0 limit has been reached. Stop at current best: ",
                this_ranking)
      }
      break
    }
    if (ctrl@SearchPi0_neighbour=="Cayley"){
      neighbours = CayleyNeighbour(curr_best_ranking)
    } else {
      neighbours = KendallNeighbour(curr_best_ranking)
    }

    testkeys = RanktoHash(neighbours)
    tested = hash::has.key(testkeys,hashtable)
    # A vector denoting whether a matrix of rankings are observed.
    observed <- hash::has.key(testkeys, obsCheckedTab)
    #denotes whether all observed rankings are checked.
    allChecked <- all(unname(values(obsCheckedTab)))
    #allChecked could be used to traverse the neighbors of the best observed ranking,
    #which is the only case that we visit the unobserved rankings.

    if (all(tested)){
      if (allChecked)
        break
      else{
        #All neighbours of the current best ranking are checked,
        #but there are still unchecked observed rankings left.
        #We should "jump" to an unchecked observed ranking.

      }
    }




    for (i in 1:nrow(neighbours)){
      # tested neighbours cannot be better
      #if(all(neighbours[i,] == c(1,2,3,4,5)))
      #  browser()
      if (tested[i]) next

      #when there are unchecked observed rankings, don't visit the unobserved rankings.
      if (!allChecked && !observed[i]) next
      this_ranking = neighbours[i,]
      if (ctrl@SearchPi0_show_message){
        message("\tNow Checking Neighbour ",this_ranking)
      }
      hash::.set(hashtable,keys=testkeys[i],
                 values=TRUE)
      if (observed[i]){
        hash::.set(obsCheckedTab, keys = testkeys[i],
                   values = TRUE)
      }
      this_solve <- EMSolver(dat, curr_model, ctrl, this_ranking,
                             init, curr_criteria, goDeep = FALSE)
      this_model <- AddInfo(this_solve,dat,this_ranking)
      this_criteria <- this_solve$criteria
      allChecked <- all(unname(values(obsCheckedTab)))

      if(this_criteria > curr_criteria && all(this_solve$param.est > 1e-6)){
        found <- TRUE
        candidates <- c(candidates, this_model)
      }
      else{
        nonCand <- c(nonCand, this_model)
      }

      this_goodness = FUN(this_model)
      if (found){
        #curr_goodness <- this_goodness
        #curr_best_ranking <- this_ranking
        curr_model <- this_model
        break
      }
      #if not found, go on.
      if (this_criteria > curr_criteria){
        curr_criteria <- this_criteria
        curr_best_ranking <- this_ranking
        curr_model <- this_model
        if (ctrl@SearchPi0_show_message){
          message("***Best changed to ",curr_best_ranking,"***")
        }
        if (ctrl@SearchPi0_fast_traversal)
          break
      }
    }
  }
  if (!found){
    curr_model <- EMSolver(dat, curr_model, ctrl, curr_best_ranking,
                           init = NULL, oldCrit = NULL, goDeep = TRUE )
  }
  curr_model$SearchPi0_step = SearchPi0_step
  return(curr_model)
}

#For running EM algorithm. Within the iterations, the rankings don't change, only
#other parameter values. Perhaps could be put in utils.r. It will output the parameter
#values, complete log likelihood given the specified ranking.
#It is very similar to function SingleClusterModel()
EMSolver <- function(dat, old_model, ctrl, modal_ranking, init, oldCrit, goDeep){
  #The only difference between EMSolver() and SingleClusterModel() is that in this
  #function, old_model replaces the init.

  param_len = max(dat@topq)
  paramCoeff = CWeightGivenPi(dat@ranking,modal_ranking)
  paramCoeff = matrix(paramCoeff, ncol = dat@ndistinct, byrow = TRUE)

  #oldRanking <- old_model$pi0.ranking
  #paramCoeffOld <- CWeightGivenPi(dat@ranking, oldRanking)
  #paramCoeffOld <- matrix(paramCoeffOld, ncol = dat@ndistinct, byrow = TRUE)

  #The iteration for EM algorithm.
  itr <- 0
  cores <- parallel::detectCores()

  while(TRUE){
    if (itr == 0){
      #No previously calculated parameters. Use the old model's values.
      if (any(abs(old_model$param.est) <= 1e-6 && !goDeep)){
        params <- init@param.init[[1]]
        alphaTau <- unlist(init@alpha.init)
      }
      else{
        params <- old_model$param.est
        alphaTau <- unlist(old_model$alpha.est)
      }
    }
    else{
      #Exist previously calculated parameters. Use the previous results.
      params <- param.est
      alphaTau <- unlist(alpha.est)
      #For each iteration, we need to calculate a new set of M and N, since both are
      #dependent on timesteps. Note that the word "old" is removed for paramCoeff.
    }
    product <- apply(params*paramCoeff, 2, sum)
    medium <- product + alphaTau
    listPDF <- NULL
    leng <- ncol(paramCoeff)
    prob <- numeric(length = leng)

    f <- function(mediumMono){
      function(x){
        adjustedPDFMono(x = x, mediumMono = mediumMono, alpha = alphaTau, param = params)
      }
    }
    listPDF <- parallel::mclapply(medium, f, mc.cores = cores)

    funcM <- NULL
    for (i in 1:leng){
      funcM <- c(funcM, local({ i <- i;
      function(x) {listPDF[[i]](x) * (-x)}
      }))
    }
    funcN <- NULL
    for (i in 1:leng){
      funcN <- c(funcN, local({ i <- i;
      function(x) {listPDF[[i]](x) * (log(x))}
      }))
    }
    f1 <- function(f, xmin, xmax, reltol){
      pracma::integral(f, xmin = xmin, xmax = xmax, reltol = reltol)
    }

    M <- unlist(parallel::mclapply(funcM, f1, xmin = 0, xmax = Inf, reltol = 3e-7,
                                   mc.cores=cores), use.names=FALSE)

    NVector <- unlist(parallel::mclapply(funcN, f1, xmin = 0, xmax = Inf, reltol = 3e-7,
                                         mc.cores = cores), use.names=FALSE)

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
        #  reference <<- a

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

      #obj() calculates the opposite of complete log likelihood of l(pi1, pi2, ..., pin)

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
      if (itr < 2)
        fctr <- 6e9
      else
        fctr <- 5e8

      opt_res <- optim(
        par = c(params, alphaTau), fn = obj, gr = gradient, method = "L-BFGS-B",
        lower = c(rep(0, length(params)), 1e-8), upper = rep(Inf, length(params) + 1),
        control = list(factr = fctr)
      )

      print("values of paramsTau:")
      print(params)
      print("values of alphaTau:")
      print(alphaTau)

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
      if (!goDeep && crit < oldCrit){
        print("No need to proceed further.")
        break
      }
    }
    if (!goDeep && any(abs(param.est) <= 1e-6)){
      print("Detected phi parameters close to 0.")
      break
    }

    if(itr > 0){
      if (log_likelihood - old_likelihood < abs(old_likelihood)*3e-6){
        #The EM algorithm failed to ascend
        print(paste("failed to ascend at iteration ", itr))
        #Use the old parameter values and likelihoods and jump out of the loop
        #param.est <- params
        #alpha.est <- alphaTau
        break
      }

      #if (itr %% 25==0)
       # browser()

      #otherwise we accept the new model and carry on to the next iteration.
      #old_goodness <- log_likelihood
    }
    old_likelihood <- log_likelihood

    itr <- itr+1
    print(itr)
    print(c(opt_res[[1]], log_likelihood))

    if (itr > ctrl@EM_limit){
      message("Algorithm did not converge in ",ctrl@EM_limit," iterations in
                  function EMSolver for modal ranking", modal_ranking)
      break
    }
  }
  #param.est = c(param.est,rep(0,dat@nobj - 1 - param_len))
  list(
    param.est = param.est, w.est = paramTow(param.est),
    alpha.est = alpha.est, log_likelihood = log_likelihood,
    criteria = crit
  )
}

# TODO need to change: does not work for d==1
t.gen <- function(d){
  t.lst = list()
  t.lst[[d]] = matrix(rep(1:d,d),ncol = d, nrow = d,byrow=T)
  left.mask = matrix(rep(0,d^2),ncol = d, nrow = d)
  left.mask[2:d,1:(d-1)] = diag(rep(1,d-1))
  t.lst[[d]][upper.tri(left.mask)] = 0
  for ( i in 1:(d-1)){
    t.lst[[d-i]] = left.mask%*%t.lst[[d-i+1]]
    diag(t.lst[[d-i]]) = c(rep(0,i),1:(d-i))
    t.lst[[d-i]][upper.tri(left.mask)] = 0
  }
  t.lst
}

GHC <- function(param,t.lst){
  d = length(param) # d = t - 1
  K = matrix(rep(0,d^2),ncol = d, nrow = d)
  for ( i in 1:d){
    K = -1 * param[i] * t.lst[[i]] + K
  }
  K = exp(K)
  K[upper.tri(K)] = 0
  gradient = numeric(d)
  ones = rep(1,d)
  denom = rowSums(K) + ones
  B = matrix(ncol=d,nrow=d)
  for (i in 1:d){
    B[,i] = rowSums(-1 * K * t.lst[[i]])
  }
  for ( i in 1:d){
    gradient[i] = sum(B[,i] / denom)
  }
  gradient
}

GHCPre <- function(param,t.lst){
  d = length(param) # d = t - 1
  K = matrix(rep(0,d^2),ncol = d, nrow = d)
  for ( i in 1:d){
    K = -1 * param[i] * t.lst[[i]] + K
  }
  K = exp(K)
  K[upper.tri(K)] = 0
  ones = rep(1,d)
  list(K, ones)
}

GHCExt <- function(param, t.lst, listPDF, f1, cores){

  grads <- function(x, K, ones, j){
    d = length(param) # d = t - 1
    K = K^x
    #browser()
    denom = rowSums(K) + ones
    #The following enclosed body is not supposed to be used. I am just being thrifty
    #To throw it away.
    if (FALSE){
      B = matrix(ncol=d,nrow=d)
      gradient <- numeric(length = d)
      for (i in 1:d){
        B[,i] = rowSums(-1 * x * K * t.lst[[i]])
      }
      for ( i in 1:d){
        gradient[i] = sum(B[,i] / denom)
      }
    }
    B <- vector(length = d)
    gradient <- numeric()
    B <- rowSums(-1 * x * K * t.lst[[j]])
    gradient <- sum(B/denom)
    gradient
  }

  #Now I know how to construct a matrix of functions. We could try

  preliminaries <- GHCPre(param, t.lst)
  K <- preliminaries[[1]]
  ones <- preliminaries[[2]]
  #Now we should start constructing the list of gradient functions.
  #Perhaps the following loop or even the whole data structure can be optimized.
  integrand <- NULL
  for (i in 1:length(listPDF)){
    for (j in 1:length(param)){
      integrand <- c(integrand, local({
        i <- i; j <- j;
        function(x){
          grads(x = x, K = K, ones = ones, j = j)*listPDF[[i]](x)
        }
      }))
    }
  }

  output <- matrix(0, nrow = length(param), ncol = length(listPDF))
  lengrow <- length(param)
  lengcol <- length(listPDF)
  #for (i in 1:lengcol){
  #  for (j in 1:lengrow){
      #output[j, i] <- pracma::integral(integrand[[(i-1)*lengrow + j]],xmin = 0, xmax = Inf, reltol = 1e-7)
      #browser()
    #}
  #}
  output1 <- unlist(parallel::mclapply(integrand, f1, xmin = 0,
                                       xmax = Inf, reltol=1e-7, mc.cores = cores),
                    use.names = FALSE)
  for (i in 1:lengcol){
    for(j in 1:lengrow){
      output[j, i] <- output1[[(i-1)*lengrow + j]]
    }
  }

  output
}

#PDF of the vector of lambda, after normalization.
adjustedPDFMono <- function (x, mediumMono, alpha, param){
  #Note that medium is a vector with length equals col number of paramCoeff.

  singleUnadPDF <- function(x, mediumMono, alpha, param){
    dgamma(x, shape = alpha, rate = mediumMono)/exp(LogC(x*param))
  }

  integra <- pracma::integral(singleUnadPDF, xmin = 0, xmax = Inf,
                              reltol = 5e-7, mediumMono = mediumMono,
                              alpha = alpha, param = param)
  output <- singleUnadPDF(x, mediumMono, alpha, param)/integra
  output
}

#This function gets the marginal conditional distribution of pi(i), by integrating out
#lambda(i).
getMarginalPiMono <- function(mediumMono, alpha, param, reltol){
  singlePDF <- function(x, mediumMono, alpha, param){
    dgamma(x, shape = alpha, rate = mediumMono)*((alpha/mediumMono)^alpha)/exp(LogC(x*param))
  }
  output <- pracma::integral(singlePDF, xmin = 0, xmax = Inf, reltol = reltol,
                             mediumMono = mediumMono, alpha = alpha, param = param)
  output
}

# find the weighted kendall distance between p1 and p2
# p1 and p2 are orderings
KwDist <- function(p1, p2,w){
  n = length(p1)
  distance = 0
  for (i in p2){
    pos1 = which(p1 == i)
    pos2 = which(p2 == i)
    relative_pos1 = (1:n - pos1)[order(p1)]
    relative_pos2 = (1:n - pos2)[order(p2)]
    Ji = which(relative_pos1 * relative_pos2 < 0)
    Ii = length(Ji)
    Li = (pos1 + pos2 + Ii)/2
    c1 = ifelse(pos1<=(Li-1), sum(w[pos1:(Li-1)]),0)
    c2 = ifelse(pos2<=(Li-1), sum(w[pos2:(Li-1)]),0)
    distance = distance + (c1 + c2)/2
  }
  distance
}

BreakTieEqualProb <- function(dat){
  ind_comp <- which(dat@topq == dat@nobj-1)
  if (max(dat@topq) == dat@nobj-1){
    ind_comp_start <- dat@q_ind[ind_comp]
    ind_comp_end <- dat@q_ind[ind_comp + 1] - 1
    comp_ranking <- dat@ranking[ind_comp_start:ind_comp_end, ]
    comp_count <- dat@count[ind_comp_start:ind_comp_end]
  } else {
    comp_ranking <- permute::allPerms(dat@nobj)
    comp_ranking <- rbind(1:dat@nobj, comp_ranking)
    comp_count <- rep(0, nrow(comp_ranking))
  }
  comp_hash <- RanktoHash(comp_ranking)

  for (i in 1:length(dat@topq)){
    if(dat@topq[i] == dat@nobj-1)
      next
    ind_start <- dat@q_ind[i]
    ind_end <- dat@q_ind[i+1] - 1
    this_q <- dat@topq[i]
    this_inc <- 1/factorial(dat@nobj - this_q)
    # generate permutations for tied group
    tie_perm <- permute::allPerms((this_q+1):dat@nobj) + this_q
    tie_perm <- rbind(tie_perm, (this_q+1):dat@nobj)
    # iterate through top-q rankings
    for (this_partial_ind in ind_start:ind_end){
      this_partial <- dat@ranking[this_partial_ind, ]
      this_count <- dat@count[this_partial_ind]
      ind_tie <- which(this_partial == this_q + 1)
      # iterate through possible tie-breakings
      for (ind_break in 1:nrow(tie_perm)){
        copy_partial <- this_partial
        this_break <- tie_perm[ind_break, ]
        ptr_break <- 1
        # iterate through tied positions
        for (this_tie_ind in ind_tie){
          copy_partial[this_tie_ind] = this_break[ptr_break]
          ptr_break <- ptr_break + 1
        }
        this_hash <- rankdist::RanktoHash(copy_partial)
        ind_incre <- which(comp_hash == this_hash)
        comp_count[ind_incre] = comp_count[ind_incre] + this_inc*this_count
      }
    }
  }
  # handle complete rankings
  ind_nonempty_count = which(comp_count != 0)
  comp_count = comp_count[comp_count > 0]
  comp_ranking = comp_ranking[ind_nonempty_count, ]
  comp_dat <- new("RankData", ranking=comp_ranking, count=comp_count)
}

