### Generate data
globalVariables(c("::", "N.sim", "filepath"))
#' Generate simulated data
#'
#' This function creates simulated correlated data. It creates both the treatment and control groups
#' with N.outcomes endpoints with certain correlation
#'
#' @param RR Params are risk ratio, incidence of outcomes, N subjects in both groups, number of outcomes and correlation among outcomes
#' @return The simulated data
#' @export
gen.data.corr <- function(RR,prop.outcome,N1, N2, N.outcomes, cor){
  sigma <- matrix(cor,N.outcomes,N.outcomes); diag(sigma)=1
  ## Simulate 10000 x-y pairs, and check that they have the specified
  ## correlation structure
  prob.Y1 <- prop.outcome
  Y1 <- bindata::rmvbin(N1, margprob = prob.Y1, sigma = sigma)
  Y2 <- bindata::rmvbin(N2, margprob = prob.Y1*RR, sigma = sigma)
  N.unvax <- nrow(Y1)
  N.vax <- nrow(Y2)
  out.data=list('Y1'=Y1,'Y2'=Y2,'N.unvax'=N.unvax,'N.vax'=N.vax)
  return(out.data)
}

###
#' Run permutation and naive functions
#'
#' This function takes the simulated data and computes the test statistics using permutation approaches
#' @param ds the simulated data and a matrix permute to indicate the permuted indeces
#' @return P-values from each method
#' @export
perm.fun <- function(ds, permute,N.outcomes){
  Y = rbind(ds$Y1, ds$Y2)
  # Y: RxC observation matrix
  # permute: NxR permutation matrix. permute[n, r] is TRUE if row r is vaxed in permutation n
  # permute %*% Y = NxC matrix. [n, c] = column sum of c column of vaxed rows of Y in permutation n
  N.vax.outcome = permute %*% Y
  N.unvax.outcome = (!permute) %*% Y
  RR <- (N.vax.outcome/ds$N.vax)/(N.unvax.outcome/ds$N.unvax)
  ### RR_varP: compute weighted perm  with 1/var as weight
  RR.log<-log(RR)
  Var_RR<-((1/N.vax.outcome+1/N.unvax.outcome)-(1/ds$N.vax+1/ds$N.unvax))
  SE_RR<-sqrt(Var_RR)
  w.j <- (1/Var_RR)/rowSums(1/Var_RR)
  RR_w_V<-rowSums(RR.log*w.j)
  ### RR_minP: take the min across RRs
  min.RR <- apply(RR,1, min)
  ### Compute test statistics over observed data
  t.obs <- test_statistics_obs(ds,N.outcomes)
  ### P_value weighted
  t.stat.V <-stats::quantile(RR_w_V,probs=0.05,na.rm = TRUE)
  sig.w.V  <- 1-mean(RR_w_V > t.obs$Obs.w.V,na.rm = TRUE)
  # P_value product for min
  #sig.RR  <- 1-colMeans(RR>t.obs$Obs.RR,na.rm = TRUE)
  #sig.RR.min <- min(sig.RR)
  sig.RR.min  <- 1-mean(min.RR>t.obs$Obs.min,na.rm = TRUE)
  #### Decision from naive
  naive<-t.obs$Obs.n
  #### Compute power for each outcome as status quo
  ind <- t.obs$Obs.ind
  out.v <- c(sig.w.V,sig.RR.min,naive)
  return(out.v)
}

###
#' Test statistics observed data
#'
#' This function computes the test statistics of the observed data
#' @param ds the simulated data
#' @return Observed test statistics
#' @export
test_statistics_obs <- function(ds,N.outcomes){
  ### Observed data
  Obs.unvax.events <- colSums(ds$Y1)
  Obs.vax.events <- colSums(ds$Y2)
  Obs.RR <- (Obs.vax.events/ds$N.vax)/(Obs.unvax.events/ds$N.unvax)
  N.vax = nrow(ds$Y2)
  N.unvax = nrow(ds$Y1)
  ### Compute test statistics 1: weighted average
  Var_RR <- (1/Obs.vax.events+1/Obs.unvax.events)-(1/N.vax+1/N.unvax)
  w.j <- (1/Var_RR)/sum(1/Var_RR)
  Obs.w.V <-sum(log(Obs.RR)*w.j)
  ### Compute test statistics 2: take the min
  Obs.min <- min(Obs.RR)
  #### Compute test statistics 3: naive approach
  naive<-RR_naive(ds,N.outcomes)
  Obs_naive<-naive$RR_m
  #### Compute test statistics for each outcome as status quo
  out.obs <- list("Obs.w.V"=Obs.w.V,"Obs.n"=Obs_naive,"Obs.min"=Obs.min)
  return(out.obs)
}

###
#' Naive method
#'
#' This function computes the p-value for each independent outcome and selects the smaller p-value
#' @param ds the simulated data
#' @return P-values from each method
#' @export
RR_naive <- function(ds,N.outcomes){
  N.unvax.outcomes <- colSums(ds$Y1)
  N.vax.outcomes <- colSums(ds$Y2)
  RR <- (N.vax.outcomes/ds$N.vax)/(N.unvax.outcomes/ds$N.unvax)
  log.RR <- log(RR)
  N.vax.prop <- (ds$N.vax-N.vax.outcomes)/N.vax.outcomes
  N.unvax.prop <- (ds$N.unvax-N.unvax.outcomes)/N.unvax.outcomes
  #CI_lower <- exp(log.RR-qnorm(1-0.05/(N.outcomes))*sqrt(N.vax.prop/ds$N.vax +N.unvax.prop/ds$N.unvax))
  CI_upper.B <- exp(log.RR+stats::qnorm(1-0.05/(N.outcomes))*sqrt(N.vax.prop/ds$N.vax +N.unvax.prop/ds$N.unvax))
  CI_upper <- exp(log.RR+stats::qnorm(1-0.05)*sqrt(N.vax.prop/ds$N.vax +N.unvax.prop/ds$N.unvax))
  RR_ind <- ifelse(CI_upper<1,1,0)
  RR_dec <- max(ifelse(CI_upper.B<1,1,0))
  result<- list("RR_ind"=RR_ind,"RR_m"=RR_dec)
  return(result)
}

###
#' Power and type I error
#'
#' This function computes the power and the type I error from each method
#' @param l.result p-values from each dataset and method
#' @return Power and type I error from each method
#' @export
compute_power<-function(l.result){
  m.result <- do.call(rbind,l.result)
  power<-rep(0,times=ncol(m.result))
  SE<-rep(0,times=ncol(m.result))
  alpha<-0.05
  power[1:2]<-colMeans(m.result[,1:2]<alpha,na.rm = TRUE)
  power[3]<-mean(m.result[,3])
  SE[1:2]<-matrixStats::colSds(1*(m.result[,1:2]<alpha),na.rm=TRUE)/sqrt(N.sim)
  SE[3] <- stats::sd(m.result[,3])/sqrt(N.sim)
  result<-rbind(power*100,SE*100)
  colnames(result)<-c("p.w.var","p.min","p.Bonf")
  return(result)
}

###
#' Main function
#'
#' This function creates the simulated data for each dataset and computes power and type I error
#' @param N.sim number of datasets, number of outcomes, RR, incidence, subjects in each category and correlation among outcomes
#' @return Final results
#' @export
main_run<-function(N.sim,N.outcomes,RR,prop.outcome,N1,N2,cor){
  sim.data.p <- pbapply::pbreplicate(N.sim, gen.data.corr(RR,prop.outcome,N1, N2, N.outcomes, cor), simplify=F)
  sim.data.t <- pbapply::pbreplicate(N.sim, gen.data.corr(RR=rep(1,length(RR)),prop.outcome,N1, N2, N.outcomes, cor), simplify=F)
  #options(future.rng.onMisuse="ignore")
  ### Repeat permutation for 999 times
  N.permute = 999
  ### Precompute permutations for better performance
  permute = plyr::raply(
    N.permute,
    # random shuffled array of N1 unvax (FALSE) and N2 vax (TRUE)
    c(rep(FALSE, N1), rep(TRUE, N2))[dqrng::dqsample.int(N1 + N2)]
  )
  d.perm.p <- plyr::llply(sim.data.p, perm.fun, permute=permute,N.outcomes=N.outcomes, .parallel=TRUE)
  d.perm.t <- plyr::llply(sim.data.t, perm.fun, permute=permute,N.outcomes=N.outcomes, .parallel=TRUE)
  r.p<-as.data.frame(compute_power(d.perm.p))
  r.t<-as.data.frame(compute_power(d.perm.t))
  r<-list("Power"=r.p,"Type I error"=r.t)
  saveRDS(r,file=filepath)
  return(r)
}

