### Generate data

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
  ### Compute weighted perm
  w.j <- (N.unvax.outcome)/rowSums(N.unvax.outcome)
  RR_w <- w.j*RR
  RR_w_s <-rowSums(RR_w)
  ### compute RR_prod/SE
  RR.log<-log(RR)
  Var_RR<-((1/N.vax.outcome+1/N.unvax.outcome)-(1/ds$N.vax+1/ds$N.unvax))
  SE_RR<-sqrt(Var_RR)
  RR_p_SE<-rowSums(RR.log/SE_RR)
  ### Comoute test statistics over observed data
  t.obs <- test_statistics_obs(ds,N.outcomes)
  ### P_value weighted
  sig.w.RR  <- 1-mean(RR_w_s > t.obs$Obs.w,na.rm = TRUE)
  t.stat.w <- quantile(RR_w_s,probs=0.05,na.rm = TRUE)
  # P_value product/SE
  sig.p.RR.SE  <- 1-mean( RR_p_SE>t.obs$Obs.p,na.rm = TRUE)
  t.stat.p <-quantile(RR_p_SE,probs=0.05,na.rm = TRUE)
  #### Decision from naive
  naive<-t.obs$Obs.n
  #### Compute power for each outcome as status quo
  ind <- t.obs$Obs.ind
  out.v <- c(t.stat.w,t.stat.p,sig.w.RR,sig.p.RR.SE,naive,ind)
  return(out.v)
}

###
#' Test statistics observed data
#'
#' This function computes the test statistcs of the observed data
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
  w_j <- colSums(ds$Y1)/sum(colSums(ds$Y1))
  Obs.w.RR <- sum(Obs.RR*w_j)
  # Compute test statistics 2: product of RRs over SE
  Var_RR <- (1/Obs.vax.events+1/Obs.unvax.events)-(1/N.vax+1/N.unvax)
  SE_RR <- sqrt(Var_RR)
  Obs.p.RR.SE<-sum(log(Obs.RR)/SE_RR)
  #### Compute test statistics 3: naive approach
  naive<-RR_naive(ds,N.outcomes)
  Obs_naive<-naive$RR_m
  #### Compute test statistics for each outcome independently as status quo
  Obs_ind <- naive$RR_ind
  out.obs <- list("Obs.w"=Obs.w.RR,"Obs.p"=Obs.p.RR.SE,"Obs.n"=Obs_naive,"Obs.ind"=Obs_ind)
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
  CI_upper.B <- exp(log.RR+qnorm(1-0.05/(N.outcomes))*sqrt(N.vax.prop/ds$N.vax +N.unvax.prop/ds$N.unvax))
  CI_upper <- exp(log.RR+qnorm(1-0.05)*sqrt(N.vax.prop/ds$N.vax +N.unvax.prop/ds$N.unvax))
  RR_ind <- ifelse(CI_upper<1,1,0)
  RR_dec <- max(ifelse(CI_upper.B<1,1,0))
  result<- list("RR_ind"=RR_ind,"RR_m"=RR_dec)
  return(result)
}

###
#' Naive method for observed data
#'
#' This function computes the p-value for each independent outcome and selects the smaller p-value using observed data
#' @param ds the simulated data
#' @return P-values from each method
#' @export
RR_naive_obs <- function(ds){
  N.unvax.outcomes <- colSums(ds$Y1)
  N.vax.outcomes <- colSums(ds$Y2)
  RR <- (N.vax.outcomes/ds$N.vax)/(N.unvax.outcomes/ds$N.unvax)
  log.RR <- log(RR)
  N.vax.prop <- (ds$N.vax-N.vax.outcomes)/N.vax.outcomes
  N.unvax.prop <- (ds$N.unvax-N.unvax.outcomes)/N.unvax.outcomes
  CI_lower.B <- exp(log.RR-qnorm(1-0.05/(N.outcomes))*sqrt(N.vax.prop/ds$N.vax +N.unvax.prop/ds$N.unvax))
  CI_upper.B <- exp(log.RR+qnorm(1-0.05/(N.outcomes))*sqrt(N.vax.prop/ds$N.vax +N.unvax.prop/ds$N.unvax))
  CI_lower<-exp(log.RR-qnorm(1-0.05)*sqrt(N.vax.prop/ds$N.vax +N.unvax.prop/ds$N.unvax))
  CI_upper<-exp(log.RR+qnorm(1-0.05)*sqrt(N.vax.prop/ds$N.vax +N.unvax.prop/ds$N.unvax))
  diff.B<-CI_upper.B-CI_lower.B
  SE.B<- diff.B/1.96
  z.B<-RR/SE.B
  p.value.B<-min(pnorm(-abs(z.B)))
  diff.ind<-CI_upper-CI_lower
  SE.ind<- diff.ind/1.96
  z.ind<-RR/SE.ind
  p.value.ind<-pnorm(-abs(z.ind))
  result<- c(p.value.B,p.value.ind)
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
    m.result<-m.result[,3:ncol(m.result)]
    power<-rep(0,times=ncol(m.result))
    SE<-rep(0,times=ncol(m.result))
    alpha<-0.05
    power[1:2]<-colMeans(m.result[,1:2]<alpha,na.rm = TRUE)
    power[3:length(power)]<-colMeans(m.result[,3:ncol(m.result)],na.rm = TRUE)
    SE[1:2]<-matrixStats::colSds(1*(m.result[,1:2]<alpha),na.rm=TRUE)/sqrt(N.sim)
    SE[3:length(power)] <- matrixStats::colSds(m.result[,3:ncol(m.result)],na.rm = TRUE)/sqrt(N.sim)
    result<-rbind(power*100,SE*100)
    colnames(result)<-c("perm.w","perm.p.SE","naive.RR",paste0("ind_",1:N.outcomes))
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
  options(future.rng.onMisuse="ignore")
  ### Repeat permutation for 999 times
  N.permute = 999
  ### Precompute permutations for better performance
  permute = plyr::raply(
    N.permute,
    # random shuffled array of N1 unvax (FALSE) and N2 vax (TRUE)
    c(rep(FALSE, N1), rep(TRUE, N2))[dqrng::dqsample.int(N1 + N2)]
  )
  d.perm.p <- plyr::llply(sim.data.p, perm.fun, permute=permute,N.outcomes=N.outcomes,.parallel=TRUE)
  d.perm.t <- plyr::llply(sim.data.t, perm.fun, permute=permute,N.outcomes=N.outcomes,.parallel=TRUE)
  r.p<-as.data.frame(compute_power(d.perm.p))
  r.t<-as.data.frame(compute_power(d.perm.t))
  r<-list("Power"=r.p,"Type I error"=r.t)
  saveRDS(r,file=filepath)
  return(r)
}

###
#' P-value real data
#'
#' This function computes the p-value for real-world data using all methods
#' @param ds the real-world data and the matrix with the permuted indeces specified
#' @return P-values from each method
#' @export
compute_pvalue_obs <- function(ds,permute){
  p.value.perm<-perm.fun(ds,permute)[1:4]
  p.value.naive<-RR_naive_obs(ds)
  p.value<-c(p.value.perm,p.value.naive)
  p.value<-setNames(p.value,c("t.stat.w","t.stat.p","pvalue.w","pvalue.p.SE","pvalue.naive",paste0("pvalue.ind_",1:N.outcomes)))
  return(p.value)
}

###
#' Boostrapping method
#'
#' This function computes the bootstrapping for each permutation method
#' @param ds the real-world data
#' @return Bootstrapped test statistics
#' @export
test_statistics_boot <- function(ds){
  ### Observed data
  Y1<-as.data.frame(matrix(NA,ncol=ncol(ds$Y1),nrow(ds$Y1)))
  for(l in 1:nrow(ds$Y1)){
    Y1[l,]<- sample(ds$Y1[l,],  replace=TRUE)
  }
  Y2<-as.data.frame(matrix(NA,ncol=ncol(ds$Y2),nrow(ds$Y2)))
  for(l in 1:nrow(ds$Y2)){
    Y2[l,]<- sample(ds$Y2[l,],  replace=TRUE)
  }
  Obs.unvax.events <- colSums(Y1)
  Obs.vax.events <- colSums(Y2)
  Obs.RR <- (Obs.vax.events/ds$N.vax)/(Obs.unvax.events/ds$N.unvax)
  N.vax = nrow(Y2)
  N.unvax = nrow(Y1)
  ### Compute test statistics 1: weighted average
  w_j <- colSums(Y1)/sum(colSums(Y1))
  Obs.w.RR <- sum(Obs.RR*w_j)
  # Compute test statistics 2: product of RRs over SE
  Var_RR <- (1/Obs.vax.events+1/Obs.unvax.events)-(1/N.vax+1/N.unvax)
  SE_RR <- sqrt(Var_RR)
  Obs.p.RR.SE<-sum(log(Obs.RR)/SE_RR)
  out.obs<-c(Obs.w.RR,Obs.p.RR.SE)
  names(out.obs) <- c("Obs.w","Obs.p")
  return(out.obs)
}

###
#' Uncertainty for bootstrapping
#'
#' This function computes the uncertainty around the bootstrapped test statistics
#' @param ds the simulated data
#' @return P-values from each method
#' @export
compute_uncertainty_boot <- function(ds){
  boot.ci<-apply(ds,1,function(x) quantile(x,probs=c(0.50,0.025,0.975)))
  return(boot.ci)
}

###
#' Load results from permutation + naive methods
#'
#' This function extracts the results from naive and permutation methods
#' @param pattern the file path
#' @return List of the results
#' @export
load_result <- function(pattern){
  file_list <- list.files(path="~/Documents/RSVtrial_simulation/results/", pattern=(pattern))
  test1 <-
    all.res <-
    pbapply::pblapply(file_list, function(x){
      print(x)
      path1 <- paste0("~/Documents/RSVtrial_simulation/results/",x)
      test1 <- readRDS(path1)
      })
  power <- dplyr::bind_rows(sapply(test1,`[`, "Power")) #combine power
  power <- power[,1:3]
  power<-stack(power)
  power.se <- power %>% filter(row_number() %% 2 == 0) ## Select even rows
  power.e <- power %>% filter(row_number() %% 2 == 1) ## Select odd rows
  typeone <- dplyr::bind_rows(sapply(test1,`[`, "Type I error")) #combine type1 error
  typeone <- typeone[,1:3]
  typeone<-stack(typeone)
  typeone.se <- typeone %>% filter(row_number() %% 2 == 0) ## Select even rows
  typeone.e <- typeone %>% filter(row_number() %% 2 == 1) ## Select odd rows
  power.l <- power.e
  power.l$values <- power.e$values-1.96*power.se$values
  power.u <- power.e
  power.u$values <- power.e$values+1.96*power.se$values
  typeone.l <-typeone.e
  typeone.l$values <- typeone.e$values-1.96*typeone.se$values
  typeone.u <- typeone.e
  typeone.u$values <- typeone.e$values+1.96*typeone.se$values
  power.df<-as.data.frame(cbind(power.e,power.l$values,power.u$values))
  typeone.df<-as.data.frame(cbind(typeone.e,typeone.l$values,typeone.u$values))
  colnames(power.df)<-c("estimate","model","lower","upper")
  colnames(typeone.df)<-c("estimate","model","lower","upper")
  result<-list("power.df"=power.df,"typeone.df"=typeone.df)
  return(result)
}

###
#' Load results from independent outcomes
#'
#' This function extracts the results from outcomes treated as independent
#' @param pattern the file path
#' @return List of the results
#' @export
load_result_ind <- function(pattern){
  file_list <- list.files(path="~/Documents/RSVtrial_simulation/results/", pattern=(pattern))
  test1 <-
    all.res <-
    pbapply::pblapply(file_list, function(x){
      print(x)
      path1 <- paste0("~/Documents/RSVtrial_simulation/results/",x)
      test1 <- readRDS(path1)
    })
  power <- dplyr::bind_rows(sapply(test1,`[`, "Power")) #combine power
  power <- power[,4:ncol(power)]
  power<-stack(power)
  power.se <- power %>% filter(row_number() %% 2 == 0) ## Select even rows
  power.e <- power %>% filter(row_number() %% 2 == 1) ## Select odd rows
  typeone <- dplyr::bind_rows(sapply(test1,`[`, "Type I error")) #combine type1 error
  typeone <- typeone[,4:ncol(typeone)]
  typeone<-stack(typeone)
  typeone.se <- typeone %>% filter(row_number() %% 2 == 0) ## Select even rows
  typeone.e <- typeone %>% filter(row_number() %% 2 == 1) ## Select odd rows
  power.l <- power.e
  power.l$values <- power.e$values-1.96*power.se$values
  power.u <- power.e
  power.u$values <- power.e$values+1.96*power.se$values
  typeone.l <-typeone.e
  typeone.l$values <- typeone.e$values-1.96*typeone.se$values
  typeone.u <- typeone.e
  typeone.u$values <- typeone.e$values+1.96*typeone.se$values
  power.df<-as.data.frame(cbind(power.e,power.l$values,power.u$values))
  typeone.df<-as.data.frame(cbind(typeone.e,typeone.l$values,typeone.u$values))
  colnames(power.df)<-c("estimate","outcome","lower","upper")
  colnames(typeone.df)<-c("estimate","outcome","lower","upper")
  result<-list("power.df"=power.df,"typeone.df"=typeone.df)
  return(result)
}

