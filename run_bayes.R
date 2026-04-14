.libPaths(c("C:/R_Cloud_Fix", "C:/Program Files/R/R-4.5.2/library"))
library(brms)
library(dplyr)

N_REP      <- 200
CHECKPOINT <- "C:/simulation_lmm/sim_bayes_checkpoint.rds"
OUTPUT     <- "C:/simulation_lmm/sim_bayes_results.rds"

grid_bayes <- expand.grid(
  J=c(3,5), n=c(50,150), icc=c(0.05,0.30), beta1=c(0,0.5),
  balance="balanced", error_dist="normal", stringsAsFactors=FALSE)
grid_bayes$n_sim        <- as.character(grid_bayes$n)
grid_bayes$condition_id <- seq_len(nrow(grid_bayes))

if (file.exists(CHECKPOINT)) {
  prev <- readRDS(CHECKPOINT)
  done_ids <- unique(prev$condition_id)
  cat("Kaldigi yerden devam:", length(done_ids), "kosul tamamlanmis
")
} else {
  prev <- NULL; done_ids <- integer(0)
  cat("Yeni Bayesian simulasyon basliyor
")
}

grid_todo <- grid_bayes[!grid_bayes$condition_id %in% done_ids, ]
cat("Kalan kosul:", nrow(grid_todo), "/", nrow(grid_bayes), "
")

prior_w <- prior(cauchy(0,5), class="sd") + prior(normal(0,10), class="b")
prior_m <- prior(normal(0,1), class="sd") + prior(normal(0,2.5), class="b")

set.seed(2025)
start_time <- proc.time()

for (ci in seq_len(nrow(grid_todo))) {
  cond <- grid_todo[ci,]
  cat(sprintf("[%s] Kosul %d/%d — J=%d n=%s ICC=%.2f beta1=%.1f
",
    format(Sys.time(),"%H:%M"), ci, nrow(grid_todo),
    cond$J, cond$n_sim, cond$icc, cond$beta1))

  results <- vector("list", N_REP)

  for (r in seq_len(N_REP)) {
    J<-cond$J; n<-as.integer(cond$n_sim); icc<-cond$icc; b1<-cond$beta1
    sigma2_u<-icc; sigma2_e<-1-icc
    nj<-rep(n,J); N<-sum(nj); cid<-rep(seq_len(J),times=nj)
    uj<-rnorm(J,0,sqrt(sigma2_u))
    eij<-rnorm(N)*sqrt(sigma2_e)
    Xij<-rnorm(N)
    dat<-data.frame(Y=b1*Xij+uj[cid]+eij, X=Xij, cluster=factor(cid))
    alpha<-0.05; row_list<-list()

    fit_w <- tryCatch(
      brm(Y~X+(1|cluster),data=dat,prior=prior_w,
          chains=4,iter=2000,warmup=1000,cores=1,refresh=0,silent=2),
      error=function(e) NULL)
    if (!is.null(fit_w)) {
      post<-as.data.frame(fit_w)$b_X
      ci_q<-quantile(post,c(alpha/2,1-alpha/2))
      rhat<-tryCatch(brms::rhat(fit_w)["b_X"],error=function(e) NA)
      row_list$bw<-data.frame(method="Bayes-W",estimate=median(post),
        se=sd(post),pval=NA_real_,ci_lo=unname(ci_q[1]),ci_hi=unname(ci_q[2]),
        rejected=as.integer(ci_q[1]>0|ci_q[2]<0),
        covered=as.integer(ci_q[1]<=b1&b1<=ci_q[2]),
        rhat=unname(rhat),converged=as.integer(!is.na(rhat)&&rhat<1.05),
        rep_id=r,condition_id=cond$condition_id,J=cond$J,n=cond$n,
        icc=cond$icc,balance=cond$balance,beta1=cond$beta1,
        error_dist=cond$error_dist)
    }

    fit_m <- tryCatch(
      brm(Y~X+(1|cluster),data=dat,prior=prior_m,
          chains=4,iter=2000,warmup=1000,cores=1,refresh=0,silent=2),
      error=function(e) NULL)
    if (!is.null(fit_m)) {
      post<-as.data.frame(fit_m)$b_X
      ci_q<-quantile(post,c(alpha/2,1-alpha/2))
      rhat<-tryCatch(brms::rhat(fit_m)["b_X"],error=function(e) NA)
      row_list$bm<-data.frame(method="Bayes-M",estimate=median(post),
        se=sd(post),pval=NA_real_,ci_lo=unname(ci_q[1]),ci_hi=unname(ci_q[2]),
        rejected=as.integer(ci_q[1]>0|ci_q[2]<0),
        covered=as.integer(ci_q[1]<=b1&b1<=ci_q[2]),
        rhat=unname(rhat),converged=as.integer(!is.na(rhat)&&rhat<1.05),
        rep_id=r,condition_id=cond$condition_id,J=cond$J,n=cond$n,
        icc=cond$icc,balance=cond$balance,beta1=cond$beta1,
        error_dist=cond$error_dist)
    }

    results[[r]] <- do.call(rbind, row_list)
    if (r%%10==0) cat(sprintf("  Rep %d/%d
",r,N_REP))
  }

  cond_res <- do.call(rbind, Filter(Negate(is.null), results))
  if (file.exists(CHECKPOINT)) {
    updated <- rbind(readRDS(CHECKPOINT), cond_res)
  } else { updated <- cond_res }
  saveRDS(updated, CHECKPOINT)
  rm(cond_res, updated); gc()

  elapsed <- round((proc.time()-start_time)["elapsed"]/60,1)
  remaining <- round((elapsed/ci)*(nrow(grid_todo)-ci),1)
  cat(sprintf("  Tamamlandi. Gecen: %.1f dk | Kalan: %.1f dk
",elapsed,remaining))
}

final_bayes <- readRDS(CHECKPOINT)
saveRDS(final_bayes, OUTPUT)
cat("BAYESIAN TAMAMLANDI! Toplam satir:", nrow(final_bayes), "
")
cat("Sure:", round((proc.time()-start_time)["elapsed"]/3600,2), "saat
")
