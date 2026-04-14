
.libPaths(c("C:/R_Cloud_Fix", "C:/Program Files/R/R-4.5.2/library"))

library(furrr); library(future); library(dplyr)
plan(multisession, workers = 3)

N_REP      <- 1000
CHECKPOINT <- "C:/simulation_lmm/sim_freq_checkpoint.rds"
OUTPUT     <- "C:/simulation_lmm/sim_freq_results.rds"

sim_grid <- expand.grid(
  J=c(3,5,7,10,15), n=c(20,50,150),
  icc=c(0.05,0.15,0.30,0.50),
  balance=c("balanced","unbalanced"),
  beta1=c(0,0.3,0.5), error_dist="normal",
  stringsAsFactors=FALSE)
sim_grid$n_sim <- ifelse(sim_grid$balance=="unbalanced","unbalanced",as.character(sim_grid$n))
grid_freq <- sim_grid
grid_freq$condition_id <- seq_len(nrow(grid_freq))

prev     <- readRDS(CHECKPOINT)
done_ids <- unique(prev$condition_id)
cat("Kaldigi yerden devam:", length(done_ids), "kosul tamamlanmis
")

grid_todo <- grid_freq[!grid_freq$condition_id %in% done_ids, ]
cat("Kalan kosul:", nrow(grid_todo), "
")

batch_size <- 10
batches    <- split(seq_len(nrow(grid_todo)),
                    ceiling(seq_len(nrow(grid_todo))/batch_size))

set.seed(2025)
start_time <- proc.time()

for (b in seq_along(batches)) {
  bidx       <- batches[[b]]
  batch_grid <- grid_todo[bidx, ]
  cat(sprintf("[%s] Batch %d/%d — kosul %d:%d
",
              format(Sys.time(),"%H:%M"), b, length(batches),
              batch_grid$condition_id[1],
              batch_grid$condition_id[nrow(batch_grid)]))

  batch_res <- future_map_dfr(
    seq_len(nrow(batch_grid)),
    function(i) {
      .libPaths(c("C:/R_Cloud_Fix","C:/Program Files/R/R-4.5.2/library"))
      library(lme4); library(lmerTest)
      cond  <- batch_grid[i,]
      alpha <- 0.05
      n_rep <- N_REP
      n_val <- if(cond$balance=="unbalanced") "unbalanced" else as.integer(cond$n_sim)
      results <- vector("list", n_rep)
      for (r in seq_len(n_rep)) {
        J<-cond$J; icc<-cond$icc; b1<-cond$beta1; edist<-cond$error_dist
        sigma2_u<-icc; sigma2_e<-1-icc
        nj<-if(identical(n_val,"unbalanced")) sample(50:150,J,replace=TRUE) else rep(n_val,J)
        N<-sum(nj); cid<-rep(seq_len(J),times=nj)
        uj<-rnorm(J,0,sqrt(sigma2_u))
        estd<-switch(edist,
          "normal"=rnorm(N),
          "skewed"={x<-rchisq(N,3);(x-3)/sqrt(6)},
          "heavy"={x<-rt(N,3);x/sqrt(3)})
        Xij<-rnorm(N); eij<-estd*sqrt(sigma2_e)
        dat<-data.frame(Y=b1*Xij+uj[cid]+eij,X=Xij,cluster=factor(cid))
        res_list<-list()
        fit<-tryCatch(lmerTest::lmer(Y~X+(1|cluster),data=dat,REML=TRUE),error=function(e)NULL)
        if(!is.null(fit)){
          cf<-summary(fit,ddf="Satterthwaite")$coefficients
          b1e<-cf["X","Estimate"]; se1<-cf["X","Std. Error"]
          pv<-2*pnorm(abs(b1e/se1),lower.tail=FALSE)
          ci<-b1e+c(-1,1)*qnorm(1-alpha/2)*se1
          res_list$reml<-data.frame(method="REML",estimate=b1e,se=se1,pval=pv,
            ci_lo=ci[1],ci_hi=ci[2],rejected=as.integer(pv<alpha),
            covered=as.integer(ci[1]<=b1&b1<=ci[2]))}
        fit_s<-tryCatch(lmerTest::lmer(Y~X+(1|cluster),data=dat,REML=TRUE),error=function(e)NULL)
        if(!is.null(fit_s)){
          cf<-summary(fit_s,ddf="Satterthwaite")$coefficients
          b1e<-cf["X","Estimate"]; se1<-cf["X","Std. Error"]
          df0<-cf["X","df"]; pv<-cf["X","Pr(>|t|)"]
          ci<-b1e+c(-1,1)*qt(1-alpha/2,df0)*se1
          res_list$satt<-data.frame(method="Satterthwaite",estimate=b1e,se=se1,pval=pv,
            ci_lo=ci[1],ci_hi=ci[2],rejected=as.integer(pv<alpha),
            covered=as.integer(ci[1]<=b1&b1<=ci[2]))}
        fit_k<-tryCatch(lmerTest::lmer(Y~X+(1|cluster),data=dat,REML=TRUE),error=function(e)NULL)
        if(!is.null(fit_k)){
          cf<-summary(fit_k,ddf="Kenward-Roger")$coefficients
          b1e<-cf["X","Estimate"]; se1<-cf["X","Std. Error"]
          df0<-cf["X","df"]; pv<-cf["X","Pr(>|t|)"]
          ci<-b1e+c(-1,1)*qt(1-alpha/2,df0)*se1
          res_list$kr<-data.frame(method="KR",estimate=b1e,se=se1,pval=pv,
            ci_lo=ci[1],ci_hi=ci[2],rejected=as.integer(pv<alpha),
            covered=as.integer(ci[1]<=b1&b1<=ci[2]))}
        rep_res<-do.call(rbind,res_list)
        if(!is.null(rep_res)){
          rep_res$rep_id<-r; rep_res$condition_id<-cond$condition_id
          rep_res$J<-cond$J; rep_res$n<-cond$n; rep_res$icc<-cond$icc
          rep_res$balance<-cond$balance; rep_res$beta1<-cond$beta1
          rep_res$error_dist<-cond$error_dist}
        results[[r]]<-rep_res}
      do.call(rbind,Filter(Negate(is.null),results))},
    .options=furrr_options(seed=TRUE))

  prev_ckpt <- readRDS(CHECKPOINT)
  updated   <- rbind(prev_ckpt, batch_res)
  saveRDS(updated, CHECKPOINT)
  rm(batch_res, prev_ckpt, updated); gc()

  elapsed   <- round((proc.time()-start_time)["elapsed"]/60,1)
  remaining <- round((elapsed/b)*(length(batches)-b),1)
  cat(sprintf("  Tamamlandi. Gecen: %.1f dk | Kalan: %.1f dk
",elapsed,remaining))
}

final_freq <- readRDS(CHECKPOINT)
saveRDS(final_freq, OUTPUT)
cat("TAMAMLANDI! Toplam satir:", nrow(final_freq), "
")
cat("Sure:", round((proc.time()-start_time)["elapsed"]/3600,2), "saat
")

