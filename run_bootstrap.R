
.libPaths(c("C:/R_Cloud_Fix", "C:/Program Files/R/R-4.5.2/library"))

library(furrr); library(future); library(dplyr)
plan(multisession, workers = 3)

N_REP      <- 1000
NSIM_BOOT  <- 99
CHECKPOINT <- "C:/simulation_lmm/sim_boot_checkpoint.rds"
OUTPUT     <- "C:/simulation_lmm/sim_boot_results.rds"

sim_grid <- expand.grid(
  J=c(3,5,7,10,15), n=c(20,50,150),
  icc=c(0.05,0.15,0.30,0.50),
  balance=c("balanced","unbalanced"),
  beta1=c(0,0.3,0.5), error_dist="normal",
  stringsAsFactors=FALSE)
sim_grid$n_sim <- ifelse(sim_grid$balance=="unbalanced","unbalanced",as.character(sim_grid$n))

grid_boot <- sim_grid[
  sim_grid$J %in% c(3,5) &
  sim_grid$error_dist == "normal" &
  sim_grid$balance == "balanced", ]
grid_boot$condition_id <- seq_len(nrow(grid_boot))

# Checkpoint kontrolü
if (file.exists(CHECKPOINT)) {
  prev     <- readRDS(CHECKPOINT)
  done_ids <- unique(prev$condition_id)
  cat("Kaldigi yerden devam:", length(done_ids), "kosul tamamlanmis
")
} else {
  prev     <- NULL
  done_ids <- integer(0)
  cat("Yeni Bootstrap simünu basliyor
")
}

grid_todo <- grid_boot[!grid_boot$condition_id %in% done_ids, ]
cat("Kalan kosul:", nrow(grid_todo), "
")

batch_size <- 6
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
      library(lme4); library(pbkrtest)

      cond      <- batch_grid[i,]
      alpha     <- 0.05
      n_rep     <- N_REP
      nsim_boot <- NSIM_BOOT
      n_val     <- as.integer(cond$n_sim)
      results   <- vector("list", n_rep)

      for (r in seq_len(n_rep)) {
        J<-cond$J; icc<-cond$icc; b1<-cond$beta1
        sigma2_u<-icc; sigma2_e<-1-icc
        nj  <- rep(n_val, J)
        N   <- sum(nj); cid <- rep(seq_len(J), times=nj)
        uj  <- rnorm(J, 0, sqrt(sigma2_u))
        eij <- rnorm(N) * sqrt(sigma2_e)
        Xij <- rnorm(N)
        dat <- data.frame(Y=b1*Xij+uj[cid]+eij, X=Xij, cluster=factor(cid))

        fit_null <- tryCatch(
          lme4::lmer(Y~1+(1|cluster), data=dat, REML=FALSE),
          error=function(e) NULL)
        fit_full <- tryCatch(
          lme4::lmer(Y~X+(1|cluster), data=dat, REML=FALSE),
          error=function(e) NULL)
        fit_reml <- tryCatch(
          lme4::lmer(Y~X+(1|cluster), data=dat, REML=TRUE),
          error=function(e) NULL)

        if (!is.null(fit_null) && !is.null(fit_full) && !is.null(fit_reml)) {
          pb <- tryCatch(
            suppressWarnings(
              pbkrtest::PBmodcomp(fit_full, fit_null,
                                  nsim=nsim_boot, seed=NULL)),
            error=function(e) NULL)

          if (!is.null(pb)) {
            pv  <- pb$test["PBtest","p.value"]
            b1e <- lme4::fixef(fit_reml)["X"]
            se1 <- sqrt(vcov(fit_reml)["X","X"])
            ci  <- b1e + c(-1,1)*qnorm(1-alpha/2)*se1
            results[[r]] <- data.frame(
              method="Bootstrap", estimate=b1e, se=se1, pval=pv,
              ci_lo=ci[1], ci_hi=ci[2],
              rejected=as.integer(pv<alpha),
              covered=as.integer(ci[1]<=b1 & b1<=ci[2]),
              rep_id=r, condition_id=cond$condition_id,
              J=cond$J, n=cond$n, icc=cond$icc,
              balance=cond$balance, beta1=cond$beta1,
              error_dist=cond$error_dist)
          }
        }
      }
      do.call(rbind, Filter(Negate(is.null), results))
    },
    .options=furrr_options(seed=TRUE))

  # Checkpoint güe
  if (file.exists(CHECKPOINT)) {
    prev_ckpt <- readRDS(CHECKPOINT)
    updated   <- rbind(prev_ckpt, batch_res)
  } else {
    updated <- batch_res
  }
  saveRDS(updated, CHECKPOINT)
  rm(batch_res, updated); gc()

  elapsed   <- round((proc.time()-start_time)["elapsed"]/60,1)
  remaining <- round((elapsed/b)*(length(batches)-b),1)
  cat(sprintf("  Tamamlandi. Gecen: %.1f dk | Kalan: %.1f dk
",
              elapsed, remaining))
}

final_boot <- readRDS(CHECKPOINT)
saveRDS(final_boot, OUTPUT)
cat("BOOTSTRAP TAMAMLANDI! Toplam satir:", nrow(final_boot), "
")
cat("Sure:", round((proc.time()-start_time)["elapsed"]/3600,2), "saat
")

