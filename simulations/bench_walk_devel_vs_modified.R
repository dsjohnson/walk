#!/usr/bin/env Rscript
###############################################################################
## Benchmark: walk 'devel' (GitHub, dsjohnson/walk) vs the modified walk
##            (GitHub, bmcclintock/walk, branch 'speed-up')
##
## 1. Clones both repositories and installs two versions side by side, each in
##    its own library (both are called 'walk', so each is run in a separate R
##    process with callr):
##      devel     : https://github.com/dsjohnson/walk     @ devel
##      modified  : https://github.com/bmcclintock/walk   @ speed-up-slim
## 2. Simulates tracks with langevinSSM::simLangevin (exampleCovs, Argos-type
##    error ellipses, irregular sampling).
## 3. For each grid resolution, in each version:
##      * proc_Langevin per animal (devel: terra::buffer + mvtnorm::pmvnorm in R;
##        modified: neighborhood = "buffer" and "disc", 1 and 2 threads), plus
##        the modified one-call version with id = "id";
##      * make_walk_data + fit_ctmc(fit = FALSE);
##      * -2 log-likelihood at several parameter vectors
##        (devel: ctmc_n2ll per animal, summed, because devel has no joint
##         likelihood; modified: ctmc_n2ll_joint per animal and in one joint
##         call with 1 and 2 threads).
##    The modified version also evaluates its likelihood on devel's own L
##    matrices, which isolates the likelihood code from the proc step.
## 4. Checks: L (devel vs modified buffer) identical sparsity pattern and values;
##    -2ll devel vs modified (buffer) equal to a relative tolerance; joint =
##    sum of per-animal; disc vs buffer reported (not expected to be identical:
##    disc also picks up the edge cells the polygonised buffer misses).
## 5. Optional: full nlminb fit of one animal in both versions.
##
## Usage:   Rscript bench_walk_devel_vs_modified.R [quick]
##          or, in an R session: source("bench_walk_devel_vs_modified.R")
##          (set QUICK <- TRUE before source() for the quick run)
## Env:     WALKBENCH_OUT  output dir for results and write-up (default: working directory)
##          WALKBENCH_DIR  build/scratch dir: clones, libraries (default: <tempdir>/walkbench)
##          WALK_MOD_REPO, WALK_MOD_REF   override the modified repo / branch
##          LANGEVIN_TAR   langevinSSM source tarball; only used if langevinSSM is not
##                         installed (the installed version is used otherwise)
## Output (in the working directory):
##          walkbench_report.md   write-up: setup, checks, speed-up tables
##          walkbench_summary.txt full printed tables and checks
##          walkbench_{proc_times,lik_times,lik_values,L_compare,fit}.csv
##          walkbench_results.rds everything above + config + sessionInfo
###############################################################################

args  <- commandArgs(trailingOnly = TRUE)
if (!exists("QUICK")) QUICK <- "quick" %in% args

CFG <- list(
  work        = Sys.getenv("WALKBENCH_DIR", file.path(tempdir(), "walkbench")),
  repo        = "https://github.com/dsjohnson/walk.git",
  devel_ref   = "origin/devel",
  mod_repo    = Sys.getenv("WALK_MOD_REPO", "https://github.com/bmcclintock/walk.git"),
  mod_ref     = Sys.getenv("WALK_MOD_REF", "speed-up-slim"),
  out         = normalizePath(Sys.getenv("WALKBENCH_OUT", here::here("simulations","results")), mustWork = TRUE),
  langevin_tar = Sys.getenv("LANGEVIN_TAR", ""),
  strip_missing = TRUE,              # drop ctmm / rARPACK from DESCRIPTION/NAMESPACE only if not installed
  # simulation (langevinSSM units: exampleCovs are 400 x 400 cells of size 1)
  seed        = 1,
  par         = list(beta = c(-4, 6, 5, -0.1), sigma = 5, gamma = 0.5),
  nbAnimals   = 4,
  obsPerAnimal = 3000,               # simulation steps (timeStep 0.01) before subsampling
  samplingRate = 10,                 # keep ~1 in 10 -> irregular dt, ~300 fixes per animal
  me          = list(smaj.sd = 3, smin.sd = 1.5, eor.lim = c(0, 180)),
  # grids: aggregation factor of exampleCovs (cell size = factor)
  grids       = c(4, 2, 1),                 
  form_m      = ~ cov1 + cov2 + cov3,
  threads     = c(1, 2),
  nrep        = 3,                    # repeats per timing (median reported)
  rel_tol     = 1e-8,                 # -2ll agreement tolerance (relative)
  L_tol       = 1e-10,                # L agreement tolerance (absolute)
  fit_one     = TRUE                  # full fit of animal 1 on the coarsest grid
)
if (QUICK) { CFG$grids <- 4; CFG$nbAnimals <- 2; CFG$obsPerAnimal <- 1000; CFG$nrep <- 1 }
if (exists("GRIDS")) CFG$grids <- GRIDS   # e.g. GRIDS <- c(4, 2) before source() to skip the 1-cell grid

dir.create(CFG$work, recursive = TRUE, showWarnings = FALSE)
res_dir <- CFG$out
outf <- function(x) file.path(res_dir, paste0("walkbench_", x))
logf <- outf("summary.txt")
t_start <- Sys.time()
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
cat("", file = logf)
run <- function(cmd, args, wd = NULL, ...) {
  old <- if (!is.null(wd)) setwd(wd); on.exit(if (!is.null(old)) setwd(old))
  out <- suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE, ...))
  st <- attr(out, "status"); if (!is.null(st) && st != 0) { cat(out, sep = "\n"); stop(cmd, " failed") }
  invisible(out)
}
for (p in c("callr", "terra", "Matrix")) if (!requireNamespace(p, quietly = TRUE)) stop("Package '", p, "' is required.")

## ---- 1. Build both versions ----------------------------------------------------
repo_dir <- file.path(CFG$work, "walk_repo")
if (!dir.exists(repo_dir)) run("git", c("clone", "-q", CFG$repo, repo_dir))
run("git", c("fetch", "-q", "origin"), wd = repo_dir)
devel_sha <- run("git", c("rev-parse", "--short", CFG$devel_ref), wd = repo_dir)

prep_src <- function(dir) {   # remove compiled objects; strip unavailable optional imports
  unlink(list.files(file.path(dir, "src"), "\\.(o|so|dll)$", full.names = TRUE))
  if (!CFG$strip_missing) return(invisible())
  d <- readLines(file.path(dir, "DESCRIPTION")); n <- readLines(file.path(dir, "NAMESPACE"))
  for (p in c("ctmm", "rARPACK")) if (!requireNamespace(p, quietly = TRUE)) {
    d <- gsub(paste0("\\b", p, ",\\s*"), "", d); d <- gsub(paste0(",\\s*", p, "\\b"), "", d)
    n <- n[!grepl(paste0("importFrom\\(", p, ","), n)]
    message("note: '", p, "' not installed; stripped from ", basename(dir), " DESCRIPTION/NAMESPACE")
  }
  writeLines(d, file.path(dir, "DESCRIPTION")); writeLines(n, file.path(dir, "NAMESPACE"))
}
install_to <- function(src, lib) {
  unlink(lib, recursive = TRUE); dir.create(lib, recursive = TRUE)
  run(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--no-test-load", "--no-docs", "--no-multiarch",
                                       "-l", shQuote(lib), shQuote(src)), env = paste0("R_LIBS=", paste(c(lib, .libPaths()), collapse = .Platform$path.sep)))
}

say("== Building walk versions (work dir: %s)", CFG$work)
src_devel <- file.path(CFG$work, "src_devel"); unlink(src_devel, recursive = TRUE)
run("git", c("worktree", "prune"), wd = repo_dir)
run("git", c("worktree", "add", "-f", "--detach", src_devel, CFG$devel_ref), wd = repo_dir)
prep_src(file.path(src_devel, "walk"))
lib_devel <- file.path(CFG$work, "lib_devel"); install_to(file.path(src_devel, "walk"), lib_devel)
say("devel    : %s @ devel (%s) -> %s", CFG$repo, devel_sha, lib_devel)

# modified: bmcclintock/walk @ speed-up, fetched into the same clone as a second remote
remotes <- run("git", c("remote"), wd = repo_dir)
if (!"modified" %in% remotes) run("git", c("remote", "add", "modified", CFG$mod_repo), wd = repo_dir) else
  run("git", c("remote", "set-url", "modified", CFG$mod_repo), wd = repo_dir)
run("git", c("fetch", "-q", "modified", CFG$mod_ref), wd = repo_dir)
mod_sha <- run("git", c("rev-parse", "--short", "FETCH_HEAD"), wd = repo_dir)
src_mod <- file.path(CFG$work, "src_mod"); unlink(src_mod, recursive = TRUE)
run("git", c("worktree", "prune"), wd = repo_dir)
run("git", c("worktree", "add", "-f", "--detach", src_mod, mod_sha), wd = repo_dir)
src_mod_pkg <- file.path(src_mod, "walk")
mod_desc <- sprintf("%s @ %s (%s)", CFG$mod_repo, CFG$mod_ref, mod_sha)
prep_src(src_mod_pkg)
lib_mod <- file.path(CFG$work, "lib_mod"); install_to(src_mod_pkg, lib_mod)
say("modified : %s -> %s", mod_desc, lib_mod)

## langevinSSM (simulation only)
# Use the installed langevinSSM; install from LANGEVIN_TAR into the work dir only if it is missing.
lib_lsm <- file.path(CFG$work, "lib_langevinSSM")
have_ver <- tryCatch(as.character(packageVersion("langevinSSM", lib.loc = c(lib_lsm, .libPaths()))), error = function(e) "")
if (!nzchar(have_ver)) {
  if (!nzchar(CFG$langevin_tar) || !file.exists(CFG$langevin_tar))
    stop("langevinSSM is not installed; install it or set LANGEVIN_TAR to the source tarball.")
  dir.create(lib_lsm, showWarnings = FALSE)
  say("installing langevinSSM from %s", CFG$langevin_tar)
  run(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--no-docs", "--no-multiarch", "-l", shQuote(lib_lsm), shQuote(CFG$langevin_tar)),
      # one make job: the TMB translation unit needs ~2-3 GB of RAM on its own
      env = c(paste0("R_LIBS=", paste(c(lib_lsm, .libPaths()), collapse = .Platform$path.sep)), "MAKEFLAGS=-j1"))
}
base_libs <- unique(c(lib_lsm, .libPaths()))

## ---- 2. Simulate tracks ----------------------------------------------------------
sim_file <- file.path(CFG$work, "sim.rds")
invisible(callr::r(function(cfg, out) {
  suppressPackageStartupMessages({ library(langevinSSM); library(terra) })
  set.seed(cfg$seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  s <- simLangevin(par = cfg$par, spatialCovs = exampleCovs, nbAnimals = cfg$nbAnimals,
                   obsPerAnimal = cfg$obsPerAnimal, measurementError = cfg$me,
                   subSample = list(samplingRate = cfg$samplingRate))
  s$id <- as.character(s$id)
  covs <- terra::wrap(c(exampleCovs$cov1, exampleCovs$cov2, exampleCovs$cov3))
  saveRDS(list(sim = s, covs = covs, version = as.character(packageVersion("langevinSSM"))), out)
  invisible(NULL)
}, args = list(cfg = CFG, out = sim_file), libpath = base_libs))
S <- readRDS(sim_file)
say("simulated (langevinSSM %s): %d animals, %d fixes (%s), median dt %.3f, median smaj %.2f",
    S$version, length(unique(S$sim$id)), nrow(S$sim), paste(table(S$sim$id), collapse = "/"),
    median(S$sim$dt[S$sim$dt > 0]), median(S$sim$smaj))

## ---- 3. Worker run inside each version ----------------------------------------
worker <- function(version, cfg, sim_file, devel_file) {
  suppressPackageStartupMessages({ library(walk); library(terra); library(Matrix) })
  S <- readRDS(sim_file); s <- S$sim; ids <- unique(s$id)
  # median of n timed runs, after one untimed warm-up run (first-call / allocation effects)
  tmed <- function(expr, n) { e <- substitute(expr); pf <- parent.frame()
  v <- eval(e, pf); tt <- numeric(n); for (k in seq_len(n)) { gc(FALSE); tt[k] <- system.time(v <- eval(e, pf))[["elapsed"]] }
  list(value = v, time = median(tt)) }
  is_mod <- version == "modified"
  dev <- if (is_mod) readRDS(devel_file) else NULL
  proc_t <- list(); lik_t <- list(); lik_v <- list(); Ls <- list(); fit <- NULL
  for (g in cfg$grids) {
    covs <- terra::unwrap(S$covs)
    if (g > 1) covs <- terra::aggregate(covs, g, fun = "mean")
    h <- terra::res(covs)[1]; G <- paste0("grid", g)
    sub_i <- function(i) { d <- s[s$id == i, ]; class(d) <- class(s); d }
    pars <- list(
      p1 = c(log(2^2 / (2 * h^2)), 0, 0, 0),
      p2 = c(log(2^2 / (2 * h^2)), 0.5, -0.5, 0.3),
      p3 = c(log(4^2 / (2 * h^2)), -1, 1, 0.5))
    mp <- ctmc_control(q_r = ctmc_model(~1, link = "log"), q_m = ctmc_model(cfg$form_m, link = "log"), norm = FALSE)
    # proc variants
    variants <- if (is_mod) {
      do.call(c, lapply(c("buffer", "disc"), function(nb) lapply(cfg$threads, function(th) list(nb = nb, th = th))))
    } else list(list(nb = "buffer", th = 1))
    walkdat <- list()
    for (v in variants) {
      key <- paste0(v$nb, "_t", v$th)
      r <- tmed(lapply(ids, function(i) {
        if (is_mod) suppressWarnings(proc_Langevin(sub_i(i), covs, neighborhood = v$nb, nthreads = v$th, time_unit = "secs"))
        else suppressWarnings(proc_Langevin(sub_i(i), covs))
      }), cfg$nrep)
      proc_t[[length(proc_t) + 1]] <- data.frame(version = version, grid = g, cells = terra::ncell(covs),
                                                 neighborhood = v$nb, threads = v$th, call = "per animal", time = r$time)
      if (v$th == 1) {
        Ls[[G]][[v$nb]] <- lapply(r$value, function(z) z$L)
        walkdat[[v$nb]] <- lapply(r$value, function(pd) suppressWarnings(make_walk_data(pd, covs)))
      }
      if (is_mod) {   # one call for all animals (id column), as used for joint fits
        rj <- tmed(suppressWarnings(proc_Langevin(s, covs, id = "id", neighborhood = v$nb, nthreads = v$th,
                                                  time_unit = "secs")), cfg$nrep)
        proc_t[[length(proc_t) + 1]] <- data.frame(version = version, grid = g, cells = terra::ncell(covs),
                                                   neighborhood = v$nb, threads = v$th, call = "one call (id)", time = rj$time)
        if (v$th == 1) walkdat[[paste0(v$nb, "_joint")]] <- suppressWarnings(make_walk_data(rj$value, covs))
      }
    }
    if (is_mod) {   # devel's own L matrices, evaluated with the modified likelihood
      walkdat$devL <- lapply(seq_along(ids), function(k) {
        pd <- list(L = dev$Ls[[G]]$buffer[[k]], times = data.frame(obs = seq_len(nrow(dev$Ls[[G]]$buffer[[k]])),
                                                                   timestamp = sub_i(ids[k])$date, dt = c(0, diff(sub_i(ids[k])$date))))
        attr(pd$times, "time_unit") <- "secs"
        suppressWarnings(make_walk_data(pd, covs)) })
    }
    dls <- lapply(walkdat, function(w) {
      one <- function(wd) suppressWarnings(fit_ctmc(wd, mp, start = list(beta_q_r = pars$p1[1], beta_q_m = pars$p1[-1]),
                                                    fit = FALSE, hessian = FALSE))$data_list
      if (inherits(w, "walk_ddl")) one(w) else lapply(w, one) })
    for (pn in names(pars)) {
      p <- pars[[pn]]
      add <- function(L_from, how, th, r) {
        lik_t[[length(lik_t) + 1]] <<- data.frame(version = version, grid = g, cells = terra::ncell(covs), par = pn,
                                                  L_from = L_from, how = how, threads = th, time = r$time)
        lik_v[[length(lik_v) + 1]] <<- data.frame(version = version, grid = g, par = pn, L_from = L_from, how = how,
                                                  threads = th, n2ll = r$value)
      }
      if (!is_mod) {
        add("devel buffer", "ctmc_n2ll per animal", 1,
            tmed(sum(vapply(dls$buffer, function(dl) as.numeric(ctmc_n2ll(p, dl)), 0)), cfg$nrep))
      } else {
        for (nb in c("buffer", "disc", "devL")) {
          lab <- c(buffer = "modified buffer", disc = "modified disc", devL = "devel buffer")[[nb]]
          add(lab, "ctmc_n2ll_joint per animal", 1,
              tmed(sum(vapply(dls[[nb]], function(dl) as.numeric(ctmc_n2ll_joint(p, dl)), 0)), cfg$nrep))
          if (nb != "devL") for (th in cfg$threads) {
            dl <- dls[[paste0(nb, "_joint")]]; dl$nthreads <- as.integer(th)
            add(lab, "ctmc_n2ll_joint one joint call", th, tmed(as.numeric(ctmc_n2ll_joint(p, dl)), cfg$nrep))
          }
        }
      }
    }
    if (cfg$fit_one && g == max(cfg$grids)) {   # full fit, animal 1, same start, same optimiser
      dl1 <- walkdat$buffer[[1]]
      tf <- system.time(f <- suppressWarnings(fit_ctmc(dl1, mp, start = list(beta_q_r = pars$p1[1], beta_q_m = pars$p1[-1]),
                                                       hessian = FALSE, method = "nlminb")))[["elapsed"]]
      fit <- data.frame(version = version, grid = g, time = tf, n2ll = -2 * f$log_lik,
                        convergence = f$opt$convergence, par = I(list(f$par)))
    }
  }
  list(version = version, walk = as.character(packageVersion("walk")), path = find.package("walk"),
       openmp = if (is_mod) walk:::omp_info()$enabled else NA,
       proc_t = do.call(rbind, proc_t), lik_t = do.call(rbind, lik_t), lik_v = do.call(rbind, lik_v),
       Ls = Ls, fit = fit)
}

say("== Running devel")
dev_file <- file.path(CFG$work, "devel_out.rds")
D <- callr::r(worker, args = list(version = "devel", cfg = CFG, sim_file = sim_file, devel_file = NULL),
              libpath = c(lib_devel, base_libs), show = TRUE)
saveRDS(D, dev_file)
say("== Running modified")
M <- callr::r(worker, args = list(version = "modified", cfg = CFG, sim_file = sim_file, devel_file = dev_file),
              libpath = c(lib_mod, base_libs), show = TRUE)
say("walk loaded from: devel %s (%s) | modified %s (%s, OpenMP %s)", D$path, D$walk, M$path, M$walk, M$openmp)

## ---- 4. Compare -----------------------------------------------------------------
cmpL <- function(A, B) {   # lists of per-animal sparse L (same dims)
  r <- lapply(seq_along(A), function(k) {
    a <- as(A[[k]], "TsparseMatrix"); b <- as(B[[k]], "TsparseMatrix")
    ka <- paste(a@i, a@j); kb <- paste(b@i, b@j); common <- intersect(ka, kb)
    va <- setNames(a@x, ka)[common]; vb <- setNames(b@x, kb)[common]
    dx <- Matrix::drop0(as(A[[k]], "CsparseMatrix") - as(B[[k]], "CsparseMatrix"))@x
    c(nnz_a = length(ka), nnz_b = length(kb), only_a = length(setdiff(ka, kb)), only_b = length(setdiff(kb, ka)),
      max_abs = if (length(common)) max(abs(va - vb)) else NA, max_abs_all = if (length(dx)) max(abs(dx)) else 0)
  })
  r <- do.call(rbind, r)
  c(colSums(r[, 1:4, drop = FALSE]), max_abs = max(r[, "max_abs"]), max_abs_all = max(r[, "max_abs_all"]))
}
Lc <- do.call(rbind, lapply(names(D$Ls), function(G) rbind(
  data.frame(grid = G, compare = "devel vs modified buffer", t(cmpL(D$Ls[[G]]$buffer, M$Ls[[G]]$buffer))),
  data.frame(grid = G, compare = "devel vs modified disc",   t(cmpL(D$Ls[[G]]$buffer, M$Ls[[G]]$disc))))))
Lc$identical_pattern <- Lc$only_a == 0 & Lc$only_b == 0
Lc$pass <- ifelse(Lc$compare == "devel vs modified buffer", Lc$identical_pattern & Lc$max_abs_all <= CFG$L_tol, NA)

V <- rbind(D$lik_v, M$lik_v)
ref <- subset(D$lik_v, select = c(grid, par, n2ll)); names(ref)[3] <- "n2ll_devel"
V <- merge(V, ref, by = c("grid", "par"))
V$abs_diff <- V$n2ll - V$n2ll_devel; V$rel_diff <- abs(V$abs_diff) / abs(V$n2ll_devel)
V$expect_equal <- V$L_from %in% c("devel buffer", "modified buffer")
V$pass <- ifelse(V$expect_equal, V$rel_diff <= CFG$rel_tol, NA)

PT <- rbind(D$proc_t, M$proc_t)
pd <- subset(D$proc_t, select = c(grid, time)); names(pd)[2] <- "devel_time"
PT <- merge(PT, pd, by = "grid"); PT$speedup <- PT$devel_time / PT$time
LT <- rbind(D$lik_t, M$lik_t)
ld <- subset(D$lik_t, select = c(grid, par, time)); names(ld)[3] <- "devel_time"
LT <- merge(LT, ld, by = c("grid", "par")); LT$speedup <- LT$devel_time / LT$time

FT <- rbind(D$fit, M$fit)
write.csv(PT, outf("proc_times.csv"), row.names = FALSE)
write.csv(LT, outf("lik_times.csv"), row.names = FALSE)
write.csv(V,  outf("lik_values.csv"), row.names = FALSE)
write.csv(Lc, outf("L_compare.csv"), row.names = FALSE)
if (!is.null(FT)) write.csv(transform(FT, par = sapply(par, function(x) paste(signif(x, 8), collapse = " "))),
                            outf("fit.csv"), row.names = FALSE)
saveRDS(list(cfg = CFG, devel_sha = devel_sha, mod = mod_desc, mod_sha = mod_sha, proc = PT, lik_time = LT, lik_value = V,
             L = Lc, fit = FT, sessionInfo = sessionInfo()), outf("results.rds"))

## ---- 5. Report -------------------------------------------------------------------
op <- options(width = 200, digits = 10)
pr <- function(title, x) { say("\n== %s", title); o <- capture.output(print(x, row.names = FALSE)); cat(o, sep = "\n"); cat(o, sep = "\n", file = logf, append = TRUE) }
pr("proc_Langevin time (s, median of nrep; all animals)",
   transform(PT[order(PT$grid, PT$version, PT$neighborhood, PT$threads), ], time = round(time, 3), devel_time = round(devel_time, 3), speedup = round(speedup, 1)))
pr("L: devel vs modified", Lc)
pr("-2 log-likelihood (summed over animals)",
   transform(V[order(V$grid, V$par, V$version, V$L_from), c("grid", "par", "version", "L_from", "how", "threads", "n2ll", "abs_diff", "rel_diff", "pass")],
             rel_diff = signif(rel_diff, 3), abs_diff = signif(abs_diff, 3)))
pr("Likelihood time (s, one evaluation, all animals)",
   transform(LT[order(LT$grid, LT$par, LT$version), c("grid", "cells", "par", "version", "L_from", "how", "threads", "time", "speedup")],
             time = round(time, 4), speedup = round(speedup, 1)))
if (!is.null(FT)) pr("Full fit, animal 1, coarsest grid (nlminb, no Hessian)",
                     transform(FT, par = sapply(par, function(x) paste(signif(x, 8), collapse = " "))))
options(op)

## ---- 6. Write-up (Markdown) -------------------------------------------------------
write_report <- function() {
  f <- outf("report.md"); L <- character()
  w <- function(...) L <<- c(L, sprintf(...))
  md_tab <- function(df) {
    df <- as.data.frame(df); df[] <- lapply(df, function(x) if (is.numeric(x)) format(x, trim = TRUE, drop0trailing = TRUE) else as.character(x))
    c(paste0("| ", paste(names(df), collapse = " | "), " |"), paste0("|", strrep("---|", ncol(df))),
      apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
  }
  rng <- function(x, d = 1) if (length(x) == 0) "NA" else if (diff(range(x)) < 10^-d) sprintf(paste0("%.", d, "f"), x[1]) else
    sprintf(paste0("%.", d, "f–%.", d, "f"), min(x), max(x))
  cells <- setNames(unique(PT[, c("grid", "cells")])$cells, unique(PT[, c("grid", "cells")])$grid)
  gs <- sort(unique(PT$grid), decreasing = TRUE)
  
  w("# walk devel vs modified (speed-up-slim): benchmark report")
  w("")
  w("Generated %s by `bench_walk_devel_vs_modified.R`%s. Run time %.1f min.", format(Sys.time(), "%Y-%m-%d %H:%M"),
    if (QUICK) " (quick mode)" else "", as.numeric(difftime(Sys.time(), t_start, units = "mins")))
  w("")
  w("## Setup")
  w("")
  w("- **devel:** %s @ `devel` (%s).", CFG$repo, devel_sha)
  w("- **modified:** %s.", mod_desc)
  w("- **Versions loaded:** walk %s (devel) and %s (modified; OpenMP %s). Each ran in its own R process.", D$walk, M$walk, M$openmp)
  w("- **Simulated data:** langevinSSM %s, `simLangevin` (underdamped, beta = (%s), sigma = %g, gamma = %g) on `exampleCovs`.",
    S$version, paste(CFG$par$beta, collapse = ", "), CFG$par$sigma, CFG$par$gamma)
  w("  - %d animals, %d fixes, subsampled to irregular times (median dt %.3f).", length(unique(S$sim$id)), nrow(S$sim),
    median(S$sim$dt[S$sim$dt > 0]))
  w("  - Argos-type error ellipses (median smaj %.2f).", median(S$sim$smaj))
  w("- **Grids:** `exampleCovs` aggregated by %s, i.e. %s cells.", paste(gs, collapse = ", "),
    paste(format(cells[as.character(gs)], big.mark = ",", trim = TRUE), collapse = ", "))
  w("- **Model:** movement `%s`, log links, `norm = FALSE`. −2ll evaluated at %d parameter vectors.",
    deparse(CFG$form_m), length(unique(V$par)))
  w("- **Timing:** median of %d run(s) after one warm-up. Machine: %s, %d cores, R %s.", CFG$nrep, Sys.info()[["sysname"]],
    parallel::detectCores(), getRversion())
  w("")
  w("## Checks")
  w("")
  w("- **L, devel vs modified `buffer` (abs tol %g): %s.**", CFG$L_tol, if (ok_L) "PASS" else "FAIL")
  for (i in which(Lc$compare == "devel vs modified buffer"))
    w("  - %s: %d vs %d non-zeros, %d/%d cells only in one version, max abs diff %.3g.", Lc$grid[i], Lc$nnz_a[i], Lc$nnz_b[i],
      Lc$only_a[i], Lc$only_b[i], Lc$max_abs_all[i])
  w("- **−2ll, devel vs modified (`buffer` L, and the modified likelihood on devel's own L), rel tol %g: %s.**",
    CFG$rel_tol, if (ok_V) "PASS" else "FAIL")
  w("  - Max relative difference %.3g over %d comparisons.", max(V$rel_diff[V$expect_equal]), sum(V$expect_equal & V$version == "modified"))
  w("- **Joint call vs sum of per-animal values (modified):** max abs diff %.3g.", max(abs(jt$n2ll.x - jt$n2ll.y)))
  if (!is.null(FT) && nrow(FT) == 2) {
    w("- **Full fit of animal 1** (%s cells, nlminb, same start):", format(cells[as.character(FT$grid[1])], big.mark = ","))
    w("  - −2ll: devel %.7f, modified %.7f (difference %.3g).", FT$n2ll[1], FT$n2ll[2], FT$n2ll[2] - FT$n2ll[1])
    w("  - Max abs parameter difference %.3g; convergence codes %s.", max(abs(FT$par[[1]] - FT$par[[2]])), paste(FT$convergence, collapse = "/"))
  }
  dsc <- Lc[Lc$compare == "devel vs modified disc", ]
  dv <- V[V$L_from == "modified disc", ]
  w("- **`disc` vs `buffer` (not expected to be identical):**")
  w("  - Extra cells in `disc` that the polygonised buffer misses: %s.",
    paste(sprintf("%s cells: %d", format(cells[sub("^grid", "", dsc$grid)], big.mark = ",", trim = TRUE), dsc$only_b), collapse = "; "))
  w("  - Its −2ll differs from devel by %.2g to %.2g (relative).", min(dv$rel_diff), max(dv$rel_diff))
  w("")
  w("## Speed")
  w("")
  w("Speed-up = devel time / modified time. All times are for all animals together.")
  w("")
  w("### `proc_Langevin`")
  w("")
  pt <- PT[PT$version == "modified", ]
  tab <- do.call(rbind, lapply(gs, function(g) {
    d <- pt[pt$grid == g, ]; dv <- PT$devel_time[PT$grid == g][1]
    cell <- function(nb, cl, th) { r <- d[d$neighborhood == nb & d$call == cl & d$threads == th, ]
    if (nrow(r)) sprintf("%.3f s (%.0f×)", r$time, r$speedup) else "—" }
    data.frame(cells = format(cells[as.character(g)], big.mark = ","), devel = sprintf("%.2f s", dv),
               `buffer, per animal, 1 thr` = cell("buffer", "per animal", 1), `buffer, one call, 2 thr` = cell("buffer", "one call (id)", max(CFG$threads)),
               `disc, per animal, 1 thr` = cell("disc", "per animal", 1), `disc, one call, 2 thr` = cell("disc", "one call (id)", max(CFG$threads)),
               check.names = FALSE)
  }))
  L <- c(L, md_tab(tab)); w("")
  w("### −2 log-likelihood (one evaluation; range over parameter vectors)")
  w("")
  lt <- LT[LT$version == "modified", ]
  tab <- do.call(rbind, lapply(gs, function(g) {
    d <- lt[lt$grid == g, ]; dv <- unique(LT[LT$grid == g & LT$version == "devel", c("par", "time")])
    cell <- function(Lf, how, th) { r <- d[d$L_from == Lf & d$how == how & d$threads == th, ]
    if (nrow(r)) sprintf("%s s (%s×)", rng(r$time, 3), rng(r$speedup, 0)) else "—" }
    data.frame(cells = format(cells[as.character(g)], big.mark = ","), devel = sprintf("%s s", rng(dv$time, 2)),
               `per animal, devel's L` = cell("devel buffer", "ctmc_n2ll_joint per animal", 1),
               `per animal, buffer` = cell("modified buffer", "ctmc_n2ll_joint per animal", 1),
               `joint, 1 thr` = cell("modified buffer", "ctmc_n2ll_joint one joint call", 1),
               `joint, 2 thr` = cell("modified buffer", "ctmc_n2ll_joint one joint call", max(CFG$threads)),
               check.names = FALSE)
  }))
  L <- c(L, md_tab(tab)); w("")
  if (!is.null(FT) && nrow(FT) == 2)
    w("**Full fit, animal 1:** devel %.1f s → modified %.2f s (%.0f×).", FT$time[1], FT$time[2], FT$time[1] / FT$time[2])
  w("")
  w("## Notes")
  w("")
  w("- **Devel has no joint likelihood,** so the devel figure is `ctmc_n2ll` per animal, summed.")
  w("- **\"per animal, devel's L\"** runs the modified likelihood on exactly the inputs devel used. It isolates the likelihood code from the `proc_Langevin` step.")
  w("- **`disc`** uses every cell that touches the exact error circle, while `buffer` uses `terra::buffer` + `terra::cells(touches = TRUE)`. Only `buffer` is expected to match devel.")
  if (any(!vapply(c("ctmm", "rARPACK"), requireNamespace, TRUE, quietly = TRUE)))
    w("- **`ctmm`/`rARPACK`:** not installed on this machine, so they were removed from DESCRIPTION/NAMESPACE of both builds before installing.")
  w("- **Files:** `walkbench_summary.txt` (all tables), `walkbench_*.csv`, `walkbench_results.rds`.")
  writeLines(L, f, useBytes = TRUE)
  invisible(f)
}

ok_L <- all(Lc$pass, na.rm = TRUE); ok_V <- all(V$pass, na.rm = TRUE)
say("\nCHECK  L identical (devel vs modified buffer, tol %g): %s", CFG$L_tol, if (ok_L) "PASS" else "FAIL")
say("CHECK  -2ll equal (devel vs modified, buffer / devel L, rel tol %g): %s  (max rel diff %.3g)",
    CFG$rel_tol, if (ok_V) "PASS" else "FAIL", max(V$rel_diff[V$expect_equal]))
jt <- merge(subset(V, version == "modified" & how == "ctmc_n2ll_joint per animal", c(grid, par, L_from, n2ll)),
            subset(V, version == "modified" & how == "ctmc_n2ll_joint one joint call", c(grid, par, L_from, threads, n2ll)),
            by = c("grid", "par", "L_from"))
say("CHECK  joint call = sum of per-animal (modified): max abs diff %.3g", max(abs(jt$n2ll.x - jt$n2ll.y)))
write_report()
say("Results and write-up (walkbench_report.md) written to %s", res_dir)
