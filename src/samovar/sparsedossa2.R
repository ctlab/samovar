#!/usr/bin/env Rscript
# SparseDOSSA2 driver for SamovaR table generation and CV scoring.
# Not an R package; only calls exported SparseDOSSA2 / future functions.

args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

parse_args <- function(argv) {
  out <- list(
    command = "simulate",
    input = NULL,
    output = NULL,
    template = "Stool",
    mode = "fit",
    n_sample = 10L,
    n_feature = NULL,
    new_features = NA,
    seed = 42L,
    workers = 1L,
    verbose = FALSE,
    cv_folds = 5L,
    lambdas = NULL,
    fit_lambda = NULL,
    maxit = NULL,
    max_eval = NULL,
    prec_bits = NULL
  )
  if (length(argv) >= 1 && !startsWith(argv[[1]], "-")) {
    out$command <- argv[[1]]
    argv <- argv[-1]
  }
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    take <- function() {
      if (i + 1 > length(argv)) stop("Missing value for ", key)
      argv[[i + 1]]
    }
    if (key %in% c("--input", "-i")) {
      out$input <- take(); i <- i + 2
    } else if (key %in% c("--output", "-o")) {
      out$output <- take(); i <- i + 2
    } else if (key == "--template") {
      out$template <- take(); i <- i + 2
    } else if (key == "--mode") {
      out$mode <- take(); i <- i + 2
    } else if (key %in% c("--n-sample", "--n_sample", "--N")) {
      out$n_sample <- as.integer(take()); i <- i + 2
    } else if (key %in% c("--n-feature", "--n_feature", "--n-features")) {
      out$n_feature <- as.integer(take()); i <- i + 2
    } else if (key %in% c("--new-features", "--new_features")) {
      nxt <- if (i + 1 <= length(argv) && !startsWith(argv[[i + 1]], "-")) {
        val <- take(); i <- i + 2; val
      } else {
        i <- i + 1; "TRUE"
      }
      out$new_features <- toupper(as.character(nxt)) %in% c("1", "TRUE", "YES", "T")
    } else if (key == "--seed") {
      out$seed <- as.integer(take()); i <- i + 2
    } else if (key %in% c("--workers", "--parallel", "--cores")) {
      out$workers <- as.integer(take()); i <- i + 2
    } else if (key %in% c("--verbose", "-v")) {
      if (i + 1 <= length(argv) && !startsWith(argv[[i + 1]], "-")) {
        out$verbose <- toupper(take()) %in% c("1", "TRUE", "YES", "T"); i <- i + 2
      } else {
        out$verbose <- TRUE; i <- i + 1
      }
    } else if (key %in% c("--cv-folds", "--K", "-K")) {
      out$cv_folds <- as.integer(take()); i <- i + 2
    } else if (key %in% c("--lambdas", "--lambda-grid")) {
      out$lambdas <- take(); i <- i + 2
    } else if (key %in% c("--lambda", "--fit-lambda")) {
      out$fit_lambda <- as.numeric(take()); i <- i + 2
    } else if (key %in% c("--maxit", "--max-it")) {
      out$maxit <- as.integer(take()); i <- i + 2
    } else if (key %in% c("--max-eval", "--max_eval")) {
      out$max_eval <- as.integer(take()); i <- i + 2
    } else if (key %in% c("--prec-bits", "--precBits", "--prec_bits")) {
      out$prec_bits <- as.integer(take()); i <- i + 2
    } else if (key == "--fit") {
      out$mode <- "fit"; i <- i + 1
    } else {
      i <- i + 1
    }
  }
  out
}

limit_blas_threads <- function() {
  for (var in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")) {
    if (!nzchar(Sys.getenv(var, unset = ""))) {
      arg <- list("1")
      names(arg) <- var
      do.call(Sys.setenv, arg)
    }
  }
}

make_control <- function(opts) {
  control <- list(verbose = isTRUE(opts$verbose))
  if (!is.null(opts$maxit) && is.finite(opts$maxit) && opts$maxit >= 1L) {
    control$maxit <- as.integer(opts$maxit)
  }
  numint <- list()
  if (!is.null(opts$max_eval) && is.finite(opts$max_eval) && opts$max_eval >= 1L) {
    numint$max_eval <- as.integer(opts$max_eval)
  }
  if (!is.null(opts$prec_bits) && is.finite(opts$prec_bits) && opts$prec_bits >= 1L) {
    numint$precBits <- as.integer(opts$prec_bits)
  }
  if (length(numint)) {
    control$control_numint <- numint
  }
  control
}

sd2_cluster <- NULL

stop_sd2_cluster <- function() {
  if (!is.null(sd2_cluster)) {
    try(parallel::stopCluster(sd2_cluster), silent = TRUE)
    sd2_cluster <<- NULL
  }
}

# EM's slope test is `if (slope1 <= slopes[2] & ...)` and aborts when a
# spline slope is NA. Treat that region as unfinished so the search continues.
patch_slope_guard <- function() {
  if (!requireNamespace("SparseDOSSA2", quietly = TRUE)) return(invisible(NULL))
  if (!requireNamespace("Rmpfr", quietly = TRUE)) return(invisible(NULL))
  ns <- asNamespace("SparseDOSSA2")
  if (isTRUE(getOption("samovar.slope_guard"))) {
    return(invisible(NULL))
  }
  orig <- get("estimate_one_error", envir = ns)
  guarded <- function(x3, y3, slopes, precBits) {
    nums <- suppressWarnings(as.numeric(slopes))
    if (length(nums) < 3L || any(!is.finite(nums))) {
      return(Rmpfr::mpfr(Inf, precBits = precBits))
    }
    orig(x3, y3, slopes, precBits)
  }
  env <- new.env(parent = environment(orig))
  env$orig <- orig
  environment(guarded) <- env
  assignInNamespace("estimate_one_error", guarded, ns = "SparseDOSSA2")
  options(samovar.slope_guard = TRUE)
  invisible(NULL)
}

# SparseDOSSA2 calls future_vapply without future.seed. The default lives on
# future_lapply; a declared seed stops the "UNRELIABLE VALUE" warning.
patch_future_seed <- function(seed) {
  if (!requireNamespace("future.apply", quietly = TRUE)) return(invisible(NULL))
  ns <- asNamespace("future.apply")
  if (isTRUE(getOption("samovar.future_seed_patched"))) {
    return(invisible(NULL))
  }
  orig <- get("future_lapply", envir = ns)
  formals(orig)$future.seed <- seed
  assignInNamespace("future_lapply", orig, ns = "future.apply")
  options(samovar.future_seed_patched = TRUE)
  invisible(NULL)
}

require_fit_size <- function(mat) {
  n_feature <- sum(rowSums(mat) > 0)
  n_sample <- sum(colSums(mat) > 0)
  if (n_feature < 2L || n_sample < 2L) {
    stop("too few features/samples for SparseDOSSA2 fit", call. = FALSE)
  }
}

fit_or_stop <- function(expr) {
  expr <- substitute(expr)
  tryCatch(
    eval(expr, parent.frame()),
    error = function(e) {
      msg <- conditionMessage(e)
      stop_sd2_cluster()
      if (grepl("missing value where TRUE/FALSE needed", msg, fixed = TRUE) ||
          grepl("leading minor", msg, fixed = TRUE)) {
        stop("too few features/samples for SparseDOSSA2 fit", call. = FALSE)
      }
      stop(e)
    }
  )
}

setup_parallel <- function(workers, cv = FALSE, seed = TRUE) {
  patch_future_seed(seed)
  patch_slope_guard()
  workers <- as.integer(workers %||% 1L)
  if (!requireNamespace("future", quietly = TRUE)) {
    return(invisible(NULL))
  }
  # Wiki uses future::plan(sequential(), ..., multisession()). On future 1.70
  # sequential() must not be called; nested CV uses a strategy list + tweak().
  if (is.na(workers) || workers <= 1L) {
    future::plan(future::sequential)
    return(invisible(NULL))
  }
  stop_sd2_cluster()
  # PSOCK workers do not see a namespace patch made in this process.
  sd2_cluster <<- parallel::makeCluster(workers)
  parallel::clusterExport(sd2_cluster, "patch_slope_guard", envir = .GlobalEnv)
  parallel::clusterEvalQ(sd2_cluster, patch_slope_guard())
  worker_plan <- future::tweak(future::cluster, workers = sd2_cluster)
  if (cv) {
    future::plan(list(future::sequential, future::sequential, worker_plan))
  } else {
    future::plan(worker_plan)
  }
  invisible(NULL)
}

read_feature_matrix <- function(path) {
  # Python always writes features (taxid) as rows and samples as columns.
  # Do not transpose when n_taxa < n_samples — that is a valid abundance table.
  raw <- utils::read.csv(path, row.names = 1, check.names = FALSE)
  mat <- as.matrix(raw)
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  mat
}

write_json <- function(path, obj) {
  txt <- paste0("{", paste(vapply(names(obj), function(k) {
    val <- obj[[k]]
    if (is.null(val) || (length(val) == 1 && is.na(val))) {
      sprintf('"%s": null', k)
    } else if (is.logical(val) && length(val) == 1) {
      sprintf('"%s": %s', k, if (isTRUE(val)) "true" else "false")
    } else if (is.numeric(val) && length(val) == 1) {
      sprintf('"%s": %s', k, format(val, digits = 16, scientific = TRUE))
    } else {
      sprintf('"%s": %s', k, json_quote(as.character(val)[[1]]))
    }
  }, character(1)), collapse = ", "), "}\n")
  writeLines(txt, path, sep = "")
}

json_quote <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  paste0("\"", x, "\"")
}

extract_cv_score <- function(fitted) {
  gof <- fitted$cv_goodness_of_fit
  if (!is.null(gof)) {
    if (is.numeric(gof) && length(gof) >= 1) {
      return(mean(as.numeric(gof), na.rm = TRUE))
    }
    if (is.list(gof)) {
      for (key in c("mean", "score", "value", "cv_goodness_of_fit")) {
        if (!is.null(gof[[key]]) && is.numeric(gof[[key]])) {
          return(mean(as.numeric(gof[[key]]), na.rm = TRUE))
        }
      }
    }
  }
  ll <- fitted$EM_fit$logLik_CV
  if (!is.null(ll)) {
    return(mean(as.numeric(ll), na.rm = TRUE))
  }
  NA_real_
}

opts <- parse_args(args)
limit_blas_threads()
if (is.null(opts$output)) stop("--output is required")
if (!requireNamespace("SparseDOSSA2", quietly = TRUE)) {
  stop("SparseDOSSA2 is not installed. Run ./install.sh SparseDOSSA2")
}
suppressPackageStartupMessages(library(SparseDOSSA2))

if (opts$command %in% c("fitcv", "cv", "score")) {
  if (is.null(opts$input)) stop("--input is required for fitcv")
  mat <- read_feature_matrix(opts$input)
  require_fit_size(mat)
  n_samp <- ncol(mat)
  k <- max(2L, min(as.integer(opts$cv_folds %||% 5L), n_samp))
  setup_parallel(opts$workers, cv = TRUE, seed = opts$seed)
  if (!is.null(opts$seed) && !is.na(opts$seed)) set.seed(opts$seed)
  control <- make_control(opts)
  lambda_vec <- NULL
  if (!is.null(opts$lambdas) && nzchar(opts$lambdas)) {
    lambda_vec <- as.numeric(strsplit(opts$lambdas, ",", fixed = TRUE)[[1]])
    lambda_vec <- lambda_vec[is.finite(lambda_vec)]
  }
  fit_args <- list(data = mat, K = k, control = control)
  if (length(lambda_vec)) fit_args$lambdas <- lambda_vec
  fitted <- fit_or_stop(do.call(fitCV_SparseDOSSA2, fit_args))
  score <- extract_cv_score(fitted)
  write_json(opts$output, list(
    cv_goodness_of_fit = score,
    K = k,
    n_feature = nrow(mat),
    n_sample = ncol(mat)
  ))
  stop_sd2_cluster()
  quit(save = "no", status = 0)
}

new_features <- opts$new_features
if (identical(opts$mode, "fit") || identical(tolower(opts$template), "fit")) {
  if (is.null(opts$input)) stop("--input is required for fit mode")
  mat <- read_feature_matrix(opts$input)
  require_fit_size(mat)
  setup_parallel(opts$workers, cv = FALSE, seed = opts$seed)
  if (!is.null(opts$seed) && !is.na(opts$seed)) set.seed(opts$seed)
  if (is.na(new_features)) new_features <- FALSE
  n_feature <- opts$n_feature %||% nrow(mat)
  control <- make_control(opts)
  lambda <- 1
  if (!is.null(opts$fit_lambda) && is.finite(opts$fit_lambda)) {
    lambda <- opts$fit_lambda
  }
  fitted <- fit_or_stop(fit_SparseDOSSA2(data = mat, lambda = lambda, control = control))
  sim <- SparseDOSSA2(
    template = fitted,
    n_sample = as.integer(opts$n_sample),
    n_feature = as.integer(n_feature),
    new_features = isTRUE(new_features),
    verbose = isTRUE(opts$verbose)
  )
} else {
  setup_parallel(opts$workers, cv = FALSE, seed = opts$seed)
  if (!is.null(opts$seed) && !is.na(opts$seed)) set.seed(opts$seed)
  if (is.na(new_features)) new_features <- TRUE
  n_feature <- opts$n_feature %||% 100L
  sim <- SparseDOSSA2(
    template = as.character(opts$template),
    n_sample = as.integer(opts$n_sample),
    n_feature = as.integer(n_feature),
    new_features = isTRUE(new_features),
    verbose = isTRUE(opts$verbose)
  )
}

counts <- sim$simulated_data
utils::write.csv(counts, opts$output, row.names = TRUE)
stop_sd2_cluster()
quit(save = "no", status = 0)
