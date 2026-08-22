#' Collapse taxIDs to an NCBI rank via a cached ``taxid|genera_taxid`` table.
#'
#' Original annotation tables are not modified. Lookups go through Python ete3
#' and are stored in ``cache_file`` (and the user cache) for reuse.
#'
#' @param data annotation table with ``taxID_*`` / ``taxid_*`` and optional ``true``.
#' @param rank NCBI rank (default ``genus``). ``none`` / ``exact`` keeps original IDs.
#' @param python optional python executable. Falls back to ``PYTHON_PATH`` then ``python3``.
#' @param cache_file path to a pipe-separated ``taxid|genera_taxid`` table.
#' @return the same table with taxIDs remapped in memory, or the original table if remapping fails.
#' @export
remap_annotation_rank <- function(data, rank = "genus", python = NULL, cache_file = NULL) {
  if (is.null(data) || is.null(rank) || !nzchar(as.character(rank)[1])) {
    return(data)
  }
  rank <- as.character(rank)[1]
  if (tolower(rank) %in% c("none", "exact", "taxid", "raw", "off", "false")) {
    return(data)
  }

  py <- python
  if (is.null(py) || !nzchar(py)) py <- Sys.getenv("PYTHON_PATH", unset = "")
  if (!nzchar(py)) py <- Sys.which("python3")
  if (!nzchar(py)) py <- Sys.which("python")
  if (!nzchar(py)) {
    warning("Python not found; using exact taxIDs instead of rank=", rank)
    return(data)
  }

  tmp_in <- tempfile(fileext = ".csv")
  tmp_out <- tempfile(fileext = ".csv")
  tmp_py <- tempfile(fileext = ".py")
  on.exit(unlink(c(tmp_in, tmp_out, tmp_py)), add = TRUE)
  utils::write.csv(data, tmp_in, row.names = FALSE)
  cache_py <- if (is.null(cache_file) || !nzchar(cache_file)) {
    "None"
  } else {
    shQuote(cache_file, type = "sh")
  }
  writeLines(
    c(
      "import pandas as pd",
      "from samovar.parse_annotators import remap_taxid_dataframe",
      sprintf("df = pd.read_csv(%s)", shQuote(tmp_in, type = "sh")),
      sprintf(
        "remap_taxid_dataframe(df, rank=%s, cache_path=%s).to_csv(%s, index=False)",
        shQuote(rank, type = "sh"),
        cache_py,
        shQuote(tmp_out, type = "sh")
      )
    ),
    tmp_py
  )
  status <- suppressWarnings(system2(py, tmp_py, stdout = TRUE, stderr = TRUE))
  ok <- is.null(attr(status, "status")) || identical(attr(status, "status"), 0L)
  if (!ok || !file.exists(tmp_out) || !file.info(tmp_out)$size) {
    warning(
      sprintf(
        "Could not remap taxIDs to rank=%s via ete3 (%s); using exact taxIDs. %s",
        rank,
        py,
        paste(status, collapse = "\n")
      )
    )
    return(data)
  }
  utils::read.csv(tmp_out, check.names = FALSE, stringsAsFactors = FALSE)
}

is_special_taxon <- function(x) {
  x <- as.character(x)
  tolower(x) %in% c("0", "other", "unclassified", "none", "na", "nan") |
    x %in% c("0")
}

normalize_taxon_token <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "0"
  x[tolower(x) == "other"] <- "other"
  x
}

lookup_taxid_names <- function(taxids, python = NULL, cache_file = NULL) {
  taxids <- unique(normalize_taxon_token(taxids))
  taxids <- taxids[!is_special_taxon(taxids)]
  if (length(taxids) == 0) {
    return(character())
  }

  py <- python
  if (is.null(py) || !nzchar(py)) py <- Sys.getenv("PYTHON_PATH", unset = "")
  if (!nzchar(py)) py <- Sys.which("python3")
  if (!nzchar(py)) py <- Sys.which("python")
  if (!nzchar(py)) {
    warning("Python not found; axis labels will stay as taxIDs")
    return(stats::setNames(taxids, taxids))
  }

  tmp_in <- tempfile(fileext = ".txt")
  tmp_out <- tempfile(fileext = ".tsv")
  tmp_py <- tempfile(fileext = ".py")
  on.exit(unlink(c(tmp_in, tmp_out, tmp_py)), add = TRUE)
  writeLines(taxids, tmp_in)
  cache_py <- if (is.null(cache_file) || !nzchar(cache_file)) {
    "None"
  } else {
    shQuote(cache_file, type = "sh")
  }
  writeLines(
    c(
      "from pathlib import Path",
      "from samovar.parse_annotators import ensure_taxid_name_map",
      sprintf("taxids = Path(%s).read_text().splitlines()", shQuote(tmp_in, type = "sh")),
      sprintf("mapping = ensure_taxid_name_map(taxids, cache_path=%s)", cache_py),
      sprintf("out = Path(%s)", shQuote(tmp_out, type = "sh")),
      "out.write_text('taxid|name\\n' + '\\n'.join(",
      "    f'{k}|{v}' for k, v in sorted(mapping.items()) if k not in {'0', ''}",
      ") + '\\n')"
    ),
    tmp_py
  )
  status <- suppressWarnings(system2(py, tmp_py, stdout = TRUE, stderr = TRUE))
  ok <- is.null(attr(status, "status")) || identical(attr(status, "status"), 0L)
  if (!ok || !file.exists(tmp_out)) {
    warning("Could not resolve taxon names via ete3; using taxIDs. ",
            paste(status, collapse = "\n"))
    return(stats::setNames(taxids, taxids))
  }
  tab <- utils::read.table(
    tmp_out, sep = "|", header = TRUE, stringsAsFactors = FALSE,
    quote = "", comment.char = "", colClasses = "character"
  )
  if (!nrow(tab)) {
    return(stats::setNames(taxids, taxids))
  }
  stats::setNames(tab$name, tab$taxid)
}

fpc_taxon_order <- function(freq_df, col_x, col_y, freq_col = "Freq") {
  x <- normalize_taxon_token(freq_df[[col_x]])
  y <- normalize_taxon_token(freq_df[[col_y]])
  f <- as.numeric(freq_df[[freq_col]])
  f[is.na(f)] <- 0
  core <- unique(c(x, y))
  core <- core[!is_special_taxon(core)]
  if (length(core) <= 1) {
    return(core)
  }
  mat <- matrix(0, length(core), length(core), dimnames = list(core, core))
  for (i in seq_along(f)) {
    a <- y[i]
    b <- x[i]
    if (a %in% core && b %in% core) {
      mat[a, b] <- mat[a, b] + f[i]
    }
  }
  mat_s <- mat + t(mat)
  diag(mat_s) <- diag(mat)
  ord <- tryCatch({
    if (nrow(mat_s) < 2 || stats::sd(as.numeric(mat_s)) == 0) {
      order(rowSums(mat_s), decreasing = TRUE)
    } else {
      pc <- stats::prcomp(mat_s, center = TRUE, scale. = FALSE)
      order(pc$x[, 1])
    }
  }, error = function(e) {
    warning("FPC ordering failed; using abundance. ", conditionMessage(e))
    order(rowSums(mat_s), decreasing = TRUE)
  })
  core[ord]
}

axis_levels_from_fpc <- function(present, fpc_core) {
  present <- unique(normalize_taxon_token(present))
  core <- fpc_core[fpc_core %in% present]
  extra <- setdiff(present[!is_special_taxon(present)], core)
  spec <- c("0", "other")
  spec <- spec[spec %in% present]
  c(core, extra, spec)
}

make_taxon_labels <- function(breaks, name_map, italic = TRUE) {
  texts <- vapply(as.character(breaks), function(b) {
    if (is_special_taxon(b)) {
      lab <- if (tolower(b) == "other") "other" else "0"
      return(sprintf("'%s'", lab))
    }
    key <- sub("\\.0$", "", b)
    nm <- unname(name_map[key])
    if (length(nm) == 0 || is.na(nm) || !nzchar(nm)) {
      nm <- unname(name_map[b])
    }
    if (length(nm) == 0 || is.na(nm) || !nzchar(nm)) {
      nm <- b
    }
    nm <- as.character(nm)[1]
    nm <- gsub("\\\\", "\\\\\\\\", nm)
    nm <- gsub("'", "\\\\'", nm)
    if (isTRUE(italic)) {
      sprintf("italic('%s')", nm)
    } else {
      sprintf("'%s'", nm)
    }
  }, character(1), USE.NAMES = FALSE)
  parse(text = texts)
}

confusion_heatmap <- function(plot_df, x_col, y_col, x_levels, y_levels,
                              name_map, italic, palette, labels_10,
                              label_tiles = TRUE) {
  lab_fun <- function(breaks) make_taxon_labels(breaks, name_map, italic)
  fill_val <- log10(plot_df$Freq + 1)
  rng <- range(fill_val, na.rm = TRUE)
  span <- diff(rng)
  fill_frac <- if (is.finite(span) && span > 0) {
    (fill_val - rng[1]) / span
  } else {
    rep(0, length(fill_val))
  }
  unclassified <- (as.character(plot_df[[x_col]]) == "0") |
    (as.character(plot_df[[y_col]]) == "0")
  # Last ~3 stops of the 11-color green ramp are very dark; use white type.
  n_pal <- max(2, length(palette))
  dark_cut <- (n_pal - 3) / (n_pal - 1)
  plot_df$label_color <- ifelse(
    fill_frac >= dark_cut,
    "white",
    ifelse(unclassified, "brown", "black")
  )
  gg <- ggplot(plot_df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_tile(aes(fill = log10(Freq + 1))) +
    geom_tile(aes(color = rect), fill = "transparent", show.legend = FALSE, linewidth = 1) +
    scale_color_manual(values = c(`TRUE` = "gray", `FALSE` = "transparent")) +
    scale_fill_gradientn(NULL, colours = palette, labels = labels_10) +
    scale_x_discrete(limits = x_levels, labels = lab_fun, drop = FALSE) +
    scale_y_discrete(limits = rev(y_levels), labels = lab_fun, drop = FALSE)
  if (label_tiles) {
    gg <- gg +
      ggnewscale::new_scale_color() +
      geom_text(
        aes(label = N, color = label_color),
        show.legend = FALSE,
        size = 2.2
      ) +
      scale_color_identity()
  }
  gg
}

#' Visualize annotation results
#'
#' @param data Processed abundance table.
##' Row names: sequence IDs,
##' Column names:
##' - annotators: (starting with `taxID_`);
##' - `true`: for true annotation
##' - `length`: length of sequence
##' - sample
#'
#' @param type character vector.
##'
##' - if present column true_annotation: could be one of
##'   - "f1",
##'   - "R2",
##'   - "confidence",
##'   - "cross-validation",
##'   - or their combination (e.g. c("f1", "R2", "cv", "conf"))
##'
#' @param show_top integer. Number of top annotations to show.
#' @param output_dir character. Directory to save the plots. If NULL, plots are not saved.
#' @param plot logical. If TRUE, plots are printed.
#' @param split logical. If TRUE, plots are split into separate files.
#' @param rank NCBI rank used to regroup taxIDs before F1 / R2 / CV (default genus).
#'   Use ``none`` to compare exact taxIDs. Lookup uses ete3 via Python and is
#'   cached as ``taxid|genera_taxid`` next to the plots; source tables are unchanged.
#' @param use_names logical. If TRUE (default), axis labels are scientific names.
#'   Genus and species names are drawn in italic; ``0`` and ``other`` are not.
#' @param reord character. Shared row/column order for F1 and CV heatmaps.
#'   ``fpc`` (default) uses the first principal component of the count matrix;
#'   ``0`` and ``other`` are excluded from that sort and appended last.
#' @return list of ggplot objects
#' @import ggplot2
#'
#' @example R/examples/check_samovar.R
#' @export

viz_annotation <- function(
    data,
    type = c("f1", "R2", "cv", "conf"),
    show_top = 10,
    output_dir = NULL,
    plot = T,
    split = T,
    rank = "genus",
    use_names = TRUE,
    reord = "fpc"
) {

  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  palette_F1 <- c(
    "#FFFFFF",
    "#FFFFE5",
    "#F7FCB9",
    "#E6E487",
    "#D9F0A3",
    "#ADDD8E",
    "#78C679",
    "#41AB5D",
    "#238443",
    "#006837",
    "#004529"
  )
  palette_taxids <- colorRampPalette(c("lightgreen","#E6E487","pink","lightblue","purple"))

  labels_10 <- function(x) {
    ifelse(x == 0, "0",
           ifelse(x == 1, "10",
                  paste0("10^", x)))
  }

  results <- list()

  cache_file <- NULL
  if (!is.null(output_dir) && !is.null(rank) &&
      !tolower(as.character(rank)[1]) %in% c("none", "exact", "taxid", "raw", "off", "false", "")) {
    rank_key <- tolower(as.character(rank)[1])
    cache_name <- if (rank_key %in% c("genus", "genera", "g")) {
      "taxid_genera_map.tsv"
    } else {
      sprintf("taxid_%s_map.tsv", rank_key)
    }
    cache_file <- file.path(output_dir, cache_name)
  }
  data <- remap_annotation_rank(data, rank = rank, cache_file = cache_file)
  rank_used <- if (is.null(rank) || tolower(as.character(rank)[1]) %in%
                   c("none", "exact", "taxid", "raw", "off", "false", "")) {
    NULL
  } else {
    as.character(rank)[1]
  }
  taxid_xlab <- if (isTRUE(use_names)) {
    if (is.null(rank_used)) "True taxon" else sprintf("True taxon (%s)", rank_used)
  } else if (is.null(rank_used)) {
    "True taxID"
  } else {
    sprintf("True taxID (%s)", rank_used)
  }
  taxid_ylab <- if (isTRUE(use_names)) {
    if (is.null(rank_used)) "Predicted taxon" else sprintf("Predicted taxon (%s)", rank_used)
  } else if (is.null(rank_used)) {
    "Predicted taxID"
  } else {
    sprintf("Predicted taxID (%s)", rank_used)
  }
  rank_key <- if (is.null(rank_used)) "" else tolower(as.character(rank_used)[1])
  italic_taxa <- isTRUE(use_names) && rank_key %in% c("genus", "genera", "g", "species", "sp")

  selected_columns <- (str_detect(colnames(data) , "^taxID_|^N_") &
                         str_detect(colnames(data) , "confidence", negate = T)) %>%
    which

  colnames(data) <- colnames(data) %>%
    str_remove("^taxID_") %>%
    str_remove("^N_")

  if(!("length" %in% colnames(data))) data$length <- NULL
  if(!("sample" %in% colnames(data))) data$sample <- NULL

  name_map <- character()
  if (isTRUE(use_names)) {
    name_cols <- colnames(data)[selected_columns]
    if ("true" %in% colnames(data)) name_cols <- c(name_cols, "true")
    all_ids <- unique(unlist(lapply(name_cols, function(cn) as.character(data[[cn]]))))
    name_cache <- NULL
    if (!is.null(output_dir)) {
      name_cache <- file.path(output_dir, "taxid_names.tsv")
    }
    name_map <- lookup_taxid_names(all_ids, cache_file = name_cache)
  }

  if ("true" %in% colnames(data)){
    for (t in type) {
      if (t %in% c("f1","F1")) {

        gglist <- list()

        for (col in selected_columns){
          tmp_name <- colnames(data)[col]

          tmp <- data[,c(tmp_name, "true")] %>%
            mutate_all(as.character)
          tmp[[tmp_name]] <- normalize_taxon_token(tmp[[tmp_name]])
          tmp$true <- normalize_taxon_token(tmp$true)

          tmp[!(unlist(tmp[,tmp_name]) %in%
                  c(unique(as.character(tmp$true)), "0")),
              tmp_name] <- "other"

          cnt <- tmp %>%
            dplyr::count(.data[[tmp_name]], true, name = "Freq")
          fpc_core <- if (identical(reord, "fpc")) {
            fpc_taxon_order(cnt, tmp_name, "true")
          } else {
            sort(unique(c(cnt[[tmp_name]], cnt$true)[
              !is_special_taxon(c(cnt[[tmp_name]], cnt$true))
            ]))
          }
          true_levels <- axis_levels_from_fpc(tmp$true, fpc_core)
          pred_levels <- axis_levels_from_fpc(tmp[[tmp_name]], fpc_core)
          tmp$true <- factor(as.character(tmp$true), levels = true_levels)
          tmp[[tmp_name]] <- factor(as.character(tmp[[tmp_name]]), levels = pred_levels)

          tmp_table <- tmp %>%
            table() %>%
            as.data.frame()

          tmp_table$rect <- as.character(tmp_table$true) == as.character(tmp_table[[tmp_name]])

          # Calculate F1 score
          true_positives <- sum(tmp_table$Freq[tmp_table$rect])
          false_positives <- sum(tmp_table$Freq[!tmp_table$rect & tmp_table$Freq > 0])
          false_negatives <- sum(tmp_table$Freq[!tmp_table$rect & tmp_table$Freq > 0])

          precision <- true_positives / (true_positives + false_positives)
          recall <- true_positives / (true_positives + false_negatives)
          f1_score <- 2 * (precision * recall) / (precision + recall)

          if(is.na(f1_score)) f1_score<- 0

          caption_text <- if (is.null(rank_used)) {
            sprintf("F1-score: %.3f", f1_score)
          } else {
            sprintf("F1-score (%s): %.3f", rank_used, f1_score)
          }

          n_true <- length(true_levels)
          n_pred <- length(pred_levels)

          plot_df <- tmp_table %>%
            mutate(
              N = ifelse(Freq > 0, as.character(Freq), ""),
              pred_is_unclassified = as.character(.data[[tmp_name]]) == "0"
            )

          n_axis <- max(length(true_levels), length(pred_levels))
          label_tiles <- n_axis <= 16
          gg <- confusion_heatmap(
            plot_df, "true", tmp_name, true_levels, pred_levels,
            name_map, italic_taxa, palette_F1, labels_10, label_tiles
          ) +
            theme_minimal() +
            theme(
              panel.grid = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
              axis.text.y = element_text(size = 7)
            ) +
            coord_equal() +
            xlab(taxid_xlab) + ylab(taxid_ylab) +
            labs(caption = caption_text) +
            ggtitle(tmp_name)

          gglist[[tmp_name]] <- gg
        }

        n_taxa_max <- max(10, length(true_levels) + 2)
        panel_w <- max(5, 0.4 * n_taxa_max)
        panel_h <- panel_w

        # Prefer per-annotator F1 panels when many taxa (ggarrange+cowplot can choke)
        results[["F1"]] <- gglist
        if (!is.null(output_dir)) {
          for (i in seq_along(gglist)) {
            w <- max(8, 0.42 * n_taxa_max)
            ggsave(
              gglist[[i]],
              filename = paste0(output_dir, "/F1_", names(gglist)[i], ".png"),
              width = w,
              height = w,
              limitsize = FALSE
            )
          }
          # Also try a stacked overview; ignore failures
          if (require(ggpubr, quietly = TRUE) && !split) {
            tryCatch({
              stacked <- ggpubr::ggarrange(plotlist = gglist, ncol = 1)
              results[["F1"]] <- stacked
              ggsave(
                stacked,
                filename = paste0(output_dir, "/F1.png"),
                width = max(8, 0.42 * n_taxa_max),
                height = max(8, 0.42 * n_taxa_max) * length(gglist),
                limitsize = FALSE
              )
            }, error = function(e) {
              warning(sprintf("Could not write stacked F1.png: %s", conditionMessage(e)))
            })
          }
        }
        if (plot) {
          for (gg in gglist) print(gg)
        }
      } else if (t %in% c("r2","R2")) {
        gglist <- list()

        for (col in selected_columns){
          tmp_name <- colnames(data)[col]

          tmp <- data[,c(tmp_name, "true")] %>%
            mutate_all(as.character)

          tmp[!(unlist(tmp[,tmp_name]) %in%
                  c(unique(tmp$true),"0")),
              tmp_name] <- "other"

          tmp_table <- tmp %>%
            table() %>%
            as.data.frame()

          tmp_table$rect <- as.character(tmp_table$true) == as.character(tmp_table[[tmp_name]])

          # Calculate R2 score
          true_totals  <- tmp_table %>%
            group_by(true) %>%
            summarise(`true taxID` = sum(Freq), .groups = "drop")

          predicted_totals  <- tmp_table %>%
            group_by(!!sym(tmp_name)) %>%
            summarise(`predicted taxID` = sum(Freq), .groups = "drop")

          colnames(predicted_totals)[1] <- "true"

          R2_table <- left_join(true_totals, predicted_totals, by = "true")

          if(sum(R2_table$`predicted taxID`, na.rm = T) > 0) {
            ss_total <- sum((R2_table$`true taxID` - mean(R2_table$`true taxID`))^2)
            ss_residual <- sum((R2_table$`predicted taxID` - R2_table$`true taxID`)^2)
            r2_score <- 1 - (ss_residual / ss_total)

            if(is.na(r2_score)) r2_score <- 0

            caption_text <- sprintf("R²-score: %.3f", r2_score)

            # Visualize
            gg <- R2_table %>%
              ggplot(aes(y = `true taxID`, `predicted taxID`)) +
              geom_smooth(alpha = .1, method = "lm", formula = 'y ~ x') +
              geom_abline(intercept = 1, color = "red", linetype = 2) +
              geom_text(aes(y = `true taxID` + max(`true taxID`)/10,
                            label = true)) +
              geom_point(aes(color = log10(`true taxID`/`predicted taxID`))) +
              scale_color_gradient2(
                NULL,
                high = "pink",
                low = "yellow",
                mid = "lightgreen",
                midpoint = 0,
              ) +
              theme_minimal() +
              theme(panel.grid.minor = element_blank()) +
              coord_equal() +
              xlab(taxid_xlab) + ylab(taxid_ylab) +
              labs(caption = caption_text) +
              ggtitle(tmp_name)

            gglist[[tmp_name]] <- gg
          } else {
            message(sprintf("No efficient annotaton for: %s", tmp_name))
          }
        }

        if(require(ggpubr,quietly = T) & !split) {
          results[["R2"]] <- ggpubr::ggarrange(plotlist = gglist, ncol = 1)
          if(plot) {
            print(results[["R2"]])
          }
          if(!is.null(output_dir)) {
            ggsave(results[["R2"]], filename = paste0(output_dir, "/R2.png"), width = 5, height = 5*length(gglist))
          }
        } else {
          results[["R2"]] <- gglist
          if(!is.null(output_dir)) {
            for(i in 1:length(gglist)) {
              ggsave(gglist[[i]], filename = paste0(output_dir, "/R2_", names(gglist)[i], ".png"), width = 10, height = 10)
              if(plot) {
                print(gglist[[i]])
              }
            }
          }
        }
      }
    }
  }
  if (("cross-validation" %in% type)|("cv" %in% type)|("CV" %in% type)) {
    gglist <- list()
    for(col1 in selected_columns){
      for(col2 in selected_columns){
        if(col1 > col2){
          tmp_name1 <- colnames(data)[col1]
          tmp_name2 <- colnames(data)[col2]

          raw <- data[, c(tmp_name1, tmp_name2)] %>%
            mutate(across(all_of(c(tmp_name1, tmp_name2)), ~normalize_taxon_token(as.character(.x))))

          if(show_top != 0) {
            s <- raw %>%
              pivot_longer(cols = everything()) %>%
              dplyr::count(value, name = "Freq") %>%
              arrange(-Freq)
            s <- s$value[!is_special_taxon(s$value)]
            s <- utils::head(s, max(1, show_top - 1))

            raw[!(raw[, tmp_name1] %in% c(s, "0")), tmp_name1] <- "other"
            raw[!(raw[, tmp_name2] %in% c(s, "0")), tmp_name2] <- "other"
          }

          cnt <- raw %>%
            dplyr::count(.data[[tmp_name1]], .data[[tmp_name2]], name = "Freq")
          fpc_core <- if (identical(reord, "fpc")) {
            fpc_taxon_order(cnt, tmp_name2, tmp_name1)
          } else {
            sort(unique(c(cnt[[tmp_name1]], cnt[[tmp_name2]])[
              !is_special_taxon(c(cnt[[tmp_name1]], cnt[[tmp_name2]]))
            ]))
          }
          lev1 <- axis_levels_from_fpc(raw[[tmp_name1]], fpc_core)
          lev2 <- axis_levels_from_fpc(raw[[tmp_name2]], fpc_core)
          raw[[tmp_name1]] <- factor(raw[[tmp_name1]], levels = lev1)
          raw[[tmp_name2]] <- factor(raw[[tmp_name2]], levels = lev2)
          tmp <- as.data.frame(table(raw))

          tmp <- tmp %>%
            mutate(N = ifelse(Freq > 0, as.character(Freq), ""))

          tmp$rect <- as.character(tmp[,tmp_name1]) == as.character(tmp[,tmp_name2])

          n_cv <- length(unique(c(as.character(lev1), as.character(lev2))))
          label_tiles <- n_cv <= 16

          gg <- confusion_heatmap(
            tmp, tmp_name2, tmp_name1, lev2, lev1,
            name_map, italic_taxa, palette_F1, labels_10, label_tiles
          ) +
            theme_minimal() +
            theme(
              panel.grid = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
              axis.text.y = element_text(size = 7)
            ) +
            coord_equal() +
            xlab(tmp_name2) + ylab(tmp_name1) +
            ggtitle(if (is.null(rank_used)) {
              paste(tmp_name1, "vs", tmp_name2)
            } else {
              paste0(tmp_name1, " vs ", tmp_name2, " (", rank_used, ")")
            })

          gglist[[paste(tmp_name1, "vs", tmp_name2)]] <- gg
        }
      }
    }

    results[["CV"]] <- gglist
    if (!is.null(output_dir)) {
      for (i in seq_along(gglist)) {
        safe <- gsub("[^A-Za-z0-9._-]+", "_", names(gglist)[i])
        ggsave(
          gglist[[i]],
          filename = paste0(output_dir, "/CV_", safe, ".png"),
          width = 12,
          height = 12,
          limitsize = FALSE
        )
      }
      if (require(ggpubr, quietly = TRUE) && !split) {
        tryCatch({
          stacked <- ggpubr::ggarrange(plotlist = gglist, ncol = 1)
          results[["CV"]] <- stacked
          ggsave(
            stacked,
            filename = paste0(output_dir, "/CV.png"),
            width = 8,
            height = 8 * max(1, length(gglist)),
            limitsize = FALSE
          )
        }, error = function(e) {
          warning(sprintf("Could not write stacked CV.png: %s", conditionMessage(e)))
        })
      }
    }
    if (plot) for (gg in gglist) print(gg)
  }

  if ((("confidence" %in% type)|("conf" %in% type))&
      (any(str_detect(colnames(data), "_conf")))) {
    gglist <- list()
    gglist2 <- list()

    selected_columns <- which(str_detect(colnames(data), "_conf"))
    match_columns <- str_remove(colnames(data)[selected_columns], "_conf.*")

    for (i in 1:length(selected_columns)) {
      tmp_name <- match_columns[i]
      tmp_conf <- selected_columns[i]

      tmp <- data.frame(
        x = data[[tmp_conf]],
        y = data[[tmp_name]] %>% as.character()
      )

      if ("true" %in% colnames(data)) {
        tmp$true <- data$true
      }

      if(show_top) {
        tmp_top <- tmp %>%
          count(y) %>%
          arrange(n) %>%
          head(show_top)

        tmp <- tmp %>%
          subset(y %in% tmp_top$y)
      }

      gg <- tmp %>%
        ggplot(aes(x = x, y = y)) +
        geom_boxplot(aes(fill = y), alpha = 0.5,
                     position = position_nudge(y = .25),
                     width = .1, show.legend = F,
                     outliers = F) +
        geom_jitter(aes(color = y), alpha = .1,
                    size = .5, height = .1, width = 0, show.legend = F) +
        theme_minimal() +
        scale_color_manual("taxID", values =
                             palette_taxids(length(unique(tmp$y)) )) +
        scale_fill_manual("taxID", values =
                             palette_taxids(length(unique(tmp$y)) )) +
        theme(panel.grid.minor = element_blank()) +
        xlab("confidence") +
        ylab("taxID") +
        scale_x_continuous(n.breaks = 3) +
        ggtitle(paste(tmp_name, "confidence"))

      gglist[[paste(tmp_name, "confidence")]] <- gg

      if ("true" %in% colnames(data)) {
        gg2 <- gg +
          facet_wrap(~true, scales = "free_y") +
          xlab("true taxID") + ylab("predicted taxID")
        gglist2[[paste(tmp_name, "confidence", "true")]] <- gg2
      }
    }
    if(require(ggpubr,quietly = T) & !split) {
      results[["confidence"]] <- ggpubr::ggarrange(plotlist = gglist, ncol = 1)
      if("true" %in% colnames(data)) {
        results[["confidence_true"]] <- ggpubr::ggarrange(plotlist = gglist2, ncol = 1)
      }
      if(!is.null(output_dir)) {
        ggsave(results[["confidence"]], filename = paste0(output_dir, "/confidence.png"), width = 5*length(gglist[[i]]$data$true), height = 10)
        if("true" %in% colnames(data)) {
          ggsave(results[["confidence_true"]], filename = paste0(output_dir, "/confidence_true.png"), width = 10, height = 10)
        }
      }
      if(plot) {
        print(results[["confidence"]])
        if("true" %in% colnames(data)) {
          print(results[["confidence_true"]])
        }
      } else {
        for(i in 1:length(gglist)) {
          ggsave(gglist[[i]], filename = paste0(output_dir, "/confidence_", names(gglist)[i], ".png"), width = 10, height = 10)
          if("true" %in% colnames(data)) {
            ggsave(gglist2[[i]], filename = paste0(output_dir, "/confidence_true_", names(gglist)[i], ".png"), width = 10, height = 10)
          }
          if(plot) {
            print(gglist[[i]])
            if("true" %in% colnames(data)) {
              print(gglist2[[i]])
            }
          }
        }
      }
    }
  }
  return(results)
}
