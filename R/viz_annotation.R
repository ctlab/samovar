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
##' - if present column true_annotation: could be one of "f1", "R2", or their combination
##'
##' - else: "cross-validation"
##'
#' @param show_top integer. Number of top annotations to show.
#' @param output_dir character. Directory to save the plots. If NULL, plots are not saved.
#' @param plot logical. If TRUE, plots are printed.
#' @return list of ggplot objects
#' @importFrom tibble tibble
#' @importFrom dplyr mutate_all
#' @importFrom dplyr mutate
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @importFrom dplyr %>%
#' @importFrom dplyr sym
#' @importFrom dplyr left_join
#' @importFrom dplyr pull
#' @importFrom dplyr arrange
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr summarise
#' @importFrom dplyr across
#' @importFrom stringr str_detect
#' @importFrom stringr str_remove
#' @import ggplot2
#' @import ggnewscale
#'
#' @example R/examples/check_samovar.R
#' @export

viz_annotation <- function(
    data,
    type = c("f1", "R2", "cv"),
    show_top = 10,
    output_dir = NULL,
    plot = T
) {

  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  palette_F1 <- c("#FFFFFF00", "#E6E487", "#90EE90")
  labels_10 <- function(x) {
    ifelse(x == 0, "0",
           ifelse(x == 1, "10",
                  parse(text = paste0("10^", x))))
  }

  results <- list()

  selected_columns <- colnames(data) %>%
    str_detect("^taxID_|^N_") %>% which

  colnames(data) <- colnames(data) %>%
    str_remove("^taxID_") %>%
    str_remove("^N_")

  if(!("length" %in% colnames(data))) data$length <- NULL
  if(!("sample" %in% colnames(data))) data$sample <- NULL

  if ("true" %in% colnames(data)){
    for (t in type) {
      if (t %in% c("f1","F1")) {

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

          # Calculate F1 score
          true_positives <- sum(tmp_table$Freq[tmp_table$rect])
          false_positives <- sum(tmp_table$Freq[!tmp_table$rect & tmp_table$Freq > 0])
          false_negatives <- sum(tmp_table$Freq[!tmp_table$rect & tmp_table$Freq > 0])

          precision <- true_positives / (true_positives + false_positives)
          recall <- true_positives / (true_positives + false_negatives)
          f1_score <- 2 * (precision * recall) / (precision + recall)

          if(is.na(f1_score)) f1_score<- 0

          caption_text <- sprintf("F1-score: %.3f", f1_score)

          # Visualize
          gg <- tmp_table %>%
            mutate(N = ifelse(Freq > 0, Freq, "")) %>%
            ggplot(aes(y = .data[[tmp_name]], true %>% fct_rev)) +
            geom_tile(aes(fill = log10(Freq+1))) +
            geom_tile(aes(color = rect), fill = "transparent", show.legend = F, linewidth = 1) +
            scale_color_manual(values = c(`TRUE` = "gray", `FALSE` = "transparent")) +
            scale_fill_gradientn(
              NULL,
              colours = palette_F1,
              labels = labels_10
            ) +

            ggnewscale::new_scale_color() +
            geom_text(aes(label = N, color = .data[[tmp_name]] == 0)) +
            scale_color_manual(values = c(`TRUE` = "brown", `FALSE` = "black")) +
            theme_minimal() +
            theme(panel.grid = element_blank()) +
            coord_equal() +
            xlab("True taxID") + ylab("Predicted taxID") +
            labs(caption = caption_text) +
            ggtitle(tmp_name)

          gglist[[tmp_name]] <- gg
        }

        if(require(ggpubr,quietly = T)) {
          results[["F1"]] <- ggpubr::ggarrange(plotlist = gglist)
          if(plot) {
            print(results[["F1"]])
          }
          if(!is.null(output_dir)) {
            ggsave(results[["F1"]], filename = "F1.png", width = 10, height = 10)
          }
        } else {
          results[["F1"]] <- gglist
          if(!is.null(output_dir)) {
            for(i in 1:length(gglist)) {
              ggsave(gglist[[i]], filename = paste0("F1_", names(gglist)[i], ".png"), width = 10, height = 10)
              if(plot) {
                print(gglist[[i]])
              }
            }
          }
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
              xlab("True taxID") + ylab("Predicted taxID") +
              labs(caption = caption_text) +
              ggtitle(tmp_name)

            gglist[[tmp_name]] <- gg
          } else {
            message(sprintf("No efficient annotaton for: %s", tmp_name))
          }
        }

        if(require(ggpubr,quietly = T)) {
          results[["R2"]] <- ggpubr::ggarrange(plotlist = gglist)
          if(plot) {
            print(results[["R2"]])
          }
          if(!is.null(output_dir)) {
            ggsave(results[["R2"]], filename = "R2.png", width = 10, height = 10)
          }
        } else {
          results[["R2"]] <- gglist
          if(!is.null(output_dir)) {
            for(i in 1:length(gglist)) {
              ggsave(gglist[[i]], filename = paste0("R2_", names(gglist)[i], ".png"), width = 10, height = 10)
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

          tmp <- data[,c(tmp_name1, tmp_name2)] %>%
            table() %>%
            as.data.frame() %>%
            mutate(across(all_of(c(tmp_name1, tmp_name2)), as.character))

          if(show_top != 0) {
            s <- tmp %>%
              pivot_longer(cols = -3) %>%
              summarise(Freq = mean(Freq), .by = value) %>%
              arrange(-Freq) %>%
              head(show_top-1) %>%
              pull(value)

            tmp[!(tmp[,tmp_name1] %in% s), tmp_name1] <- "other"
            tmp[!(tmp[,tmp_name2] %in% s), tmp_name2] <- "other"

            tmp <- tmp %>%
              summarise(Freq = sum(Freq), .by = c(all_of(c(tmp_name1, tmp_name2)))) %>%
              mutate(N = ifelse(Freq > 0, Freq, ""))
          }

          tmp$rect <- as.character(tmp[,tmp_name1]) == as.character(tmp[,tmp_name2])

          gg <- tmp %>%
            ggplot(aes(y = .data[[tmp_name1]], .data[[tmp_name2]] %>% fct_rev)) +
            geom_tile(aes(fill = log10(Freq+1))) +
            geom_tile(aes(color = rect), fill = "transparent", show.legend = F, linewidth = 1) +
            geom_text(aes(label = N)) +
            scale_color_manual(values = c(`TRUE` = "gray", `FALSE` = "transparent")) +
            scale_fill_gradientn(
              NULL,
              colours = palette_F1,
              labels = labels_10
            ) +
            ggnewscale::new_scale_color() +
            geom_text(aes(label = N, color =
                            (.data[[tmp_name1]] == 0) |
                            (.data[[tmp_name2]] == 0)
                          )) +
            scale_color_manual(values = c(`TRUE` = "brown", `FALSE` = "black")) +
            theme_minimal() +
            theme(panel.grid = element_blank()) +
            coord_equal() +
            xlab(tmp_name2) + ylab(tmp_name1) +
            ggtitle(paste(tmp_name1, "vs", tmp_name2))

          gglist[[paste(tmp_name1, "vs", tmp_name2)]] <- gg
        }
      }
    }

    if(require(ggpubr,quietly = T)) {
      results[["CV"]] <- ggpubr::ggarrange(plotlist = gglist)
      if(plot) {
        print(results[["CV"]])
      }
      if(!is.null(output_dir)) {
        ggsave(results[["CV"]], filename = "cross-validation.png", width = 10, height = 10)
      }
    } else {
      results[["CV"]] <- gglist
      if(!is.null(output_dir)) {
        for(i in 1:length(gglist)) {
          ggsave(gglist[[i]], filename = paste0("CV_", names(gglist)[i], ".png"), width = 10, height = 10)
          if(plot) {
            print(gglist[[i]])
          }
        }
      }
    }
  }
}
