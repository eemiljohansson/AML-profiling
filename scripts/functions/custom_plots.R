library(ggplot2)
library(dplyr)
library(ggbeeswarm)

# For DE analysis

plot_protein_boxplots <- function(data, metadata, top_proteins, disease_filter = NULL, title = "Top differentially expressed proteins", subtitle = "Acute VTE vs. Acute VTE (controls)") {
  # If disease_filter is NULL, use all unique diseases from metadata
  if (is.null(disease_filter)) {
    disease_filter <- unique(metadata$Disease)
  }
  
  # Join the metadata with the main data
  plot_data <- data %>%
    right_join(metadata %>%
                 filter(Disease %in% disease_filter), by = "DAid") %>%
    filter(Assay %in% top_proteins$Assay)
  
  # Prepare labels from top_proteins
  labels <- top_proteins %>%
    mutate(
      LogFC = round(logFC, 4),
      UnadjustedPValue = formatC(P.Value, format = "e", digits = 2),
      AdjustedPValue = formatC(adj.P.Val, format = "e", digits = 2),
      label = paste("LogFC:", LogFC, "\nUnadjusted p-value:", UnadjustedPValue, "\nAdjusted p-value:", AdjustedPValue)
    ) %>%
    select(Assay, label)
  
  # Create the plot
  plot <- ggplot(plot_data, aes(x = Disease, y = NPX, fill = Disease)) +
    geom_quasirandom(alpha = 0.5, show.legend = FALSE) +
    geom_boxplot(color = "black", outlier.color = NA, alpha = 0.8) +
    facet_wrap(~Assay, ncol = 5, scales = "free_y") +
    geom_text(data = labels, aes(label = label, x = Inf, y = Inf), hjust = 1.1, vjust = 1.1, inherit.aes = FALSE) +  # specify inherit.aes = FALSE
    labs(x = "Disease", y = "NPX", title = title, subtitle = subtitle) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "right")
  
  return(plot)
}

# For ML analysis

plot_protein_boxplots_ml <- function(data, metadata, top_proteins, disease_filter = c("Acute venous thromboembolism", "Acute venous thromboembolism (controls)"), title = "Top differentially expressed proteins", subtitle = "Acute VTE vs. Acute VTE (controls)") {
  # Join the metadata with the main data
  plot_data <- data %>%
    right_join(metadata %>%
                 filter(Disease %in% disease_filter), by = "DAid") %>%
    filter(Assay %in% top_proteins$Assay)
  
  options(repr.plot.width = 15, repr.plot.height = 5) 
  
  # Create the plot
  plot <- ggplot(plot_data, aes(x = Disease, y = NPX, fill = Disease)) +
    geom_quasirandom(alpha = 0.5, show.legend = FALSE) +
    geom_boxplot(color = "black", outlier.color = NA, alpha = 0.8) +
    facet_wrap(~Assay, ncol = 5, scales = "free_y") +
    labs(x = "Disease", y = "NPX", title = title, subtitle = subtitle) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(plot)
}


### Volcano plot ###

library(ggplot2)
library(ggrepel)

plot_volcano <- function(data, logFC_column = "logFC", p_value_column = "adj.P.Val", assay_column = "Assay", top_n = 20, sig_up_color = "#FF7176", sig_down_color = "#92C9DA", not_sig_color = "#D3D3D3") {
  ggplot(data, aes_string(x = logFC_column, y = paste0("-log10(", p_value_column, ")"), label = assay_column, fill = "sig")) +
    geom_point(size = 3, alpha = 0.8, shape = 21, stroke = 0.5, color = "black") +
    geom_text_repel(
      data = data %>% arrange_(p_value_column) %>% head(n = top_n),
      size = 4,
      box.padding = 0.8,
      segment.size = 0.5,
      min.segment.length = 0,
      max.overlaps = Inf,
      arrow = arrow(length = unit(0.010, "npc")),
      nudge_x = 0.15,
      nudge_y = 0.5,
      color = "grey50"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("significant up" = sig_up_color, "significant down" = sig_down_color, "not significant" = not_sig_color)) +
    theme_light() +
    theme(
      axis.text = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      x = "log2 fold change",
      y = "-log10(adjusted p-value)",
      fill = "Significance"
    ) +
    guides(color = FALSE, fill = guide_legend(title = "Significance"))
}

### Confusion matrix ###

create_conf_plot <- function(conf_matrix, true_label_A, true_label_B) {
  class_codes <- c("0", "1")
  class_labels <- c(true_label_A, true_label_B)
  
  if (length(dim(conf_matrix)) != 2L ||
      !identical(dim(conf_matrix), c(2L, 2L)) ||
      !setequal(rownames(conf_matrix), class_codes) ||
      !setequal(colnames(conf_matrix), class_codes)) {
    stop("conf_matrix must be a 2 x 2 table with rows and columns named '0' and '1'.")
  }
  
  axis_names <- names(dimnames(conf_matrix))
  if (!is.null(axis_names) &&
      !identical(tolower(axis_names), c("prediction", "truth"))) {
    stop("Expected Prediction in the rows and Truth in the columns.")
  }
  
  # Read the counts using the matrix's class codes.
  conf_df <- expand.grid(
    Prediction = class_codes,
    Truth = class_codes,
    stringsAsFactors = FALSE
  )
  conf_df$n <- as.numeric(conf_matrix[cbind(
    match(conf_df$Prediction, rownames(conf_matrix)),
    match(conf_df$Truth, colnames(conf_matrix))
  )])
  conf_df$Correct <- conf_df$Prediction == conf_df$Truth
  
  conf_df$Prediction <- factor(
    conf_df$Prediction, levels = class_codes, labels = class_labels
  )
  conf_df$Truth <- factor(
    conf_df$Truth, levels = class_codes, labels = class_labels
  )
  
  # Match the order on both axes when read left to right and top to bottom.
  display_order <- c(true_label_B, true_label_A)
  
  ggplot2::ggplot(
    conf_df,
    ggplot2::aes(x = Prediction, y = Truth, fill = Correct)
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = n), color = "black", size = 6) +
    ggplot2::scale_x_discrete(limits = display_order) +
    ggplot2::scale_y_discrete(limits = rev(display_order)) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#C5E0B3", "FALSE" = "#FB8875"),
      guide = "none"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Predicted Condition",
      y = "True Condition",
      title = paste("Confusion Matrix for", true_label_B, "vs.", true_label_A)
    )
}
