library(ggplot2)
library(dplyr)
library(ggbeeswarm)  # for geom_quasirandom

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
  # Convert matrix to a tidy data frame
  conf_df <- as.data.frame.matrix(conf_matrix) %>%
    tidyr::pivot_longer(everything(), names_to = "Predicted_class", values_to = "n") %>%
    mutate(Truth = rep(c(true_label_A, true_label_B), each = 2)) %>%
    mutate(
      Annotation_color = case_when(
        Predicted_class == "1" & Truth == true_label_B ~ "#C5E0B3",  # TP
        Predicted_class == "1" & Truth == true_label_A ~ "#FB8875",  # FP
        Predicted_class == "0" & Truth == true_label_B ~ "#FB8875",  # FN
        Predicted_class == "0" & Truth == true_label_A ~ "#C5E0B3"   # TN
      )
    )
  
  # Adjusting the factor levels for clarity and correct plotting
  conf_df$Predicted_class <- factor(conf_df$Predicted_class, levels = c("0", "1"), labels = c(true_label_A, true_label_B))
  conf_df$Truth <- factor(conf_df$Truth, levels = c(true_label_A, true_label_B))
  
  # Plot the confusion matrix
  plot <- conf_df %>%
    ggplot(aes(x = Predicted_class, y = Truth, fill = Annotation_color)) +
    geom_tile(color = "white") +  # Add white borders for clarity
    geom_text(aes(label = n), color = "black", size = 6) +  # Show count
    scale_fill_identity(name = "Outcome", labels = c("TP", "FP", "FN", "TN"), guide = "none") +
    coord_fixed() +  # Ensure square tiles
    theme_minimal() +  # Minimal theme
    labs(
      x = "Predicted Condition",
      y = "True Condition",
      title = paste("Confusion Matrix for", true_label_B, "vs.", true_label_A, "Classification"),
      fill = "Outcome"
    )
  
  return(plot)
}
