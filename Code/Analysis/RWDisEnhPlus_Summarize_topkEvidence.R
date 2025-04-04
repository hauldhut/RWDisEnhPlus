setwd("~/Manuscripts/98RWDisEnhPlus/Results")

# Load required libraries
library(ggplot2)
library(dplyr)

# Step 1: Define the topk values (10 to 100, in increments of 10)
topk_values <- seq(10, 100, by = 10)

# Step 2: Initialize vectors to store the number of evidence associations
evidence_RWDisEnh <- numeric(length(topk_values))
evidence_RWDisEnhPlus <- numeric(length(topk_values))

# Step 3: Loop through topk values to read files and count evidence
for (i in seq_along(topk_values)) {
  topk <- topk_values[i]
  
  # File names for both methods
  file_RWDisEnh <- paste0("RWDisEnh_predict_top", topk, "_EvidenceSum.txt")
  file_RWDisEnhPlus <- paste0("RWDisEnhPlus_predict_top", topk, "_EvidenceSum.txt")
  
  # Read the files and count the number of rows (evidence associations)
  # For RWDisEnh
  if (file.exists(file_RWDisEnh)) {
    data_RWDisEnh <- read.table(
      file_RWDisEnh,
      header = TRUE,
      sep = "\t",
      fill = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    evidence_RWDisEnh[i] <- nrow(data_RWDisEnh)
  } else {
    warning(paste("File not found:", file_RWDisEnh))
    evidence_RWDisEnh[i] <- 0
  }
  
  # For RWDisEnhPlus
  if (file.exists(file_RWDisEnhPlus)) {
    data_RWDisEnhPlus <- read.table(
      file_RWDisEnhPlus,
      header = TRUE,
      sep = "\t",
      fill = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    evidence_RWDisEnhPlus[i] <- nrow(data_RWDisEnhPlus)
  } else {
    warning(paste("File not found:", file_RWDisEnhPlus))
    evidence_RWDisEnhPlus[i] <- 0
  }
}

# Step 4: Create a data frame for plotting
summary_df <- data.frame(
  topk = rep(topk_values, times = 2),
  Evidence = c(evidence_RWDisEnh, evidence_RWDisEnhPlus),
  Method = rep(c("RWDisEnh", "RWDisEnh+"), each = length(topk_values))
)

fontsize = 16
# Step 5: Draw the bar plot using ggplot2
p <- ggplot(summary_df, aes(x = factor(topk), y = Evidence, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(
    x = "Top ranked enhancer (k)",
    y = "Number of evidence associations"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = fontsize),
    axis.text = element_text(size = fontsize),
    legend.position = "top",               # Legend at the top (outside the plot)
    legend.title = element_blank(),        # Remove the legend title ("Method")
    legend.text = element_text(size = fontsize),
    legend.background = element_blank(),   # Remove the legend background
    legend.box.background = element_blank()  # Remove the legend box background
  ) +
  scale_fill_manual(values = c("RWDisEnh" = "#1f77b4", "RWDisEnh+" = "#ff7f0e"))

# Step 6: Display the plot
print(p)

# Step 7: Save the plot (optional)
ggsave("topkEvidence.pdf", plot = p, width = 10, height = 5, dpi = 300)
ggsave("topkEvidence.png", plot = p, width = 10, height = 5, dpi = 300)

