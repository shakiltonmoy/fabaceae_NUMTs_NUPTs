setwd("D:/NUMT/R Analysis/Trifolium pratense")

# Load necessary libraries
library(readr)
library(tidyverse)
library(ggbreak)
library(scales)

# Step 1: Read the data
df <- read_csv("nupts_assembled_columns_Trifolium.csv")

# Step 2: Combine simulated values (columns 3 to 1002)
simulated_values <- unlist(df[, 3:1002])
simulated_values <- as.numeric(simulated_values)

# Step 3: Calculate cutoffs
lower_cutoff <- quantile(simulated_values, 0.025, na.rm = TRUE)
upper_cutoff <- quantile(simulated_values, 0.975, na.rm = TRUE)

# Step 4: Density using 99.9% trimmed values
filtered_values <- simulated_values[simulated_values < quantile(simulated_values, 0.999, na.rm = TRUE)]
density_data <- density(filtered_values, na.rm = TRUE, n = 4096, adjust = 3)
density_df <- data.frame(x = density_data$x, y = density_data$y)

# Step 5: Mark inside/outside regions
density_df <- density_df %>%
  mutate(region = ifelse(x >= lower_cutoff & x <= upper_cutoff, "inside", "outside"))

# Add exact cutoff points to density_df
cutoff_points <- data.frame(
  x = c(lower_cutoff, upper_cutoff),
  y = approx(density_df$x, density_df$y, xout = c(lower_cutoff, upper_cutoff))$y,
  region = "inside"
)
density_df <- rbind(density_df, cutoff_points) %>% arrange(x)

# Step 6: Tail extension
max_coverage <- max(simulated_values, na.rm = TRUE)
tail_start <- max(density_df$x)
tail_y <- density_df$y[which.max(density_df$x)]
tail_df <- data.frame(
  x = c(tail_start, max_coverage),
  y = c(tail_y, 0)
)

# Step 7: Plot
ggplot(density_df, aes(x = x, y = y)) +
  geom_line(color = "black", size = 0.7) +
  geom_line(data = tail_df, aes(x = x, y = y), inherit.aes = FALSE, color = "black", size = 0.7) +
  geom_area(data = filter(density_df, region == "inside"), fill = "#00B050", alpha = 0.8) +
  geom_vline(xintercept = c(lower_cutoff, upper_cutoff), color = "red", linetype = "dashed") +
  annotate("text", x = lower_cutoff, y = 0.02, label = paste0("2.5% = ", round(lower_cutoff, 2)),
           color = "black", size = 3, hjust = 1) +
  annotate("text", x = upper_cutoff + 1, y = 0.02, label = paste0("97.5% = ", round(upper_cutoff, 2)),
           color = "black", size = 3, hjust = 0) +
  geom_point(aes(x = max_coverage, y = 0), color = "red", shape = 17, size = 3) +
  annotate("text", x = max_coverage, y = 0.01, label = paste0("Max = ", round(max_coverage, 0)),
           hjust = 1.1, color = "red", size = 3) +
  labs(title = expression(italic("Trifolium pratense")~"(NUPTs)"), x = "Coverage", y = "Density") +
  theme_classic(base_size = 14, base_family = "Times New Roman") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "italic", size = 16, margin = margin(b = 15)),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(color = "black"),
    axis.line.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.text.x.top = element_blank()
  ) +
  scale_x_break(c(120, max_coverage * 0.9), scales = 0.5) +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0.01)),
    labels = label_number(scale_cut = cut_short_scale())
  )

# Step 8: Flag start/stop regions
df$start_flag <- ifelse(df$nupts_start_region < lower_cutoff | df$nupts_start_region > upper_cutoff, "-", "+")
df$stop_flag  <- ifelse(df$nupts_stop_region < lower_cutoff | df$nupts_stop_region > upper_cutoff, "-", "+")

# Step 9: Mark authenticity
df$check_status <- ifelse(df$start_flag == "+" & df$stop_flag == "+", "Authentic", "Possible Artifact")

# Step 10: Count and print summaries
authentic_count <- sum(df$check_status == "Authentic")
artifact_count  <- sum(df$check_status == "Possible Artifact")
print(paste("Authentic count:", authentic_count))
print(paste("Possible Artifact count:", artifact_count))

# Step 11: Add summary to row 2
df$authentic_count <- ""
df$artifact_count <- ""
df$authentic_count[1] <- authentic_count
df$artifact_count[1] <- artifact_count

# Step 12: Outlier quantiles
below_2.5 <- df$nupts_start_region < lower_cutoff | df$nupts_stop_region < lower_cutoff
above_97.5 <- df$nupts_start_region > upper_cutoff | df$nupts_stop_region > upper_cutoff

df$below_2.5_percentile <- ""
df$above_97.5_percentile <- ""
df$below_2.5_percentile[1] <- sum(below_2.5, na.rm = TRUE)
df$above_97.5_percentile[1] <- sum(above_97.5, na.rm = TRUE)

# Step 13: Save results
write_csv(df, "Trifolium_nupts_outlier_results.csv")
ggsave("Trifolium_nupts_outlier.png", width = 6, height = 4, dpi = 300)
