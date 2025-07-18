setwd("D:/NUMT/R Analysis/Pisum sativum")

# Load necessary libraries
library(readr)
library(tidyverse)
library(ggbreak)
library(scales)

# Step 1: Read the data
df <- read_csv("numts_assembled_columns_Pisum.csv")

# Step 2: Combine all simulated values (columns 3 to 1002)
simulated_values <- unlist(df[, 3:1002])
simulated_values <- as.numeric(simulated_values)

# Step 3: Calculate 2.5% and 97.5% cutoffs
lower_cutoff <- quantile(simulated_values, 0.025, na.rm = TRUE)
upper_cutoff <- quantile(simulated_values, 0.975, na.rm = TRUE)

# Step 4: Filter values to smooth density
filtered_values <- simulated_values[simulated_values < quantile(simulated_values, 0.999, na.rm = TRUE)]
density_data <- density(filtered_values, na.rm = TRUE, n = 4096, adjust = 3)
density_df <- data.frame(x = density_data$x, y = density_data$y)

# Mark inside vs outside regions
density_df <- density_df %>%
  mutate(region = ifelse(x >= lower_cutoff & x <= upper_cutoff, "inside", "outside"))

# Add exact cutoff points for complete fill
cutoff_points <- data.frame(
  x = c(lower_cutoff, upper_cutoff),
  y = approx(density_df$x, density_df$y, xout = c(lower_cutoff, upper_cutoff))$y,
  region = "inside"
)
density_df <- rbind(density_df, cutoff_points) %>% arrange(x)

# Step 5: Prepare tail line
max_coverage <- max(simulated_values, na.rm = TRUE)
tail_start <- max(density_df$x)
tail_y <- density_df$y[which.max(density_df$x)]
tail_df <- data.frame(
  x = c(tail_start, max_coverage),
  y = c(tail_y, 0)
)

# Step 6: Plot
ggplot(density_df, aes(x = x, y = y)) +
  geom_line(color = "black", size = 0.7) +
  geom_line(data = tail_df, aes(x = x, y = y), inherit.aes = FALSE, color = "black", size = 0.7) +
  geom_area(data = filter(density_df, region == "inside"), fill = "#002060", alpha = .8) +
  geom_vline(xintercept = c(lower_cutoff, upper_cutoff), color = "red", linetype = "dashed") +
  annotate("text", x = lower_cutoff, y = 0.01, label = paste0("2.5% = ", round(lower_cutoff, 2)),
           color = "black", size = 3, hjust = 1.2) +
  annotate("text", x = upper_cutoff + 1, y = 0.01, label = paste0("97.5% = ", round(upper_cutoff, 2)),
           color = "black", size = 3, hjust = 0) +
  geom_point(aes(x = max_coverage, y = 0), color = "red", shape = 17, size = 3) +
  annotate("text", x = max_coverage, y = 0.005, label = paste0("Max = ", round(max_coverage, 0)),
           hjust = 1.1, color = "red", size = 3) +
  labs(title = expression(italic("Pisum sativum")~ "(NUMTs)"), x = "Coverage", y = "Density") +
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
  scale_x_break(c(250, max_coverage * 0.9), scales = 0.5) +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0.01)),
    labels = label_number(scale_cut = cut_short_scale())
  )

# Step 7: Flag start/stop values
df$start_flag <- ifelse(df$numts_start_region < lower_cutoff | df$numts_start_region > upper_cutoff, "-", "+")
df$stop_flag <- ifelse(df$numts_stop_region < lower_cutoff | df$numts_stop_region > upper_cutoff, "-", "+")

# Step 8: Mark as Authentic or Possible Artifact
df$check_status <- ifelse(df$start_flag == "+" & df$stop_flag == "+", "Authentic", "Possible Artifact")

# Step 9: Count classification
authentic_count <- sum(df$check_status == "Authentic")
artifact_count  <- sum(df$check_status == "Possible Artifact")

print(paste("Authentic count:", authentic_count))
print(paste("Possible Artifact count:", artifact_count))

# Step 10: Write counts to first row
df$authentic_count <- ""
df$artifact_count <- ""
df$authentic_count[1] <- authentic_count
df$artifact_count[1] <- artifact_count

# Step 11: Count entries outside cutoffs (either end)
below_2.5 <- df$numts_start_region < lower_cutoff | df$numts_stop_region < lower_cutoff
above_97.5 <- df$numts_start_region > upper_cutoff | df$numts_stop_region > upper_cutoff
below_2.5_count <- sum(below_2.5, na.rm = TRUE)
above_97.5_count <- sum(above_97.5, na.rm = TRUE)

df$below_2.5_percentile <- ""
df$above_97.5_percentile <- ""
df$below_2.5_percentile[1] <- below_2.5_count
df$above_97.5_percentile[1] <- above_97.5_count

# Step 12: Save
write_csv(df, "Pisum_numts_outlier_results.csv")
ggsave("Pisum_numts_outlier.png", width = 6, height = 4, dpi = 300)
