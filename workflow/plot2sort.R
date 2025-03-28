options(warn = -1)  # Suppress warnings
suppressMessages({  # Suppress messages
  library(arrow)
  library(tidyverse)
  library(ggsci)
  library(ggalluvial)
  library(ggdist)
})

futurama_unlimited_palette <- function(input_colour_number = 12, output_colour_number) {
  planet_express <- pal_futurama("planetexpress")(input_colour_number)
  output_colour <- colorRampPalette(planet_express)(output_colour_number)
  return (output_colour)
}

print("Generating plots...")

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Please provide exactly three arguments: input folder, output directory, and filtering threshold")
}

input.dir <- args[1]
output.dir <- args[2]
plot_filter_threshold <- as.numeric(args[3])

# Find parquet files
main_files <- list.files(input.dir, pattern = "_main\\.parquet$", full.names = TRUE)
accessory_files <- list.files(input.dir, pattern = "_accessory\\.parquet$", full.names = TRUE)

if (length(main_files) == 0) {
  stop("No main parquet files found in the input folder.")
}
if (length(accessory_files) == 0) {
  stop("No accessory parquet files found in the input folder.")
}

# Read and combine main vs. accessory
all.main <- purrr::map_dfr(main_files, arrow::read_parquet)
all.accessory <- purrr::map_dfr(accessory_files, arrow::read_parquet)

all.full <- bind_rows(all.main, all.accessory)

# Distributions
mean_bit.main <- mean(all.main$mean_bitscore)
q1_bit.main <- quantile(all.main$mean_bitscore, 0.25)
median_bit.main <- quantile(all.main$mean_bitscore, 0.5)
q3_bit.main <- quantile(all.main$mean_bitscore, 0.75)

mean_bit.accessory <- mean(all.accessory$mean_bitscore)
q1_bit.accessory <- quantile(all.accessory$mean_bitscore, 0.25)
median_bit.accessory <- quantile(all.accessory$mean_bitscore, 0.5)
q3_bit.accessory <- quantile(all.accessory$mean_bitscore, 0.75)

all.counted_probe <- all.full %>%
  group_by(species, virus, probe, family, abbreviation) %>%
  summarise(count = dplyr::n(), .groups = "keep") %>%
  ungroup()

main.counted_probe <- all.main %>%
  group_by(species, virus, probe, family, abbreviation) %>%
  summarise(count = dplyr::n(), .groups = "keep") %>%
  ungroup()

accessory.counted_probe <- all.accessory %>%
  group_by(species, virus, probe, family, abbreviation) %>%
  summarise(count = dplyr::n(), .groups = "keep") %>%
  ungroup()

#########################################
# Density plots
#########################################
density_bitscore_plot <- function(data, q1_bit, median_bit, q3_bit) {
  manual_colours = futurama_unlimited_palette(3, length(unique(data$probe)))

  ggplot(data, aes(
    x = mean_bitscore,
    fill = probe,
    colour = probe,
    weight = mean_identity
  )) +
    geom_density(alpha = 0.4, adjust = 3) +
    scale_fill_manual(values = manual_colours) +
    scale_color_manual(values = manual_colours) +
    labs(
      x = "HSP Bitscore",
      y = "Density"
    ) +
    geom_vline(
      xintercept = c(q1_bit, median_bit, q3_bit),
      color = "black",
      linetype = "dashed",
      linewidth = 0.2
    ) +
    annotate("text",
             x = q1_bit + 20,
             y = 1e-4,
             label = "Q1") +
    annotate("text",
             x = median_bit + 20,
             y = 1e-4,
             label = "Q2") +
    annotate("text",
             x = q3_bit + 20,
             y = 1e-4,
             label = "Q3") +
    scale_x_continuous(breaks = seq(0, max(data$mean_bitscore, na.rm = TRUE), by = 100)) +
    theme_minimal() +
    theme(
      text = element_text(face = "bold"),
      axis.title.y = element_text(
        size = 15L,
        face = "bold",
        vjust = 5,
        margin = margin(l = 10)
      ),
      axis.title.x = element_text(
        size = 15L,
        face = "bold",
        vjust = -5,
        margin = margin(b = 15)
      ),
      axis.text.x = element_text(size = 12L, face = "bold"),
      axis.text.y = element_text(size = 12L, face = "bold"),
      margin = margin(t = 10, l = 10, b = 10, r = 10)
    )
}

all.main.density.plot <- density_bitscore_plot(all.main, q1_bit.main, median_bit.main, q3_bit.main)
all.accessory.density.plot <- density_bitscore_plot(all.accessory, q1_bit.accessory, median_bit.accessory, q3_bit.accessory)

#########################################
# Bar plots
#########################################
bar_plot <- function(data) {

  data <- data %>%
    mutate(species = str_replace_all(species, "_", " "))ick %>%
    group_by(species, family) %>%
    summarise(count = sum(count), .groups = 'drop')

  ggplot(data) +
    aes(
      x = species,
      y = count,
      fill = family,
      alpha = count
    ) +
    geom_col(color = "black", linewidth = 0.2) +
    scale_fill_futurama() +
    scale_alpha_continuous(range = c(0.5, 1)) +
    theme_void() +
    labs(fill = "Group",
         x = NULL,
         y = "Count") +
    theme(
      text = element_text(face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold", angle = 90, vjust = 3, margin = margin(l = 10)),
      axis.text.x = element_text(size = 11, face = "bold.italic", angle = 45, vjust = 1, hjust = 1, margin = margin(t = -5)),
      axis.text.y = element_text(size = 11, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 10),
      plot.margin = margin(t = 10, l = 30, b = 60)
    ) +
    guides(alpha = "none")
}

all.full.bar.plot <- bar_plot(all.counted_probe)
all.main.bar.plot <- bar_plot(main.counted_probe)
all.accessory.bar.plot <- bar_plot(accessory.counted_probe)

#########################################
# Raincloud plots
#########################################
raincloud_bitscore_plot <- function(data) {
  required_cols <- c("probe", "mean_bitscore")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Error: Missing columns in data:", paste(missing_cols, collapse = ", ")))
  }
  
  num_colours <- length(unique(data$probe))
  manual_colours <- futurama_unlimited_palette(3, num_colours)
  
  plot <- ggplot(data, aes(y = probe, x = mean_bitscore, fill = probe, color = probe)) +
    ggdist::stat_halfeye(
      adjust = 0.5, 
      justification = 0, 
      .width = c(0.5, 0.95),
      point_colour = NA, 
      alpha = 0.5
    ) +
    geom_boxplot(
      width = 0.15, 
      outlier.shape = NA, 
      alpha = 0.7
    ) +
    ggdist::stat_dots(
      side = "right", 
      dotsize = 0.1,
      alpha = 0.01,
      justification = 0,
      binwidth = 0.2
    ) +
    scale_fill_manual(values = manual_colours) +
    scale_color_manual(values = manual_colours) +
    theme_minimal() +
    labs(
      x = "HSP Bitscore",
      y = "Probe"
    ) +
    theme(
      text = element_text(face = "bold"),
      plot.margin = margin(t = 10, l = 30, b = 10, r = 30),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  return(plot)
}

main.full.raincloud.bitscore_probe.plot <- raincloud_bitscore_plot(all.main)
accessory.full.raincloud.bitscore_probe.plot <- raincloud_bitscore_plot(all.accessory)


##################################################
# Sankey plots
##################################################
sankey_species_probe_plot <- function(data, filter_threshold = plot_filter_threshold) {
  grand_total <- sum(data$count)
  
  species_sums <- data %>%
    group_by(species) %>%
    summarise(species_count = sum(count)) %>%
    arrange(desc(species_count))
  
  probe_sums <- data %>%
    group_by(probe) %>%
    summarise(probe_count = sum(count)) %>%
    arrange(desc(probe_count))
  
  # 3) Apply threshold-based filtering using those totals.
  data_filtered <- data %>%
    inner_join(
      probe_sums %>% 
        filter(probe_count >= filter_threshold * grand_total) %>% 
        select(probe),
      by = "probe"
    ) %>%
    inner_join(
      species_sums %>% 
        filter(species_count >= filter_threshold * grand_total) %>% 
        select(species),
      by = "species"
    )
  
  # 4) Convert species and probe to factors ordered by descending total count
  data_filtered <- data_filtered %>%
    mutate(
      species = factor(species, levels = species_sums$species),
      probe   = factor(probe,   levels = probe_sums$probe)
    )
  
  # Check required columns
  required_cols <- c("species", "probe", "count")
  missing_cols <- setdiff(required_cols, names(data_filtered))
  if (length(missing_cols) > 0) {
    stop(paste("Error: Missing columns in data:", paste(missing_cols, collapse = ", ")))
  }
  
  # Create palette + Sankey
  num_colours <- max(length(unique(data_filtered$species)), length(unique(data_filtered$probe)))
  manual_colours <- futurama_unlimited_palette(12, num_colours)
  
  ggplot(data_filtered, aes(axis1 = species, axis2 = probe, y = count)) +
    geom_alluvium(aes(fill = species, alpha = count)) +
    geom_stratum(alpha = 0, color = "black", linewidth = 0.2) +
    geom_text(
      aes(
        size = after_stat(count),
        label = sprintf("%s [%d]", after_stat(stratum), after_stat(count))
      ),
      nudge_x = -0.15,
      fontface = "bold",
      hjust = 0,
      direction = "y",
      stat = "stratum"
    ) +
    scale_size_continuous(range = c(2, 5)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    theme_void() +
    theme(legend.position = "none") +
    scale_fill_manual(values = manual_colours)
}

sankey_family_probe_plot <- function(data, filter_threshold = plot_filter_threshold) {
  grand_total <- sum(data$count)
  
  family_sums <- data %>%
    group_by(family) %>%
    summarise(family_count = sum(count)) %>%
    arrange(desc(family_count))
  
  probe_sums <- data %>%
    group_by(probe) %>%
    summarise(probe_count = sum(count)) %>%
    arrange(desc(probe_count))
  
  data_filtered <- data %>%
    inner_join(
      probe_sums %>% 
        filter(probe_count >= filter_threshold * grand_total) %>% 
        select(probe),
      by = "probe"
    ) %>%
    inner_join(
      family_sums %>% 
        filter(family_count >= filter_threshold * grand_total) %>% 
        select(family),
      by = "family"
    )
  
  # Factor re-level
  data_filtered <- data_filtered %>%
    mutate(
      family = factor(family, levels = family_sums$family),
      probe  = factor(probe,  levels = probe_sums$probe)
    )
  
  required_cols <- c("family", "probe", "count")
  missing_cols <- setdiff(required_cols, names(data_filtered))
  if (length(missing_cols) > 0) {
    stop(paste("Error: Missing columns in data:", paste(missing_cols, collapse = ", ")))
  }
  
  num_colours <- max(length(unique(data_filtered$family)), length(unique(data_filtered$probe)))
  manual_colours <- futurama_unlimited_palette(12, num_colours)
  
  ggplot(data_filtered, aes(axis1 = family, axis2 = probe, y = count)) +
    geom_alluvium(aes(fill = family, alpha = count)) +
    geom_stratum(alpha = 0, color = "black", linewidth = 0.2) +
    geom_text(
      aes(
        size = after_stat(count),
        label = sprintf("%s [%d]", after_stat(stratum), after_stat(count))
      ),
      nudge_x = -0.15,
      fontface = "bold",
      hjust = 0,
      direction = "y",
      stat = "stratum"
    ) +
    scale_size_continuous(range = c(2, 5)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    theme_void() +
    theme(legend.position = "none") +
    scale_fill_manual(values = manual_colours)
}

sankey_species_family_plot <- function(data, filter_threshold = plot_filter_threshold) {
  grand_total <- sum(data$count)
  
  family_sums <- data %>%
    group_by(family) %>%
    summarise(family_count = sum(count)) %>%
    arrange(desc(family_count))
  
  species_sums <- data %>%
    group_by(species) %>%
    summarise(species_count = sum(count)) %>%
    arrange(desc(species_count))
  
  data_filtered <- data %>%
    inner_join(
      family_sums %>% 
        filter(family_count >= filter_threshold * grand_total) %>% 
        select(family),
      by = "family"
    ) %>%
    inner_join(
      species_sums %>% 
        filter(species_count >= filter_threshold * grand_total) %>% 
        select(species),
      by = "species"
    )
  
  # Factor re-level
  data_filtered <- data_filtered %>%
    mutate(
      species = factor(species, levels = species_sums$species),
      family  = factor(family,  levels = family_sums$family)
    )
  
  required_cols <- c("species", "family", "count")
  missing_cols <- setdiff(required_cols, names(data_filtered))
  if (length(missing_cols) > 0) {
    stop(paste("Error: Missing columns in data:", paste(missing_cols, collapse = ", ")))
  }
  
  num_colours <- max(length(unique(data_filtered$species)), length(unique(data_filtered$family)))
  manual_colours <- futurama_unlimited_palette(12, num_colours)
  
  ggplot(data_filtered, aes(axis1 = species, axis2 = family, y = count)) +
    geom_alluvium(aes(fill = species, alpha = count)) +
    geom_stratum(alpha = 0, color = "black", linewidth = 0.2) +
    geom_text(
      aes(
        size = after_stat(count),
        label = sprintf("%s [%d]", after_stat(stratum), after_stat(count))
      ),
      nudge_x = -0.15,
      fontface = "bold",
      hjust = 0,
      direction = "y",
      stat = "stratum"
    ) +
    scale_size_continuous(range = c(2, 5)) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    theme_void() +
    theme(legend.position = "none") +
    scale_fill_manual(values = manual_colours)
}

#########################################
# Balloon plot
#########################################
balloon_virus_species_plot <- function(data) {
  required_cols <- c("abbreviation", "species", "probe", "count")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Error: Missing columns in data:", paste(missing_cols, collapse = ", ")))
  }
  
  num_colours <- length(unique(data$probe))
  manual_colours <- futurama_unlimited_palette(5, num_colours)
  
  plot <- ggplot(data, aes(x = abbreviation, y = factor(species, levels = rev(unique(species))))) +
    geom_point(aes(
      color = probe,
      fill = probe,
      alpha = count,
      size = count ^ 3
    ), shape = 15) +
    scale_color_manual(values = manual_colours) +
    scale_fill_manual(values = manual_colours) +
    facet_wrap(~ probe) +
    theme_minimal() +
    labs(
      x = "Virus",
      y = "Species"
    ) +
    scale_alpha_continuous(range = c(0.4, 1)) +
    theme(
      text = element_text(face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 9, face = "bold"),
      axis.text.x = element_text(size = 9, face = "bold", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9, face = "bold.italic"),
      strip.text = element_text(size = 12, face = "bold"),
    ) +
    guides(size = "none")
  
  return(plot)
}

main.full.balloon.virus_species.plot <- balloon_virus_species_plot(main.counted_probe)
accessory.full.balloon.virus_species.plot <- balloon_virus_species_plot(accessory.counted_probe)

#########################################
# Generate & save the Sankey plots
#########################################
main.full.sankey.species_probe.plot  <- sankey_species_probe_plot(main.counted_probe)
accessory.full.sankey.species_probe.plot <- sankey_species_probe_plot(accessory.counted_probe)

main.full.sankey.family_probe.plot <- sankey_family_probe_plot(main.counted_probe)
accessory.full.sankey.family_probe.plot <- sankey_family_probe_plot(accessory.counted_probe)

main.full.sankey.species_family.plot <- sankey_species_family_plot(main.counted_probe)
accessory.full.sankey.species_family.plot <- sankey_species_family_plot(accessory.counted_probe)

# Finally, save all plots
ggsave(
  filename = file.path(output.dir, "main_density.png"),
  plot = all.main.density.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "accessory_density.png"),
  plot = all.accessory.density.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "main_raincloud.png"),
  plot = main.full.raincloud.bitscore_probe.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "accessory_raincloud.png"),
  plot = accessory.full.raincloud.bitscore_probe.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "full_bar.png"),
  plot = all.full.bar.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "main_bar.png"),
  plot = all.main.bar.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "accessory_bar.png"),
  plot = all.accessory.bar.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "main_sankey_a.png"),
  plot = main.full.sankey.species_probe.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "accessory_sankey_a.png"),
  plot = accessory.full.sankey.species_probe.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "main_sankey_b.png"),
  plot = main.full.sankey.family_probe.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "accessory_sankey_b.png"),
  plot = accessory.full.sankey.family_probe.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "main_sankey_c.png"),
  plot = main.full.sankey.species_family.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "accessory_sankey_c.png"),
  plot = accessory.full.sankey.species_family.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "main_balloon.png"),
  plot = main.full.balloon.virus_species.plot,
  width = 15,
  height = 12,
  dpi = 300
)
ggsave(
  filename = file.path(output.dir, "accessory_balloon.png"),
  plot = accessory.full.balloon.virus_species.plot,
  width = 15,
  height = 15,
  dpi = 300
)
