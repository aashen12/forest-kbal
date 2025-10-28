library(tidyverse)

df <- read_csv("results/kernel_comparison.csv")

andy_theme <- function() {
  theme_minimal() + 
    theme(text = element_text(size = 16))
}



kij_plot = df %>% 
  # dplyr::filter(rf > 0) %>% 
  ggplot(aes(x = bart, y = rf)) +
  geom_point(alpha = 0.6, color = "forestgreen", size = 2.5) +
  # geom_smooth(method = "lm", se = FALSE, color = "blue") + 
  labs(x = "BART", y = "RF") + 
  andy_theme()

ggsave(
  filename = "paper-figs/kij_plot.pdf",
  plot     = kij_plot,
  device   = "pdf",      # base grDevices::pdf()
  width    = 6.5,          # double-column width
  height   = 5,        # balanced height
  units    = "in",
  useDingbats = FALSE
)


# Overlapping histograms of BART vs RF (assumes `df` in memory)

df_hist <- df %>%
  dplyr::select(bart, rf) %>%                     # only the two columns
  tidyr::pivot_longer(everything(),
                      names_to = "method",
                      values_to = "value") %>%
  dplyr::filter(is.finite(value))

# common binwidth for aligned histograms
rng <- range(df_hist$value)
binwidth <- diff(rng) / 40   # ~40 bins; adjust as needed

ggplot(df_hist, aes(x = value, fill = method)) +
  geom_histogram(
    position = "identity",
    alpha = 0.45,
    binwidth = binwidth,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = c(bart = "#1B9E77", rf = "#D95F02"),
    name = NULL,
    labels = c("BART", "RF")
  ) +
  labs(x = "Value", y = "Count") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 11)
  )


