df <- read_csv("results/kernel_pcs.csv")

andy_theme <- function() {
  theme_minimal() + 
    theme(text = element_text(size = 16))
}


pkbal <- df %>% 
  ggplot(aes(x = KBAL_PC1, y = KBAL_PC2)) +
  geom_point(alpha = 0.8, size = 2, color = "#d62728") + 
  andy_theme() + 
  labs(x = "PC1 (Gauss. Kernel)", y = "PC2 (Gauss. Kernel)")
pkbal

ggsave(
  filename = "paper-figs/kbal_pc.pdf",
  plot     = pkbal,
  device   = "pdf",      # base grDevices::pdf()
  width    = 6,          # double-column width
  height   = 5,        # balanced height
  units    = "in",
  useDingbats = FALSE
)

pbart = data %>% 
  ggplot(aes(x = BART_PC1, y = BART_PC2)) +
  geom_point(alpha = 0.8, size = 2, color = "#ff7f0e") + 
  andy_theme() + 
  labs(x = "PC1 (BART)", y = "PC2 (BART)")
pbart

ggsave(
  filename = "paper-figs/bart_pc.pdf",
  plot     = pbart,
  device   = "pdf",      # base grDevices::pdf()
  width    = 6,          # double-column width
  height   = 5,        # balanced height
  units    = "in",
  useDingbats = FALSE
)


prf = data %>% 
  ggplot(aes(x = RF_PC1, y = RF_PC2)) +
  geom_point(alpha = 0.8, size = 2, color = "#1f77b4") + 
  andy_theme() + 
  labs(x = "PC1 (RF)", y = "PC2 (RF)")
prf

# ggsave(
#   filename = "paper-figs/rf_pc.pdf",
#   plot     = prf,
#   device   = "pdf",      # base grDevices::pdf()
#   width    = 6,          # double-column width
#   height   = 5,        # balanced height
#   units    = "in",
#   useDingbats = FALSE
# )
# 