# survalis_colors.R
survalis_palette <- c(
  "primary" = "#3E64FF",      # blue
  "secondary" = "#00C49A",    # turquoise
  "accent" = "#FF6B6B",       # coral red
  "neutral" = "#ADB5BD",      # gray
  "background" = "#F8F9FA",   # light background
  "text" = "#212529"          # dark text
)

theme_survalis <- function(base_size = 12, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = element_text(face = "bold", size = base_size + 2, color = survalis_palette["text"]),
      axis.title = element_text(face = "bold", color = survalis_palette["text"]),
      axis.text = element_text(color = survalis_palette["text"]),
      legend.position = "right",
      panel.grid.major = element_line(color = survalis_palette["neutral"], size = 0.2),
      panel.grid.minor = element_blank()
    )
}



