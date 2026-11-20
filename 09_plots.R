library(ggplot2)
library(patchwork)
library(cowplot)

GSE63952_wang <- readRDS("0_RT_GSE63952_wang.rds")
GSE63952_li <- readRDS("0_RT_GSE63952_li.rds")
GSE63952_zhao <- readRDS("0_RT_GSE63952_zhao.rds")

GSE102971_wang <- readRDS("0_RT_GSE102971_wang.rds")
GSE102971_li <- readRDS("0_RT_GSE102971_li.rds")
GSE102971_zhao <- readRDS("0_RT_GSE102971_zhao.rds")

GSE23515_wang <- readRDS("0_RT_GSE23515_wang.rds")
GSE23515_li <- readRDS("0_RT_GSE23515_li.rds")
GSE23515_zhao <- readRDS("0_RT_GSE23515_zhao.rds")

plot_dose_predictions <- function(averaged_testing_matrix, title) {
  
  ggplot(averaged_testing_matrix, aes(x = factor(Dose), y = CatBoost)) +
    stat_summary(aes(y = Dose, fill = "Real Dose"), 
                 fun = mean, 
                 geom = "bar", 
                 width = 0.6, 
                 alpha = 1,
                 color = "steelblue4",
                 linewidth = 0.8) +
    geom_jitter(aes(color = "Predicted Dose"),
                width = 0.15, 
                height = 0, 
                size = 2.5, 
                alpha = 0.8) +
    scale_fill_manual(name = NULL, 
                      values = c("Real Dose" = "honeydew3")) +
    scale_color_manual(name = NULL, 
                       values = c("Predicted Dose" = "darkorange3")) +
    labs(title = title,
         x = "Dose (Gy)",
         y = "Predicted Dose (Gy)") +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 10),
      legend.position = "top"
    )
}


dataset_name <- "GSE63952"
plot_GSE63952 <- (plot_dose_predictions(get(paste0(dataset_name, "_wang")) , paste0("Dataset ", dataset_name, " with the Wang signature")) +
                    plot_dose_predictions(get(paste0(dataset_name, "_li")) , paste0("Dataset ", dataset_name, " with the Li signature")) + 
                    plot_dose_predictions(get(paste0(dataset_name, "_zhao")) , paste0("Dataset ", dataset_name, " with the Zhao signature"))) +
  plot_layout(ncol = 1, guides = "collect")

dataset_name <- "GSE102971"
plot_GSE102971 <- (plot_dose_predictions(get(paste0(dataset_name, "_wang")) , paste0("Dataset ", dataset_name, " with the Wang signature")) +
                     plot_dose_predictions(get(paste0(dataset_name, "_li")) , paste0("Dataset ", dataset_name, " with the Li signature")) + 
                     plot_dose_predictions(get(paste0(dataset_name, "_zhao")) , paste0("Dataset ", dataset_name, " with the Zhao signature"))) +
  plot_layout(ncol = 1, guides = "collect")

dataset_name <- "GSE23515"
plot_GSE23515 <- (plot_dose_predictions(get(paste0(dataset_name, "_wang")) , paste0("Dataset ", dataset_name, " with the Wang signature")) +
                    plot_dose_predictions(get(paste0(dataset_name, "_li")) , paste0("Dataset ", dataset_name, " with the Li signature")) + 
                    plot_dose_predictions(get(paste0(dataset_name, "_zhao")) , paste0("Dataset ", dataset_name, " with the Zhao signature"))) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "top")


library(patchwork)

plots <- (
  plot_dose_predictions(get("GSE63952_wang"), "GSE63952 – Wang") +
    plot_dose_predictions(get("GSE63952_li"), "GSE63952 – Li") +
    plot_dose_predictions(get("GSE63952_zhao"), "GSE63952 – Zhao") +
    plot_dose_predictions(get("GSE102971_wang"), "GSE102971 – Wang") +
    plot_dose_predictions(get("GSE102971_li"), "GSE102971 – Li") +
    plot_dose_predictions(get("GSE102971_zhao"), "GSE102971 – Zhao") +
    plot_dose_predictions(get("GSE23515_wang"), "GSE23515 – Wang") +
    plot_dose_predictions(get("GSE23515_li"), "GSE23515 – Li") +
    plot_dose_predictions(get("GSE23515_zhao"), "GSE23515 – Zhao")
) +
  plot_layout(ncol = 3, guides = "collect")

big_plot <- plots +
  plot_annotation(
    title = "Predicted Radiation Dose Across Independent Datasets",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5
      )
    )
  ) &
  theme(
    legend.position = "top",
    legend.box = "horizontal"
  )

legend_plot <- plot_dose_predictions(
  get("GSE63952_wang"),
  "GSE63952 – Wang"
) +
  theme(legend.position = "top")

shared_legend <- cowplot::get_legend(legend_plot)

plots <- (
  plot_dose_predictions(get("GSE63952_wang"), "GSE63952 – Wang") +
    plot_dose_predictions(get("GSE63952_li"), "GSE63952 – Li") +
    plot_dose_predictions(get("GSE63952_zhao"), "GSE63952 – Zhao") +
    plot_dose_predictions(get("GSE102971_wang"), "GSE102971 – Wang") +
    plot_dose_predictions(get("GSE102971_li"), "GSE102971 – Li") +
    plot_dose_predictions(get("GSE102971_zhao"), "GSE102971 – Zhao") +
    plot_dose_predictions(get("GSE23515_wang"), "GSE23515 – Wang") +
    plot_dose_predictions(get("GSE23515_li"), "GSE23515 – Li") +
    plot_dose_predictions(get("GSE23515_zhao"), "GSE23515 – Zhao")
) +
  plot_layout(ncol = 3) &
  theme(legend.position = "none")

big_plot <- cowplot::plot_grid(
  cowplot::ggdraw() +
    cowplot::draw_label(
      "Predicted Radiation Dose Across Independent Datasets",
      fontface = "bold",
      size = 16,
      x = 0.5,
      hjust = 0.5
    ),
  shared_legend,
  cowplot::plot_grid(plotlist = list(plots), ncol = 1),
  ncol = 1,
  rel_heights = c(0.08, 0.08, 1)
)

big_plot