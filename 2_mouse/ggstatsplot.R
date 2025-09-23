library(ggstatsplot)
library(ggplot2)

gfi1_data <- readRDS("2_mouse/dynchro_bootstrap_gfi1_data.rds")
irf8_data <- readRDS("2_mouse/dynchro_bootstrap_irf8_data.rds")

p1 <- ggwithinstats(gfi1_data, x = branch, y = normalised_distance, title = "GFI1 knockout", xlab = "Branch", ylab = "Normalised distance",
        # boxplot.args = list(width = 0),
        boxplot.args = list(width = 0),
        results.subtitle = FALSE,
        centrality.plotting = FALSE,
        pairwise.display = "significant"
        ) + theme_classic() + theme(legend.position="none")
ggsave("2_mouse/gfi1_knockout.pdf", p1, width = 4.5, height = 4.5, dpi = 300)

p2 <- ggwithinstats(irf8_data, x = branch, y = normalised_distance, title = "IRF8 knockout", xlab = "Branch", ylab = "Normalised distance",
        boxplot.args = list(width = 0),
        results.subtitle = FALSE,
        centrality.plotting = FALSE,
        pairwise.display = "significant"
        )+ theme_classic() + theme(legend.position="none")

ggsave("2_mouse/irf8_knockout.pdf", p2, width = 4.5, height = 4.5, dpi = 300)
