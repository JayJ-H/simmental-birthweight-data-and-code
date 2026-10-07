# Preserve the submitted compositions; update only the revised analyses.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(tail(args, 1), winslash = "/")
base <- file.path(root, "revision_analysis")
project <- file.path(root, "original_analysis")
out <- file.path(root, "reproduced_figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
library(ggplot2)
library(patchwork)
stopifnot(requireNamespace("ragg"), requireNamespace("ggrepel"))
FIGURE_DIR <- file.path(project, "outputs", "figures")
read_csv_safe <- function(path) read.csv(path, check.names = FALSE)

# Evaluate only style definitions and the four unchanged upper panels.
src <- readLines(file.path(project, "scripts", "17_make_final_grouped_figures.R"),
                 warn = FALSE, encoding = "UTF-8")
eval(parse(text = src[grep("^pal <-", src):(grep("^save_final <-", src)-1)]))
tail_cols <- c("Lower tail" = pal[["lower"]], "Central range" = pal[["normal"]],
               "Upper tail" = pal[["upper"]])
chunk <- src[grep("^tail_stack_year <-", src):(grep("^qcoef <-", src)-1)]
chunk <- gsub("factor(tail_class, levels = names(tail_cols))",
              "factor(sub('Reference range', 'Central range', tail_class), levels = names(tail_cols))",
              chunk, fixed = TRUE)
set.seed(20261006)
eval(parse(text = chunk))
tail_stack_parity$parity_group <- factor(tail_stack_parity$parity_group,
                                        levels = c(">=7", "1", "2", "3", "4", "5", "6"))
p4B$data <- tail_stack_parity

qcoef <- read.csv(file.path(base, "results", "quantile_regression_cluster_bootstrap.csv"))
qcoef <- subset(qcoef, grepl("calf_sex|parity_group|calf_birth_season", term))
qcoef$label <- sub("parity_group", "Parity ", qcoef$term)
qcoef$label <- sub("calf_birth_season", "", qcoef$label)
qcoef$label[qcoef$term == "calf_sexmale"] <- "Male vs female"
qcoef$label[grepl("^Parity", qcoef$label)] <- paste0(qcoef$label[grepl("^Parity", qcoef$label)], " vs 1")
qcoef$label[grepl("^calf_birth_season", qcoef$term)] <-
  paste0(qcoef$label[grepl("^calf_birth_season", qcoef$term)], " vs spring")
qcoef$label <- factor(qcoef$label, levels = rev(unique(qcoef$label)))
p4E <- ggplot(qcoef, aes(x = estimate, y = label, colour = factor(tau))) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.25, colour = pal[["mid"]]) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y",
                width = 0, linewidth = 0.36, position = position_dodge(width = 0.5)) +
  geom_point(size = 1.1, position = position_dodge(width = 0.5)) +
  labs(x = "Quantile coefficient (kg)", y = NULL, colour = "Tau",
       title = "E  Quantile-regression effects") +
  theme_final(base_size = 5.7)

sens <- read.csv(file.path(base, "results", "effect_comparisons.csv"))
sens$scenario <- sub("_(lower|upper)$", "", sens$model)
scenario_labels <- c(p5p95 = "P5/P95", p15p85 = "P15/P85",
                     restricted = "25-65 kg only", complete = "2023-2025 only")
print(unique(sens$scenario))
keep <- grepl("p5p95|p15p85|weight_25_65|complete", sens$scenario)
sens <- sens[keep, ]
sens$scenario_label <- ifelse(grepl("p5p95", sens$scenario), "P5/P95",
  ifelse(grepl("p15p85", sens$scenario), "P15/P85",
  ifelse(grepl("complete|2023", sens$scenario), "2023-2025 only", "25-65 kg only")))
stopifnot(length(unique(sens$scenario_label)) == 4, nrow(sens) == 94)
sens$scenario_label <- factor(sens$scenario_label,
                             levels = c("2023-2025 only", "25-65 kg only", "P15/P85", "P5/P95"))
sens$tail <- ifelse(sens$outcome == "lower", "Lower tail", "Upper tail")
p4F <- ggplot(sens, aes(x = OR_ratio, y = scenario_label, colour = tail, shape = tail)) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.25, colour = pal[["mid"]]) +
  geom_point(size = 1.15, alpha = 0.8,
             position = position_jitterdodge(jitter.width = 0, jitter.height = 0.08,
                                             dodge.width = 0.5, seed = 20261006)) +
  scale_colour_manual(values = c("Lower tail" = pal[["lower"]], "Upper tail" = pal[["upper"]])) +
  scale_shape_manual(values = c("Lower tail" = 16, "Upper tail" = 17)) +
  labs(x = "Sensitivity OR / primary OR", y = NULL,
       title = "F  Sensitivity effect comparisons") +
  theme_final(base_size = 5.7)
fig <- (p4A | p4B) / (p4C | p4D) / (p4E | p4F) +
  plot_layout(heights = c(1, 1, 1.15))
name <- "Figure_3_Distribution_quantile_regression_threshold_stability.png"
ragg::agg_png(file.path(out, name), width = 183, height = 185, units = "mm", res = 450)
print(fig)
dev.off()
ragg::agg_png(file.path(out, "Figure_3_preview.png"), width = 1098, height = 1110, res = 152.4)
print(fig)
dev.off()
write.csv(qcoef, file.path(out, "Figure_3E_source.csv"), row.names = FALSE)
write.csv(sens, file.path(out, "Figure_3F_source.csv"), row.names = FALSE)

original <- file.path(root, "figures")
for (n in c("Figure_2_Temporal_distribution_and_tail_dynamics.png",
            "Figure_4_Exploratory_response_and_label_associations.png")) {
  stopifnot(file.copy(file.path(original, n), file.path(out, n), overwrite = TRUE))
}
stopifnot(file.copy(file.path(root, "figures", "Figure_1_Study_design_and_analytical_workflow.png"),
                   file.path(out, "Figure_1_Study_design_and_analytical_workflow.png"), overwrite = TRUE))
message("Four main figures restored; original Figure 1, 2 and 4 image bytes preserved.")
