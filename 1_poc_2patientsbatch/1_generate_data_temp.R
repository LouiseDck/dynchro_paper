library(tidyverse)
library(dyngen)
library(dyno)

source("dyngen_pooled/diverging_kinetics.R")

set.seed(42)
DEBUG <- FALSE

#############################
# 2. batch effect
#############################

generate_batch_effect <- function(index) {
  backbone <- backbone_bifurcating()

  config <-
    initialise_model(
      backbone = backbone,
      num_cells = 1000,
      num_tfs = nrow(backbone$module_info),
      num_targets = 250,
      num_hks = 250,
      simulation_params = simulation_default(
        census_interval = 10,
        ssa_algorithm = ssa_etl(tau = 300 / 3600),
        experiment_params = simulation_type_wild_type(num_simulations = 100)
      )
    )

  if (DEBUG) {
    config <-
      initialise_model(
        backbone = backbone,
        num_cells = 1000,
        num_tfs = nrow(backbone$module_info),
        num_targets = 50,
        num_hks = 50,
        verbose = interactive(),
        download_cache_dir = tools::R_user_dir("dyngen", "data"),
        simulation_params = simulation_default(
          census_interval = 5,
          ssa_algorithm = ssa_etl(tau = .01),
          experiment_params = simulation_type_wild_type(num_simulations = 10)
        )
      )
  }

  model_common <-
    config %>%
    generate_tf_network() %>%
    generate_feature_network()

  model_a <- model_common %>%
    generate_kinetics() %>%
    generate_gold_standard() %>%
    generate_cells()

  model_b <- model_common %>%
    generate_kinetics() %>%
    generate_gold_standard() %>%
    generate_cells()

  model_between <- generate_diverging_kinetics(model_a, model_b, 0.4)

  model_a <- model_a %>% generate_experiment()
  dataset_a <- model_a %>% as_dyno()

  a1 <- dataset_a$milestone_percentages %>% filter(milestone_id == "sA")
  a1 <- a1[order(a1[, "percentage"]), ][1, 1]$cell_id
  adata_ds <- as_anndata(model_a)
  adata_ds$obs[["model"]] <- dataset_a$cell_info[["model"]]
  adata_ds$obs[["milestones"]] <- dataset_a$progressions$to
  adata_ds$uns[["iroot"]] <- a1
  anndata::write_h5ad(
    adata_ds,
    paste0("dyngen_pooled/data/batcheffect_dataseta", index, ".h5ad")
  )

  saveRDS(
    dataset_a,
    paste0("dyngen_pooled/data/batcheffect_dataseta", index, ".rds")
  )

  model_between2 <- generate_diverging_kinetics(model_a, model_b, 0.99)

  dataset_between <- model_between2 %>% as_dyno()

  b1 <- dataset_between$milestone_percentages %>% filter(milestone_id == "sA")
  b1 <- b1[order(b1[, "percentage"]), ][1, 1]$cell_id

  adata_btwn <- as_anndata(model_between2)
  adata_btwn$obs[["model"]] <- dataset_between$cell_info[["model"]]
  adata_btwn$obs[["milestones"]] <- dataset_between$progressions$to
  adata_btwn$uns[["iroot"]] <- b1
  anndata::write_h5ad(
    adata_btwn,
    paste0("dyngen_pooled/data/batcheffect99_datasetb", index, ".h5ad")
  )
  saveRDS(
    dataset_between,
    paste0("dyngen_pooled/data/batcheffect99_datasetb", index, ".rds")
  )
}

generate_batch_effect(0)
generate_batch_effect(1)
generate_batch_effect(2)


model_comb2 <-
  combine_models(list(common = model_a, KO = model_between2)) %>%
  generate_experiment()
data_comb2 <- as_dyno(model_comb2)
plot_heatmap(data_comb2, features_oi = 50)


saveRDS(dataset_between, "dyngen_pooled/data/batcheffect99_datasetb.rds")
# dataset_between <- readRDS("dyngen_pooled/data/batcheffect9_datasetb.rds")
thing <- anndata::read_h5ad("dyngen_pooled/data/batcheffect99_datasetb.h5ad")

dataset_between <- readRDS("dyngen_pooled/data/batcheffect99_datasetb.rds")
plot_heatmap(dataset_between, features_oi = 50)


# ####
# # RUN TI methods
# ####
# part <- dataset_all$milestone_percentages %>% filter(milestone_id == "left_sA")
# id <- part[order(part[,"percentage"]),][1,1]$cell_id
# dataset_all <- add_prior_information(dataset_all, start_id = "cell838")
# dataset_all <- infer_trajectory(dataset_all, ti_paga())
# plot_dimred(dataset_all, dimred=dyndimred::dimred_umap, label_milestones = T)

# run_ti <- function(dataset, mst_id){
#   part <- dataset$milestone_percentages %>% filter(milestone_id == mst_id)
#   id <- part[order(part[,"percentage"]),][1,1]$cell_id
#   dataset <- add_prior_information(dataset, start_id = id)
#   infer_trajectory(dataset, ti_paga())
# }

# dataset_all <- run_ti(dataset_all, "left_sA")
# plot_dimred(dataset_all, dimred=dyndimred::dimred_umap, label_milestones = T)
# dataset_a <- run_ti(dataset_a, "sA")
# plot_dimred(dataset_a, dimred=dyndimred::dimred_umap, label_milestones = T)

## This batch effect is insufficient
# Trying a combo of:
#   1. switching genes in the feature_network
#   2. knockout certain genes (not completely)
DEBUG <- TRUE
backbone <- backbone_bifurcating()

config <-
  initialise_model(
    backbone = backbone,
    num_cells = 1000,
    num_tfs = nrow(backbone$module_info) * 3,
    num_targets = 250,
    num_hks = 250,
    simulation_params = simulation_default(
      census_interval = 10,
      ssa_algorithm = ssa_etl(tau = 300 / 3600),
      experiment_params = simulation_type_wild_type(num_simulations = 100)
    )
  )

if (DEBUG) {
  config <-
    initialise_model(
      backbone = backbone,
      num_cells = 1000,
      num_tfs = nrow(backbone$module_info),
      num_targets = 50,
      num_hks = 50,
      verbose = interactive(),
      download_cache_dir = tools::R_user_dir("dyngen", "data"),
      simulation_params = simulation_default(
        census_interval = 5,
        ssa_algorithm = ssa_etl(tau = 300 / 3600),
        experiment_params = simulation_type_wild_type(num_simulations = 10)
      )
    )
}

model_common <-
  config %>%
  generate_tf_network() %>%
  generate_feature_network() %>%
  generate_kinetics() %>%
  generate_gold_standard()

model_wt <- model_common %>%
  generate_cells()

plot_gold_mappings(model_wt, do_facet = FALSE)

backbone$module_info
plot_backbone_modulenet(model_common)
b3_genes <- model_common$feature_info %>%
  filter(module_id == "B3") %>%
  pull(feature_id)
plot_backbone_modulenet(model_common)

model_ko <- model_common
model_ko$simulation_params$experiment_params <- simulation_type_knockdown(
  num_simulations = 100L,
  timepoint = 0,
  genes = c("B6"),
  num_genes = 1,
  multiplier = 0
)
model_ko <- model_ko %>%
  generate_cells() %>%
  generate_experiment()
plot_gold_mappings(model_ko, do_facet = FALSE)

library(dyno)
dataset <- as_dyno(model_ko)
plot_dimred(dataset)

model_common <- model_common %>% generate_cells()

model_test1 <- model_common
ft_ids <- model_test1$feature_info[c(24, 26), ]$feature_id
ft_ids2 <- model_test1$feature_info[c(25, 28), ]$feature_id
model_test1$feature_info[c(25, 28), ]$feature_id <- ft_ids
model_test1$feature_info[c(24, 26), ]$feature_id <- ft_ids2

model_test1$simulation_params$experiment_params <- simulation_type_knockdown(
  num_simulations = 100L,
  timepoint = 0,
  genes = c("C1_TF1", "C2_TF1", "C3_TF1", "C4_TF1", "C5_TF1"),
  num_genes = 5,
  multiplier = c(0.3, 0.3, 0.3, 0.3, 0.3)
)

model_test1 <- model_test1 %>%
  # config %>%
  # generate_tf_network() %>%
  # generate_feature_network() %>%
  generate_kinetics() %>%
  generate_gold_standard() %>%
  generate_cells()

model_comb <-
  combine_models(list(common = model_common, KO = model_test1)) %>%
  generate_experiment()

plot_gold_mappings(model_comb, do_facet = FALSE)
data_comb <- as_dyno(model_comb)
plot_dimred(data_comb, color_cells = "milestone")

plot_heatmap(data_comb, features_oi = 50)
