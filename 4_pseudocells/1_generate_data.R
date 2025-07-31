library(dyngen)
library(magrittr)
library(dynwrap)
library(dplyr)
library(wrapr)
library(anndata)

source("1_poc_2patientsbatch/0_diverging_kinetics.R")

set.seed(584)
DEBUG <- FALSE

#############################
# This generates 2 datasets with a 99 procent batch effect
#############################

generate_batch_effect <- function(index) {
    backbone <- backbone_linear()

    config <-
        initialise_model(
            backbone = backbone,
            num_cells = 5000,
            num_tfs = nrow(backbone$module_info),
            num_targets = 250,
            num_hks = 250,
            simulation_params = simulation_default(
                # total_time = 1000,
                # census_interval = 1, 
                # ssa_algorithm = ssa_etl(tau = 300/3600),
                experiment_params = simulation_type_wild_type(num_simulations = 45)
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

    # Generate the model up to a common part
    model_common <-
        config %>%
        generate_tf_network() %>%
        generate_feature_network()

    # Generate kinetics and gold standard and cells twice, resulting in slightly different cells
    model_a <- model_common %>%
        generate_kinetics() %>%
        generate_gold_standard() %>%
        generate_cells()

    model_b <- model_common %>%
        generate_kinetics() %>%
        generate_gold_standard() %>%
        generate_cells()

    # Generate the experiment for model a and convert to dyno dataset
    model_a <- model_a %>% generate_experiment()
    dataset_a <- model_a %>% as_dyno()

    # Set iroot, keep model info and milestones info and convert to anndata format
    a1 <- dataset_a$milestone_percentages %>% filter(milestone_id == "sA")
    a1 <- a1[orderv(a1[, "percentage"]), ][1, 1]$cell_id
    adata_ds <- as_anndata(model_a)
    adata_ds$obs[["model"]] <- dataset_a$cell_info[["model"]]
    adata_ds$obs[["milestones"]] <- dataset_a$progressions$to
    adata_ds$uns[["iroot"]] <- a1

    anndata::write_h5ad(
        adata_ds,
        paste0("4_pseudocells/data/dataseta", index, ".h5ad")
    )

    saveRDS(
        dataset_a,
        paste0("4_pseudocells/data/dataseta", index, ".rds")
    )

    # Generate the experiment for model b and convert to dyno dataset
    # Use diverging kinetics to introduce batch effect
    model_between2 <- generate_diverging_kinetics(model_a, model_b, 0.99)
    dataset_between <- model_between2 %>% as_dyno()

    # Set iroot, keep model info and milestones info and convert to anndata format
    b1 <- dataset_between$milestone_percentages %>% filter(milestone_id == "sA")
    b1 <- b1[orderv(b1[, "percentage"]), ][1, 1]$cell_id

    adata_btwn <- as_anndata(model_between2)
    adata_btwn$obs[["model"]] <- dataset_between$cell_info[["model"]]
    adata_btwn$obs[["milestones"]] <- dataset_between$progressions$to
    adata_btwn$uns[["iroot"]] <- b1
    anndata::write_h5ad(
        adata_btwn,
        paste0("4_pseudocells/data/datasetb", index, ".h5ad")
    )

    saveRDS(
        dataset_between,
        paste0("4_pseudocells/data/datasetb", index, ".rds")
    )
}

generate_batch_effect(2)
