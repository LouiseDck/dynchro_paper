library(dyngen)
library(dynwrap)
library(anndata)
library(dplyr)
library(wrapr)

source("1_poc_2patientsbatch/0_diverging_kinetics.R")

set.seed(584)
DEBUG <- FALSE

#' Generates two dyngen datasets with a batch effect introduced by modifying the kinetics of the second dataset.
#' @param index A string to append to the filenames when saving the datasets.
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

  # Generate the model up to the feature network, the common parts
  model_common <-
    config |>
    generate_tf_network() |>
    generate_feature_network()

  # Generate kinetics and gold standard and cells twice, resulting in two models with slightly different cells
  model_a <- model_common |>
    generate_kinetics() |>
    generate_gold_standard() |>
    generate_cells()

  model_b <- model_common |>
    generate_kinetics() |>
    generate_gold_standard() |>
    generate_cells()

  # Generate the experiment for model_a and convert to dyno dataset
  model_a <- model_a |> generate_experiment()
  dataset_a <- model_a |> as_dyno()

  # convert to anndata and save objects
  adata_ds_a <- convert_anndata(model_a, dataset_a)
  save_objects(dataset_a, adata_ds_a, paste0("1_poc_2patientsbatch/data/dataseta", index))

  # Generate the experiment for model b and convert to dyno dataset
  # Use diverging kinetics to introduce batch effect
  model_between2 <- generate_diverging_kinetics(model_a, model_b, 0.99)
  dataset_between <- model_between2 |> as_dyno()

  adata_ds_btwn <- convert_anndata(model_between2, dataset_between)
  save_objects(dataset_between, adata_ds_btwn, paste0("1_poc_2patientsbatch/data/datasetb", index))
}

#' Helper function to convert dyngen model and dataset to anndata format, keeping some metadata.
#' @param dyngen_model A dyngen model object.
#' @param dyngen_dataset A dyngen dataset object.
#' @return An anndata object with the expression data and metadata, such as model, milestones and root cell.
convert_anndata <- function(dyngen_model, dyngen_dataset) {
  print(dyngen_dataset$milestone_percentages)
  sA_cells <- dyngen_dataset$milestone_percentages |> dplyr::filter(milestone_id == "sA")
  root <- sA_cells[wrapr::orderv(sA_cells[, "percentage"]), ][1, 1]$cell_id

  adata_ds <- dyngen::as_anndata(dyngen_model)
  adata_ds$obs[["model"]] <- dyngen_dataset$cell_info[["model"]]
  adata_ds$obs[["milestones"]] <- dyngen_dataset$progressions$to
  adata_ds$uns[["iroot"]] <- root
  adata_ds
}

save_objects <- function(dataset, anndata, path) {
  saveRDS(
    dataset,
    paste0(path, ".rds")
  )

  anndata::write_h5ad(
    anndata,
    paste0(path, ".h5ad")
  )
}

generate_batch_effect(2)
