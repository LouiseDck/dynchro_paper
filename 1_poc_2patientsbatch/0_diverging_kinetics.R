library(dyngen)

#' Generates a batch effect by modifying the kinetics of model 2, based on the alpha parameter.
#' @param model1 The first dyngen model.
#' @param model2 The second dyngen model.
#' @param alpha A numeric value between 0 and 1 indicating the weight of model2 in the kinetics modification. The closer alpha is to 1, the more the kinetics will resemble those of model2.
#' @return A new dyngen model with modified kinetics.
generate_diverging_kinetics <- function(model1, model2, alpha) {
  model_between <- model1
  fii <- model1$feature_info
  fni <- model1$feature_network

  genes <- c()

  for (i in seq(from = 1, to = nrow(fni))) {
    if (grepl("D[0-9]*", fni[i, "from"])) {
      genes <- c(genes, fni[i, ]$from)
      genes <- c(genes, fni[i, ]$to)
      genes <- unique(genes)
      for (col in colnames(fii)) {
        if (
          is.numeric(model1$feature_network[[col]]) &&
            !all(model1$feature_network[[col]] == model2$feature_network[[col]])
        ) {
          fii[[i, col]] <- model1$feature_network[[i, col]] *
            (1 - alpha) +
            model2$feature_network[[i, col]] * alpha
        }
      }
    }
  }

  for (i in seq(from = 1, to = nrow(fii))) {
    if (fii[i, ]$feature_id %in% genes) {
      for (col in colnames(fii)) {
        if (
          is.numeric(model1$feature_info[[col]]) &&
            !all(model1$feature_info[[col]] == model2$feature_info[[col]])
        ) {
          fii[[i, col]] <- model1$feature_info[[i, col]] * (1 - alpha) +
                           model2$feature_info[[i, col]] * alpha
        }
      }
    }
  }

  model_between$feature_info <- fii
  model_between$feature_network <- fni

  model_between <- model_between |>
    generate_gold_standard() |>
    generate_cells() |>
    generate_experiment()

  model_between
}
