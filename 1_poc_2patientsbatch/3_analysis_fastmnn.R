library(batchelor)
library(SingleCellExperiment)

library(anndataR)

d1 <- anndataR::read_h5ad("1_poc_2patientsbatch/data/fastmnn1.h5ad")
d2 <- anndataR::read_h5ad("1_poc_2patientsbatch/data/fastmnn2.h5ad")

result <- batchelor::fastMNN(t(d1$X), t(d2$X))
assay(result, "reconstructed") <- as.matrix(assay(result, "reconstructed"))

# anndataR::write_h5ad(result, "1_poc_2patientsbatch/data/fastmnn_result2.h5ad", mode = "w")

ad <- anndataR::from_SingleCellExperiment(
  result,
  x_mapping = "reconstructed",
  layers_mapping = c(),
  var_mapping = c(),
  uns_mapping = c()
)
ad$var <- c()
ad$uns <- c()
ad$write_h5ad("1_poc_2patientsbatch/data/fastmnn_result.h5ad", mode = "w")
