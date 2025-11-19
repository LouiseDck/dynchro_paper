library(batchelor)
library(SingleCellExperiment)

library(anndataR)

d1 <- anndataR::read_h5ad("1_poc_2patientsbatch/data/fastmnn1.h5ad", as = "HDF5AnnData")
d2 <- anndataR::read_h5ad("1_poc_2patientsbatch/data/fastmnn2.h5ad", as = "HDF5AnnData")
rownames_1 <- rownames(d1$obs)
rownames_2 <- paste0(rownames(d2$obs), "_2")


result <- batchelor::fastMNN(t(d1$X), t(d2$X))
colnames(result) <- c(rownames_1, rownames_2)

assay(result, "reconstructed") <- as.matrix(assay(result, "reconstructed"))



ad <- anndataR::as_AnnData(
  result,
  x_mapping = "reconstructed",
  layers_mapping = c(),
  var_mapping = c(),
  uns_mapping = c()
)
ad$var <- c()
ad$uns <- c()
ad$write_h5ad("1_poc_2patientsbatch/data/fastmnn_result.h5ad", mode = "w")
