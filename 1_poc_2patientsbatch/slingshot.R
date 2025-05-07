library('slingshot')
library('anndataR')

scanorama <- anndataR::read_h5ad("1_poc_2patientsbatch/data/scanorama.h5ad")
sce <- anndataR::to_SingleCellExperiment(scanorama)
scanorama_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = scanorama$X))

scanorama_result <- slingshot::slingshot(sce, cl = scanorama$obs$milestones, reducedDim = "X_scanorama", start.clus = "sB", end.clus = c("sEndD", "sEndC"))

colData(scanorama_result)$slingshot <- c()
sc_ad_res <- anndataR::from_SingleCellExperiment(scanorama_result)

anndataR::write_h5ad(scanorama_result, "1_poc_2patientsbatch/data/scanorama_slingshot.h5ad", mode = "w")
