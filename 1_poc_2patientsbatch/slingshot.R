library('slingshot')
library('anndataR')
library('SingleCellExperiment')


#############
# SCANORAMA #
#############

scanorama <- anndataR::read_h5ad("1_poc_2patientsbatch/data/scanorama.h5ad")
sce <- anndataR::to_SingleCellExperiment(scanorama)
scanorama_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = scanorama$X))

# TODO don't use milestones, use the clusters
scanorama_result <- slingshot::slingshot(sce, cl = scanorama$obs$milestones, reducedDim = "X_scanorama", start.clus = "sB", end.clus = c("sEndD", "sEndC"))

colData(scanorama_result)$slingshot <- c()
sc_ad_res <- anndataR::as_AnnData(scanorama_result)

anndataR::write_h5ad(scanorama_result, "1_poc_2patientsbatch/data/scanorama_slingshot.h5ad", mode = "w")

###########
# SCGEN 1 #
###########

scgen1 <- anndataR::read_h5ad("1_poc_2patientsbatch/data/scgen1.h5ad")
sce_scgen1 <- anndataR::to_SingleCellExperiment(scgen1)

scgen1_result <- slingshot::slingshot(sce_scgen1, cl = scgen1$obs$leiden1, reducedDim = "corrected_latent", start.clus = 1, end.clus = c(7, 0))
colData(scgen1_result)$slingshot <- c()
scgen1_ad_res <- anndataR::from_SingleCellExperiment(scgen1_result)
anndataR::write_h5ad(scgen1_result, "1_poc_2patientsbatch/data/scgen1_slingshot.h5ad", mode = "w")

###########
# SCGEN 2 #
###########

scgen2 <- anndataR::read_h5ad("1_poc_2patientsbatch/data/scgen2.h5ad")
sce_scgen2 <- anndataR::to_SingleCellExperiment(scgen2)

scgen2_result <- slingshot::slingshot(sce_scgen2, cl = scgen2$obs$leiden1, reducedDim = "corrected_latent", start.clus = 0, end.clus = c(1, 5))
colData(scgen2_result)$slingshot <- c()
# scgen2_ad_res <- anndataR::from_SingleCellExperiment(scgen2_result)
anndataR::write_h5ad(scgen2_result, "1_poc_2patientsbatch/data/scgen2_slingshot.h5ad", mode = "w")


###########
# FASTMNN #
###########


fastmnn <- anndataR::read_h5ad("1_poc_2patientsbatch/data/fastmnn.h5ad")
sce_fastmnn <- anndataR::to_SingleCellExperiment(fastmnn)

# change of dist method cf https://github.com/kstreet13/slingshot/issues/87
fastmnn_result <- slingshot::slingshot(sce_fastmnn, cl = fastmnn$obs$leiden1, reducedDim = "X_fastmnn", start.clus = 3, end.clus = c(4, 1), dist.method = "mnn")    
colData(fastmnn_result)$slingshot <- c()

# fastmnn_ad_res <- anndataR::from_SingleCellExperiment(fastmnn_result)
anndataR::write_h5ad(fastmnn_result, "1_poc_2patientsbatch/data/fastmnn_slingshot.h5ad", mode = "w")

###########
# dynchro #
###########

dynchro1 <- anndataR::read_h5ad("1_poc_2patientsbatch/data/dataseta0_processed.h5ad")
sce_dynchro1 <- anndataR::to_SingleCellExperiment(dynchro1)
SingleCellExperiment::reducedDim(sce_dynchro1, "X_pca") <- dynchro1$obsm$X_pca
dynchro1_result <- slingshot::slingshot(sce_dynchro1, cl = dynchro1$obs$leiden, reducedDim = "X_pca", start.clus = 3, end.clus = c(4, 5))
colData(dynchro1_result)$slingshot <- c()
dynchro1_ad <- anndataR::as_AnnData(dynchro1_result)
anndataR::write_h5ad(dynchro1_result, "1_poc_2patientsbatch/data/dynchro1_slingshot.h5ad", mode = "w")

dynchro2 <- anndataR::read_h5ad("1_poc_2patientsbatch/data/datasetb0_processed.h5ad")
sce_dynchro2 <- anndataR::to_SingleCellExperiment(dynchro2)
SingleCellExperiment::reducedDim(sce_dynchro2, "X_pca") <- dynchro2$obsm$X_pca
dynchro2_result <- slingshot::slingshot(sce_dynchro2, cl = dynchro2$obs$leiden, reducedDim = "X_pca", start.clus = 1, end.clus = c(4, 3))
colData(dynchro2_result)$slingshot <- c()
anndataR::write_h5ad(dynchro2_result, "1_poc_2patientsbatch/data/dynchro2_slingshot.h5ad", mode = "w")
