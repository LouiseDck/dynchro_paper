library('slingshot')
library('anndataR')
library('SingleCellExperiment')

run_slingshot <- function(sce, ad, reducedDim) {
  ss_result <- slingshot::slingshot(
    sce,
    cl = ad$obs$leiden,
    reducedDim = reducedDim,
    start.clus = ad$uns$slingshot_start,
    end.clus = ad$uns$slingshot_end
  )
  colData(ss_result)$slingshot <- c()
  return(ss_result)
}


#############
# SCANORAMA #
#############

scanorama <- anndataR::read_h5ad("1_poc_2patientsbatch/data/scanorama.h5ad")
sce_scanorama <- anndataR::to_SingleCellExperiment(scanorama)

scanorama_res <- run_slingshot(sce_scanorama, scanorama, "X_scanorama")
anndataR::write_h5ad(
  scanorama_res,
  "1_poc_2patientsbatch/data/scanorama_slingshot.h5ad",
  mode = "w"
)

###########
# SCGEN 1 #
###########

scgen1 <- anndataR::read_h5ad("1_poc_2patientsbatch/data/scgen1.h5ad")
sce_scgen1 <- scgen1$as_SingleCellExperiment(scgen1)

scgen1_result <- slingshot::slingshot(
  sce_scgen1,
  cl = scgen1$obs$leiden1,
  reducedDim = "corrected_latent",
  start.clus = 1,
  end.clus = c(4, 3)
)
colData(scgen1_result)$slingshot <- c()
scgen1_ad_res <- anndataR::as_AnnData(scgen1_result)
anndataR::write_h5ad(
  scgen1_ad_res,
  "1_poc_2patientsbatch/data/scgen1_slingshot.h5ad",
  mode = "w"
)

###########
# SCGEN 2 #
###########

scgen2 <- anndataR::read_h5ad("1_poc_2patientsbatch/data/scgen2.h5ad")
sce_scgen2 <- scgen2$as_SingleCellExperiment()

scgen2_res <- run_slingshot(sce_scgen2, scgen2, "corrected_latent")

anndataR::write_h5ad(
  scgen2_res,
  "1_poc_2patientsbatch/data/scgen2_slingshot.h5ad",
  mode = "w"
)

###########
# FASTMNN #
###########

fastmnn <- anndataR::read_h5ad("1_poc_2patientsbatch/data/fastmnn.h5ad")
sce_fastmnn <- fastmnn$as_SingleCellExperiment()

# change of dist method cf https://github.com/kstreet13/slingshot/issues/87
fastmnn_result <- slingshot::slingshot(
  sce_fastmnn,
  cl = fastmnn$obs$leiden1,
  reducedDim = "X_fastmnn",
  start.clus = 3,
  end.clus = c(4, 1),
  dist.method = "mnn"
)
colData(fastmnn_result)$slingshot <- c()

# fastmnn_ad_res <- anndataR::from_SingleCellExperiment(fastmnn_result)
anndataR::write_h5ad(
  fastmnn_result,
  "1_poc_2patientsbatch/data/fastmnn_slingshot.h5ad",
  mode = "w"
)

##############
# preprocess #
##############

dynchro1 <- anndataR::read_h5ad(
  "1_poc_2patientsbatch/data/dataseta0_processed.h5ad"
)
slingshot_start <- dynchro1$uns$slingshot_start
slingshot_end <- dynchro1$uns$slingshot_end

sce_dynchro1 <- anndataR::to_SingleCellExperiment(dynchro1)
SingleCellExperiment::reducedDim(sce_dynchro1, "X_pca") <- dynchro1$obsm$X_pca
dynchro1_result <- slingshot::slingshot(
  sce_dynchro1,
  cl = dynchro1$obs$leiden,
  reducedDim = "X_pca",
  start.clus = slingshot_start,
  end.clus = slingshot_end
)
colData(dynchro1_result)$slingshot <- c()

anndataR::write_h5ad(
  dynchro1_result,
  "1_poc_2patientsbatch/data/dataset1_slingshot.h5ad",
  mode = "w"
)

dynchro2 <- anndataR::read_h5ad(
  "1_poc_2patientsbatch/data/datasetb0_processed.h5ad"
)
slingshot_start <- dynchro2$uns$slingshot_start
slingshot_end <- dynchro2$uns$slingshot_end

sce_dynchro2 <- anndataR::to_SingleCellExperiment(dynchro2)
SingleCellExperiment::reducedDim(sce_dynchro2, "X_pca") <- dynchro2$obsm$X_pca
dynchro2_result <- slingshot::slingshot(
  sce_dynchro2,
  cl = dynchro2$obs$leiden,
  reducedDim = "X_pca",
  start.clus = slingshot_start,
  end.clus = slingshot_end
)
colData(dynchro2_result)$slingshot <- c()
anndataR::write_h5ad(
  dynchro2_result,
  "1_poc_2patientsbatch/data/dataset2_slingshot.h5ad",
  mode = "w"
)
