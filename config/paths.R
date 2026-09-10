# Single definition of the project data root for all R code in this repository.
#
# The scripts here read and write large data files (Giotto objects, Seurat/DESeq2
# inputs) that are not distributed with the code. Point ACC_DATA_ROOT at a local
# copy of the data tree, or edit the fallback below.
#
#     export ACC_DATA_ROOT=/path/to/acc_data

ACC_DATA_ROOT <- Sys.getenv("ACC_DATA_ROOT", path.expand("~/acc_data"))

# Reference resources the pipelines expect underneath the data root.
MSIGDB_GMT  <- file.path(ACC_DATA_ROOT, "reference", "msigdb_gmt")
CELLPHONEDB <- file.path(ACC_DATA_ROOT, "reference", "cellphonedb")
