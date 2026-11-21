
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(qs)
})

set.seed(123)

# ==================== 1) Input / Output ====================
sample_tbl <- data.frame(
  sampleNames=c("Sample1","Sample2"),
  file=c(
    "path/to/Sample1/genexspot.tsv",
    "path/to/Sample2/genexspot.tsv"
  ),
  stringsAsFactors=FALSE
)

OUTDIR <- "path/to/clustering_output"
dir.create(OUTDIR,showWarnings=FALSE,recursive=TRUE)

# ==================== 2) Read one sample ====================
read_one_sample <- function(tsv_path,sample_name) {
  dt <- fread(tsv_path,sep="\t",header=TRUE,check.names=FALSE)
  genes <- as.character(dt[[1]])
  mtx <- as.matrix(dt[,-1,with=FALSE])
  storage.mode(mtx) <- "numeric"
  rownames(mtx) <- genes
  mtx <- Matrix::Matrix(mtx,sparse=TRUE)
  
  obj <- CreateSeuratObject(counts=mtx,assay="RNA",min.cells=1,min.features=0)
  obj <- RenameCells(obj,new.names=paste0(sample_name,"#",colnames(obj)))
  obj$sampleNames <- sample_name
  obj[["percent.mt"]] <- PercentageFeatureSet(obj,pattern="^[Mm][Tt]-")
  
  spot_raw <- sub("^.*?#","",colnames(obj))
  sp <- strsplit(spot_raw,"x",fixed=TRUE)
  obj$X <- suppressWarnings(as.integer(vapply(sp,`[`,"",1)))
  obj$Y <- suppressWarnings(as.integer(vapply(sp,`[`,"",2)))
  obj
}

# ==================== 3) Read and merge samples ====================
objs <- lapply(seq_len(nrow(sample_tbl)),function(i) read_one_sample(sample_tbl$file[i],sample_tbl$sampleNames[i]))
names(objs) <- sample_tbl$sampleNames

obj_merged <- merge(x=objs[[1]],y=objs[-1])
rm(objs)
gc()

# ==================== 4) SCT + PCA ====================
DefaultAssay(obj_merged) <- "RNA"
obj_merged <- SCTransform(obj_merged,vars.to.regress="percent.mt",verbose=FALSE)
DefaultAssay(obj_merged) <- "SCT"
obj_merged <- RunPCA(obj_merged,assay="SCT",npcs=50,verbose=FALSE)

p_elbow <- ElbowPlot(obj_merged,ndims=50)
ggsave(file.path(OUTDIR,"ElbowPlot.png"),plot=p_elbow,width=6,height=4,dpi=600)

# ==================== 5) Harmony integration ====================
obj_merged <- IntegrateLayers(
  object=obj_merged,
  method=HarmonyIntegration,
  group.by.vars="sampleNames",
  normalization.method="SCT",
  orig.reduction="pca",
  new.reduction="harmony"
)

# ==================== 6) UMAP + Clustering ====================
obj_merged <- RunUMAP(obj_merged,reduction="harmony",dims=1:20)
obj_merged <- FindNeighbors(obj_merged, reduction = "harmony", dims = 1: 20)
obj_merged <- FindClusters(obj_merged,graph.name="SCT_snn",resolution=0.5)

# ==================== 7) UMAP plots ====================
p1 <- DimPlot(obj_merged,reduction="umap",group.by="sampleNames",label=TRUE,repel=TRUE,pt.size=0.1)
p2 <- DimPlot(obj_merged,reduction="umap",group.by="seurat_clusters",label=TRUE,repel=TRUE,pt.size=0.1)

ggsave(file.path(OUTDIR,"umap_by_sample.png"),plot=p1,width=7,height=6,dpi=600)
ggsave(file.path(OUTDIR,"umap_by_cluster.png"),plot=p2,width=7,height=6,dpi=600)

# ==================== 8) Export ====================
qsave(obj_merged,file.path(OUTDIR,"host_clustering.qs"))

