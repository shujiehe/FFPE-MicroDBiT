
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

# ==================== 1) Input / Output ====================
sample_files <- c(
  "Sample1"="path/to/Sample1/Filtered_expmat.tsv",
  "Sample2"="path/to/Sample2/Filtered_expmat.tsv"
)
species_list_file <- "path/to/species_list.txt"
out_dir <- "path/to/output"
dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)

# ==================== 2) Read species list ====================
species_list <- fread(species_list_file,header=FALSE,sep="\t",data.table=FALSE)[,1]
species_list <- trimws(as.character(species_list))
species_list <- unique(species_list[!is.na(species_list) & species_list!=""])
species_list <- species_list[order(nchar(species_list),decreasing=TRUE)]

# ==================== 3) Read Gene/Contig x Spot matrix ====================
read_tsv_matrix <- function(path) {
  dt <- fread(path,sep="\t",header=TRUE,check.names=FALSE)
  contig_ids <- as.character(dt[[1]])
  dt[[1]] <- NULL
  mat <- as.matrix(dt)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- contig_ids
  Matrix(mat,sparse=TRUE)
}

# ==================== 4) Map Gene/Contig to species ====================
assign_species <- function(contig_ids,species_list) {
  uniq_contigs <- unique(contig_ids)
  mapped <- rep(NA_character_,length(uniq_contigs))
  for(sp in species_list) {
    hit <- is.na(mapped) & startsWith(uniq_contigs,paste0(sp,"_"))
    mapped[hit] <- sp
  }
  mapped[match(contig_ids,uniq_contigs)]
}

# ==================== 5) Shannon index ====================
calc_shannon <- function(x) {
  x <- as.numeric(x)
  x <- x[x>0]
  p <- x/sum(x)
  -sum(p*log(p))
}

# ==================== 6) Per-sample analysis ====================
all_spot_res <- list()
all_sample_res <- list()

for(sample_name in names(sample_files)) {
  mat <- read_tsv_matrix(sample_files[[sample_name]])
  cat("[INFO]",sample_name,"| Features:",nrow(mat),"| Spots:",ncol(mat),"\n")
  
  species_of_contig <- assign_species(rownames(mat),species_list)
  keep <- !is.na(species_of_contig)
  mat <- mat[keep,,drop=FALSE]
  species_of_contig <- species_of_contig[keep]
  
  species_levels <- unique(species_of_contig)
  group_index <- match(species_of_contig,species_levels)
  
  G <- sparseMatrix(
    i=group_index,
    j=seq_along(group_index),
    x=1,
    dims=c(length(species_levels),length(group_index)),
    dimnames=list(species_levels,rownames(mat))
  )
  
  species_by_spot <- G %*% mat
  spot_total_reads <- Matrix::colSums(species_by_spot)
  
  spot_shannon <- vapply(seq_len(ncol(species_by_spot)),function(j) calc_shannon(species_by_spot[,j]),numeric(1))
  
  spot_dt <- data.table(
    sample=sample_name,
    spot=colnames(species_by_spot),
    total_bacterial_reads=as.numeric(spot_total_reads),
    shannon=as.numeric(spot_shannon)
  )
  
  spot_dt[,x:=as.integer(sub("x.*$","",spot))]
  spot_dt[,y:=as.integer(sub("^.*x","",spot))]
  
  sample_species_reads <- Matrix::rowSums(species_by_spot)
  
  sample_dt <- data.table(
    sample=sample_name,
    total_bacterial_reads=as.numeric(sum(sample_species_reads)),
    shannon=calc_shannon(sample_species_reads)
  )
  
  fwrite(spot_dt,file.path(out_dir,paste0(sample_name,"_spot_shannon.tsv")),sep="\t")
  all_spot_res[[sample_name]] <- spot_dt
  all_sample_res[[sample_name]] <- sample_dt
}

# ==================== 7) Export combined results ====================
all_spot_dt <- rbindlist(all_spot_res,use.names=TRUE)
all_sample_dt <- rbindlist(all_sample_res,use.names=TRUE)

fwrite(all_spot_dt,file.path(out_dir,"all_samples_spot_shannon.tsv"),sep="\t")
fwrite(all_sample_dt,file.path(out_dir,"all_samples_sample_shannon.tsv"),sep="\t")

