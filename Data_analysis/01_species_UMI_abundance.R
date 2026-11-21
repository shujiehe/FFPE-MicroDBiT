
suppressPackageStartupMessages({
  library(data.table)
})

# 1) Input / Output
MTX_TSV <- "path/to/Expression_matrix.tsv"
SPECIES_LIST <- "path/to/species_list.txt"
OUTDIR <- "path/to/output"
dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)

# 2) Load matrix and species list
dt <- fread(MTX_TSV,sep="\t",header=TRUE,check.names=FALSE)
rn <- trimws(as.character(dt[[1]]))
sp_list <- unique(trimws(readLines(SPECIES_LIST,warn=FALSE)))
sp_list <- sp_list[nzchar(sp_list)]
sp_list <- sp_list[order(nchar(sp_list),decreasing=TRUE)]

# 3) Map genes to species
sp_hit <- rep(NA_character_,length(rn))
for(sp in sp_list) {
  idx <- is.na(sp_hit) & startsWith(rn,paste0(sp,"_"))
  sp_hit[idx] <- sp
}

# 4) Calculate species-level UMI abundance
gene_umi <- rowSums(as.matrix(dt[,-1,with=FALSE]))
res <- data.table(species=sp_hit,gene_umi=gene_umi)[,.(total_umi=sum(gene_umi)),by=species]
grand_total <- sum(res$total_umi)
res[,umi_fraction:=total_umi/grand_total]
res[,umi_percent:=umi_fraction*100]
setorder(res,-total_umi,species)

# 5) Export
fwrite(res,file.path(OUTDIR,"species_UMI_abundance.tsv"),sep="\t")
