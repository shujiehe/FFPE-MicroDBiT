
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ==================== 1) Input / Output ====================
sample_tbl <- data.frame(
  sampleNames=c("Sample1","Sample2","Sample3","Sample4"),
  file=c(
    "path/to/Sample1/genexspot.tsv",
    "path/to/Sample2/genexspot.tsv",
    "path/to/Sample3/genexspot.tsv",
    "path/to/Sample4/genexspot.tsv"
  ),
  stringsAsFactors=FALSE
)

OUTDIR <- "path/to/output"
dir.create(OUTDIR,showWarnings=FALSE,recursive=TRUE)

# ==================== 2) Calculate genes per spot ====================
calc_meta_one_sample <- function(tsv_path,sample_name) {
  dt <- fread(tsv_path,sep="\t",header=TRUE,check.names=FALSE)
  Genes <- colSums(as.matrix(dt[,-1,with=FALSE])>0)
  data.frame(sample=sample_name,Genes=as.numeric(Genes),log10Genes=log10(as.numeric(Genes)+1))
}

meta <- bind_rows(lapply(seq_len(nrow(sample_tbl)),function(i) calc_meta_one_sample(sample_tbl$file[i],sample_tbl$sampleNames[i])))
meta$sample <- factor(meta$sample,levels=sample_tbl$sampleNames)

# ==================== 3) Plot parameters ====================
cols <- setNames(c("#EFA5A5","#98BED9","#8DB18D","#E6D58C"),sample_tbl$sampleNames)
median_df <- meta %>% group_by(sample) %>% summarise(med=median(log10Genes),.groups="drop") %>% mutate(x=as.numeric(sample))

fig_theme <- function() {
  theme_classic(base_size=10) +
    theme(
      axis.text.x=element_text(size=10,angle=45,hjust=1,vjust=1),
      axis.text.y=element_text(size=10),
      axis.title.y=element_text(size=10),
      axis.line=element_line(linewidth=0.5),
      axis.ticks=element_line(linewidth=0.5),
      legend.position="none"
    )
}

# ==================== 4) Log10 Genes violin ====================
p_genes <- ggplot(meta,aes(x=sample,y=log10Genes,fill=sample)) +
  geom_violin(scale="width",trim=TRUE,width=0.6,linewidth=0.4,alpha=0.8) +
  geom_segment(data=median_df,aes(x=x-0.45,xend=x+0.45,y=med,yend=med),inherit.aes=FALSE,linetype="22",linewidth=0.5,colour="black") +
  geom_boxplot(width=0.075,outlier.size=0.5,alpha=1,linewidth=0.4,fill="white",colour="black") +
  scale_fill_manual(values=cols) +
  scale_y_continuous(breaks=scales::breaks_pretty(n=3)) +
  labs(x=NULL,y="Log10(Genes per pixel + 1)") +
  fig_theme()

ggsave(file.path(OUTDIR,"violin_log10_Genes.png"),plot=p_genes,width=4,height=4,dpi=600)

# ==================== 5) Genes statistics ====================
stats_tbl <- meta %>% group_by(sample) %>% summarise(n_spots=n(),Genes_mean=mean(Genes),Genes_median=median(Genes),.groups="drop")
write.csv(stats_tbl,file.path(OUTDIR,"per_sample_Genes_stats.csv"),row.names=FALSE)
