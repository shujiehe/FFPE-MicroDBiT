
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
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

# ==================== 2) Calculate UMI per spot ====================
calc_meta_one_sample <- function(tsv_path,sample_name) {
  dt <- fread(tsv_path,sep="\t",header=TRUE,check.names=FALSE)
  UMI <- colSums(as.matrix(dt[,-1,with=FALSE]))
  data.frame(sample=sample_name,UMI=as.numeric(UMI),log10UMI=log10(as.numeric(UMI)+1))
}

meta <- bind_rows(lapply(seq_len(nrow(sample_tbl)),function(i) calc_meta_one_sample(sample_tbl$file[i],sample_tbl$sampleNames[i])))
meta$sample <- factor(meta$sample,levels=sample_tbl$sampleNames)

# ==================== 3) Plot parameters ====================
cols <- setNames(c("#EFA5A5","#98BED9","#8DB18D","#E6D58C"),sample_tbl$sampleNames)
y_breaks_log_umi <- seq(floor(min(meta$log10UMI)),ceiling(max(meta$log10UMI)),by=1)

fig_theme <- function() {
  theme_classic(base_size=10) +
    theme(
      axis.text.x=element_text(size=10,angle=45,hjust=1,vjust=1),
      axis.text.y=element_text(size=10),
      axis.title.y=element_text(size=10),
      axis.line=element_line(linewidth=0.7),
      axis.ticks=element_line(linewidth=0.7),
      legend.position="none"
    )
}

# ==================== 4) Log10 UMI violin ====================
median_df <- meta %>% group_by(sample) %>% summarise(median_log10UMI=median(log10UMI),.groups="drop")
median_df$x <- as.numeric(median_df$sample)

p_log <- ggplot(meta,aes(x=sample,y=log10UMI,fill=sample)) +
  geom_violin(scale="width",trim=TRUE,width=0.6,linewidth=0.4,alpha=0.8) +
  geom_segment(data=median_df,aes(x=x-0.45,xend=x+0.45,y=median_log10UMI,yend=median_log10UMI),inherit.aes=FALSE,linetype="22",linewidth=0.5,colour="black") +
  geom_boxplot(width=0.075,outlier.size=0.5,alpha=1,linewidth=0.4,fill="white",colour="black") +
  scale_fill_manual(values=cols) +
  scale_y_continuous(breaks=y_breaks_log_umi) +
  labs(x=NULL,y="Log10(UMIs per pixel + 1)") +
  fig_theme()

ggsave(file.path(OUTDIR,"violin_log10_UMI.png"),plot=p_log,width=4,height=4,dpi=600)

# ==================== 5) UMI statistics ====================
stats_tbl <- meta %>% group_by(sample) %>% summarise(n_spots=n(),UMI_mean=mean(UMI),UMI_median=median(UMI),.groups="drop")
write.csv(stats_tbl,file.path(OUTDIR,"per_sample_UMI_stats.csv"),row.names=FALSE)
