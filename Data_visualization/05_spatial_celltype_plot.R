suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(ragg)
})

# ==================== 1) Input / Output ====================
INPUT_QS <- "path/to/annotation_output/host_annotation.qs"
OUTDIR <- "path/to/spatial_output"
dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)

obj_merged <- qread(INPUT_QS)

# ==================== 2) Extract spatial metadata ====================
df_all <- obj_merged@meta.data %>%
  tibble::as_tibble(rownames="spot") %>%
  dplyr::transmute(spot,sample=sampleNames,X=as.numeric(.data$X),Y=as.numeric(.data$Y),ctype=.data[["celltype"]]) %>%
  dplyr::filter(is.finite(X) & is.finite(Y) & !is.na(sample) & !is.na(ctype))

# ==================== 3) Cell type colors ====================
gut_palette <- c(
  "TA"         ="#E47B78",
  "Colonocyte" ="#A688C2",
  "Neuron"     ="#7DBB82",
  "SMC"        ="#69A8CF",
  "B"          ="#D9BE62"
)

ctype_order <- c("Colonocyte","TA","SMC","Neuron")
df_all$ctype <- factor(as.character(df_all$ctype),levels=ctype_order)

# ==================== 4) Spatial plotting function ====================
plot_one_sample <- function(dat,sample_name) {
  ggplot(dat,aes(x=X,y=Y,color=ctype)) +
    geom_point(shape=16,size=1.3,alpha=1) +
    scale_color_manual(values=gut_palette,breaks=ctype_order,drop=TRUE,name="Cell type") +
    scale_y_reverse() +
    coord_equal() +
    labs(title=sample_name) +
    theme_minimal(base_size=9) +
    theme(
      plot.title=element_text(hjust=0.5,face="bold",size=9),
      legend.position="right",
      legend.title=element_text(size=9,face="bold"),
      legend.text=element_text(size=9),
      axis.text=element_blank(),
      axis.title=element_blank(),
      axis.ticks=element_blank(),
      panel.grid=element_blank(),
      panel.background=element_blank(),
      plot.background=element_rect(fill="white",colour=NA)
    )
}

# ==================== 5) Export PNG ====================
samples <- sort(unique(df_all$sample))

for(s in samples) {
  dat_s <- dplyr::filter(df_all,sample==s)
  p <- plot_one_sample(dat_s,s)
  ggsave(file.path(OUTDIR,sprintf("spatial_%s.png",s)),plot=p,width=4,height=4,units="in",dpi=600,device=ragg::agg_png)
  message("[DONE] ",s," | PNG")
}
