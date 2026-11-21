
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ==================== 1) Input / Output ====================
MTX_TSV <- "path/to/genexspot.tsv"
OUTDIR <- "path/to/output"
TARGET_GENE <- "GeneA"
dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)

# ==================== 2) Plot parameters ====================
V_MAX_EXPRESSION <- NULL  # NULL = automatic; e.g. 100 = fixed upper limit
POINT_SIZE <- 1.58
BG_COLOR <- "white"

# ==================== 3) Extract target gene expression ====================
mtx <- fread(MTX_TSV,sep="\t",header=TRUE,check.names=FALSE)
spots <- names(mtx)[-1]
expression <- as.numeric(unlist(mtx[mtx[[1]]==TARGET_GENE,-1,with=FALSE],use.names=FALSE))

sp <- strsplit(spots,"x",fixed=TRUE)
df <- data.frame(
  spot=spots,
  A=as.numeric(vapply(sp,`[[`,"",1)),
  B=as.numeric(vapply(sp,`[[`,"",2)),
  Expression=expression,
  stringsAsFactors=FALSE
)

# ==================== 4) Color scale ====================
vmax <- if(is.null(V_MAX_EXPRESSION)) max(df$Expression,na.rm=TRUE) else V_MAX_EXPRESSION

lab_k <- function(x) {
  ifelse(x==0,"0",ifelse(x%%1000==0,paste0(x/1000,"k"),format(x,trim=TRUE,scientific=FALSE)))
}

# ==================== 5) Spatial gene expression plot ====================
df$Expression_plot <- pmax(0,pmin(vmax,df$Expression))

p <- ggplot(df,aes(x=A,y=B,color=Expression_plot)) +
  geom_point(shape=16,size=POINT_SIZE) +
  scale_color_gradientn(
    colours=c("#2A7F9E","#D8D9D4","#B44A3C"),
    limits=c(0,vmax),
    breaks=c(0,vmax),
    labels=lab_k,
    oob=scales::squish,
    name=NULL
  ) +
  scale_y_reverse() +
  coord_equal(expand=TRUE,clip="off") +
  labs(x=NULL,y=NULL,title=TARGET_GENE) +
  theme_void() +
  theme(
    plot.background=element_rect(fill=BG_COLOR,colour=NA),
    panel.background=element_rect(fill=BG_COLOR,colour=NA),
    plot.title=element_text(hjust=0.5,size=12),
    legend.title=element_blank(),
    legend.text=element_text(size=12),
    legend.key.height=grid::unit(10,"pt"),
    legend.key.width=grid::unit(6,"pt"),
    legend.position="right",
    plot.margin=margin(t=10,r=18,b=10,l=10)
  ) +
  guides(color=guide_colorbar(barheight=grid::unit(38,"mm"),barwidth=grid::unit(3,"mm")))

# ==================== 6) Export PNG ====================
ggsave(file.path(OUTDIR,paste0("spatial_",TARGET_GENE,".png")),plot=p,width=4.2,height=3.6,units="in",dpi=600)
