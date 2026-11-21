
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ==================== 1) Input / Output ====================
MTX_TSV <- "path/to/genexspot.tsv"
OUTDIR <- "path/to/output"
dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)

# ==================== 2) Plot parameters ====================
V_MAX_UMI <- NULL  # NULL = automatic; e.g. 20000 = fixed upper limit
POINT_SIZE <- 1.58
BG_COLOR <- "white"

# ==================== 3) Calculate UMI per spot ====================
mtx <- fread(MTX_TSV,sep="\t",header=TRUE,check.names=FALSE)
spots <- names(mtx)[-1]
umi <- as.numeric(colSums(as.matrix(mtx[,-1,with=FALSE]),na.rm=TRUE))

sp <- strsplit(spots,"x",fixed=TRUE)
df <- data.frame(
  spot=spots,
  A=as.numeric(vapply(sp,`[`,"",1)),
  B=as.numeric(vapply(sp,`[`,"",2)),
  UMI=umi,
  stringsAsFactors=FALSE
)

# ==================== 4) Determine color scale ====================
vmax <- if(is.null(V_MAX_UMI)) max(df$UMI,na.rm=TRUE) else V_MAX_UMI
lab_k <- function(x) {
  ifelse(x==0,"0",ifelse(x%%1000==0,paste0(x/1000,"k"),format(x,trim=TRUE,scientific=FALSE)))
}

# ==================== 5) Spatial UMI plot ====================
df$UMI_plot <- pmax(0,pmin(vmax,df$UMI))

p <- ggplot(df,aes(x=A,y=B,color=UMI_plot)) +
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
  labs(x=NULL,y=NULL) +
  theme_void() +
  theme(
    plot.background=element_rect(fill=BG_COLOR,colour=NA),
    panel.background=element_rect(fill=BG_COLOR,colour=NA),
    legend.title=element_blank(),
    legend.text=element_text(size=12),
    legend.key.height=grid::unit(10,"pt"),
    legend.key.width=grid::unit(6,"pt"),
    legend.position="right",
    plot.margin=margin(t=10,r=18,b=10,l=10)
  ) +
  guides(color=guide_colorbar(barheight=grid::unit(38,"mm"),barwidth=grid::unit(3,"mm")))

# ==================== 6) Export PNG ====================
ggsave(file.path(OUTDIR,"ALL_spots_UMI_spatial.png"),plot=p,width=4.2,height=3.6,units="in",dpi=600)

