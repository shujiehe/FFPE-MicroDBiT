
suppressPackageStartupMessages({
  library(ggplot2)
})

# ==================== 1) Input / Output ====================
tsv1 <- "path/to/metagenome_relative_abundance.tsv"
tsv2 <- "path/to/spatial_relative_abundance.tsv"
OUTDIR <- "path/to/output"
dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE)

label_x <- "Metagenome"
label_y <- "Spatial"

# ==================== 2) Read data ====================
df1 <- read.delim(tsv1,header=TRUE,sep="\t",stringsAsFactors=FALSE)
df2 <- read.delim(tsv2,header=TRUE,sep="\t",stringsAsFactors=FALSE)

df1 <- df1[,c("species","relative_abundance")]
df2 <- df2[,c("species","relative_abundance")]

df1$species <- trimws(df1$species)
df2$species <- trimws(df2$species)
df1$relative_abundance <- as.numeric(df1$relative_abundance)
df2$relative_abundance <- as.numeric(df2$relative_abundance)

# ==================== 3) Match species ====================
species_order <- df1$species
df2 <- df2[match(species_order,df2$species),c("species","relative_abundance")]

corr_df <- data.frame(
  species=species_order,
  metagenome=df1$relative_abundance,
  spatial=df2$relative_abundance,
  stringsAsFactors=FALSE
)

# ==================== 4) Pearson correlation ====================
ct <- cor.test(corr_df$metagenome,corr_df$spatial,method="pearson")

stats_df <- data.frame(
  n_species=nrow(corr_df),
  pearson_r=unname(ct$estimate),
  pearson_p=ct$p.value
)

write.csv(stats_df,file.path(OUTDIR,"metagenome_vs_spatial_correlation.stats.csv"),row.names=FALSE)

# ==================== 5) Correlation plot ====================
p <- ggplot(corr_df,aes(x=metagenome,y=spatial)) +
  geom_point(size=3.2,shape=16,colour="#B2182B") +
  geom_smooth(method="lm",se=FALSE,linewidth=0.8,linetype="dashed",colour="grey40") +
  scale_x_continuous(labels=percent_format(accuracy=1),expand=expansion(mult=c(0.02,0.08))) +
  scale_y_continuous(labels=percent_format(accuracy=1),expand=expansion(mult=c(0.02,0.08))) +
  labs(x=paste0(label_x," relative abundance"),y=paste0(label_y," relative abundance")) +
  theme_classic(base_size=12) +
  theme(
    plot.background=element_rect(fill="white",colour=NA),
    panel.background=element_rect(fill="white",colour=NA),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
    axis.title.x=element_text(size=12,margin=margin(t=10)),
    axis.title.y=element_text(size=12,margin=margin(r=10)),
    axis.line=element_line(linewidth=0.5),
    axis.ticks=element_line(linewidth=0.5),
    axis.ticks.length=grid::unit(3.5,"pt"),
    plot.margin=margin(t=10,r=18,b=10,l=14)
  )

# ==================== 6) Export PNG ====================
ggsave(file.path(OUTDIR,"metagenome_vs_spatial_correlation.png"),plot=p,width=4.8,height=4.4,units="in",dpi=600)

