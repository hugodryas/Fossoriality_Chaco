############################################################################################################################
### Functional and phylogenetic data explain distribution of snakes in the world's largest dry forest ecoregion, the Gran Chaco ###

# Hugo Cabral, Thaís B. Guedes, Diego Santana

#Programa de Pós-Graduação em Biologia Animal, Universidade Estadual Paulista, São José do Rio Preto, SP, Brazil
#Instituto de Investigación Biológica del Paraguay. Del Escudo 1607, Asunción, Paraguay
#Mapinguari - Laboratório de Biogeografia e Sistemática de Anfíbios e Répteis, Universidade Federal de Mato Grosso do Sul, 79002-970, Campo Grande, MS, Brazil
#Universidade Estadual do Maranhão, Programa de Pós-Graduação em Biodiversidade, Ambiente e Saúde, Caxias, Maranhão, 65604-380, Brazil
#University of Gothenburg, Gothenburg Global Biodiversity Center and Department of Biological and Environmental Sciences, Göteborg, SE 405 30, Sweden
#
# This script is arranged in the same order that the results are presented in the article
#Este script esta ordenado de la misma forma que los resultados son presentados en el articulo
# Metrics and PCA  # Metricas y PCA

# Apply the Null Model # Aplicamos un modelo nulo a nuestros resultados

# Analyze the geographical Patterns # Analizamos los patrones geograficos

# Drivers of vertical stratification # Analizamos las variables que influyen en la estratificación vertical

# Phylogenetic regionalization # Regionazliación filogenetica


rm(list=ls()); gc()
setwd("")

# Install and load R packages needed to run the analysis: # Instalamos y cargamos los paquetes 
needed_packages<-c("data.table", "stringr", "plyr", "readr", "viridis", "tidyverse", "ggnewscale", "foreach", "doParallel",
                   "ggpubr", "RColorBrewer", "ggplot2", "foreach", "cowplot", "dplyr" ,"egg", "gridExtra", 
                   "rgdal", "raster", "tidyr", "maptools", "lattice", "scales", "SpatialPack")
new.packages<-needed_packages[!(needed_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(needed_packages, require, character.only = TRUE)
rm(needed_packages,new.packages)

# STEP 1 - GET THE METRICS OF VERTICALITY PER SPECIES
############################################################################################################################
# STEP 1 - OBTENEMOS LAS METRICAS DE VERTICALIDAD PARA LAS ESPECIES
rm(list=ls())

# Load snake trait data and prepare verticality-related variables:
trait_data<-fread("trait_data.csv", stringsAsFactors=TRUE, encoding="UTF-8")

# Make 'Habitat' an ordinal variable
trait_data$HabitatNum<-trait_data$Habitat
trait_data$HabitatNum<-factor(trait_data$HabitatNum, levels = c("Arboreal", "Fossorial", "Semiarboreal", "Terrestrial", "aquatic", "semifossorial"),
                              labels = c(3, 1, 2.5, 2, 2, 1.5))
trait_data$HabitatNum<-as.numeric(as.character(trait_data$HabitatNum))

# Make 'EyeSize' an ordinal variable:
trait_data$EyeSizeNum<-trait_data$EyeSize
trait_data$EyeSizeNum<-factor(trait_data$EyeSizeNum, levels = c("large", "medium", "small"),
                              labels = c(3, 2, 1))
trait_data$EyeSizeNum<-as.numeric(as.character(trait_data$EyeSizeNum))

# Put 'BodyMass' in the log10 scale:
trait_data$BodyMass<-as.numeric(as.character(trait_data$BodyMass))
trait_data$LogBodyMass<-log10(trait_data$BodyMass)

# Get the proportion of total length that corresponds to the tail length:
trait_data$BodyLength<-as.numeric(as.character(trait_data$BodyLength))
trait_data$SVL<-as.numeric(as.character(trait_data$SVL))
trait_data$TL<-as.numeric(as.character(trait_data$TL))
trait_data$TailProp<-trait_data$TL/trait_data$BodyLength

# Get a dataframe with only the variables of interest:
trait_data_Num<-trait_data[, .(Binomial, LogBodyMass, EyeSizeNum, HabitatNum, TailProp)]
trait_data_Num_PCA<-trait_data_Num[,2:5]

# Get ortogonal axes to represent the 'trait space':
pca_scores<-(prcomp(trait_data_Num_PCA, center=T, scale=T))
# Get the percentage of variation explained by the first two PCA axis:
eigs <- pca_scores$sdev^2
axis1 <- round(eigs[1] / sum(eigs), 3)*100
axis2 <- round(eigs[2] / sum(eigs), 3)*100

# Get the loadings for the variables (arrows in the PCA plot):
var_loadings<-as.data.frame(pca_scores[[2]])
var_loadings$variable <- rownames(var_loadings)

# Get the scores for the species (points in the PCA plot):
pca_scores<-as.data.frame(pca_scores[[5]])
pca_scores$Genus<-gsub("([A-Za-z]+).*", "\\1", trait_data_Num$Binomial) # get only the genus from the binomial name
pca_scores$Epithet<-stringr::word(trait_data_Num$Binomial, 2)
pca_scores$Label<-paste0(stringr::str_extract(pca_scores$Genus, "^.{2}"), stringr::str_extract(pca_scores$Epithet, "^.{2}"))

# Get the range of the first two axes:
ylim<-c(min(pca_scores[,2]), max(pca_scores[,2]))
xlim<-c(min(pca_scores[,1]), max(pca_scores[,1]))

# Get a colour pallete to represent different levels of verticality:
myColors<-viridis_pal(direction=-1, option="cividis")(6) # 6 colours
myColors<-myColors[2:6] # discard the yellow shades (low contrast)

PCA_plot
# Plot PCA results for the two first unconstrained axes:
PCA_plot <- ggplot() + coord_cartesian(x = xlim, y = ylim) + # Adjust the xlim and ylim if necessary.
  geom_point(data=pca_scores, aes(x=PC1, y=PC2), shape=16, colour="gray50", size=0.5) +
  geom_text(data=pca_scores, aes(x=PC1, y=PC2, label=Label, colour=PC1), check_overlap = TRUE, size=3) +
  geom_segment(data=var_loadings, aes(x=0, xend=PC1*3, y=0, yend=PC2*3), arrow=arrow(length=unit(0.25, "cm")), colour="black") +
  geom_text(data=var_loadings, aes(x=PC1*4.75, y=PC2*3.2, label=variable), colour="black", size=4) + # adjust the arrows fit
  scale_colour_manual(values=myColors) +
  scale_colour_gradientn(colours=viridis_pal(direction=-1, option="plasma")(10), na.value = NA) +
  theme(panel.background = element_rect(fill="white"),
        axis.text = element_text(size=12),
        panel.border = element_rect(fill=NA, colour="black"),
        axis.title = element_text(size=14),
        panel.grid.minor = element_blank(), # remove minor gridlines
        panel.grid.major = element_blank(), # remove major gridlines
        legend.position="none") +
  labs(x=paste0("PCA 1", " (", axis1, "%)"), y=paste0("PCA 2", " (", axis2, "%)")); PCA_plot

# Export to disk:
ggsave("PCA_biplot.png", plot=PCA_plot, width=6, height=5, units="in", bg = "transparent", limitsize=F)
ggsave("PCA_biplot.pdf", plot=PCA_plot, width=6, height=5, units="in", bg = "transparent", limitsize=F)

# Add to the trait data, the score each species received in the first PCA axis.
trait_data_Num$PCA1<-pca_scores$PC1

#####

# STEP 2 - APPLY A NULL MODEL TO CONTROL THE INFLUENCE OF SPECIES RICHNESS ON THE VERTICALITY METRICS
############################################################################################################################
# STEP 2 - APLICAMOS EL MODELO NULO PARA CONTROLAR LA INFLUENCIA DE LA RIQUEZA DE ESPECIES EN NUESTRAS METRICAS DE VERTICALIDAD

# Remove all objects from the workspace, but 'trait_data_Num':
rm(list=setdiff(ls(), "trait_data_Num"))
trait_data_Num<-trait_data_Num[,c(1:6)]

# Create a function to get the geometric mean:
geometric.mean <- function(x,na.rm=TRUE){
  exp(mean(log(x), na.rm=na.rm)) }

# Load the assemblage-level data on snake species at 110 km of spatial resolution and filter the Chaco-related cells:
assemblage_data<-fread("", stringsAsFactors=TRUE, encoding="UTF-8")#presence/absence data.
selected_cells<-fread("", stringsAsFactors=TRUE, encoding="UTF-8")#grid to be used, a table with cell Id of each grid

# Filter the assemblage-data to reflect only the assemblages intersecting with the Chaco ecoregion:
assemblage_data<-assemblage_data[Cell_Id55 %in% as.character(selected_cells$Cell_Id55)]
assemblage_data<-droplevels(assemblage_data)

# Calculate the number of occurrences per species and add the information to the trait data:
assemblage_data<-as.data.table(assemblage_data)
gr_cells_per_spp<-assemblage_data[, .(n_occ = length(unique(Cell_Id55))), by = Binomial]
trait_data_Num<-merge(trait_data_Num, gr_cells_per_spp, by = 'Binomial', all.x = TRUE)
rm(gr_cells_per_spp)

trait_data_Num<-trait_data_Num[-which(is.na(trait_data_Num$n_occ)),]
trait_data_Num<-droplevels(trait_data_Num)

# Filter the assemblage-data to contain only species in the 'trait_data_Num' object: 
trait_data_Num<-droplevels(trait_data_Num)
assemblage_data<-assemblage_data[Binomial %in% as.character(trait_data_Num$Binomial)]
assemblage_data<-droplevels(assemblage_data)

# Configure the system to run null models in parallel:
cl<-makePSOCKcluster(detectCores()-1)
registerDoParallel(cl)
getDoParWorkers()

# Get null model assemblages and extract the average verticality metric in each assemblage:
NullValues_HabitatNum<-foreach(i = 1:1000, .combine = 'cbind', .packages = c("data.table"))  %dopar% {
  trait_data_Num<-trait_data_Num[, c(1:7)]
  trait_data_Num$NullMetric<-sample(trait_data_Num$HabitatNum, replace=F, prob=(1/(trait_data_Num$n_occ)))
  trait_data_Null<-trait_data_Num[, c(1,8)]
  null_assemblage_data<-merge(assemblage_data, trait_data_Null, by='Binomial', all.x=TRUE, allow.cartesian=TRUE)
  null_assemblage_data<-null_assemblage_data[, .(geometric.mean(NullMetric, na.rm=T)), by = .(Cell_Id55)]
  null_assemblage_data$V1 # null metric value for each assemblage
}

# Same as above, but for the TailProp metric:
NullValues_TailProp<-foreach(i = 1:1000, .combine = 'cbind', .packages = c("data.table"))  %dopar% {
  trait_data_Num<-trait_data_Num[, c(1:7)]
  trait_data_Num$NullMetric<-sample(trait_data_Num$TailProp, replace=F, prob=(1/(trait_data_Num$n_occ)))
  trait_data_Null<-trait_data_Num[, c(1,8)]
  null_assemblage_data<-merge(assemblage_data, trait_data_Null, by='Binomial', all.x=TRUE, allow.cartesian=TRUE)
  null_assemblage_data<-null_assemblage_data[, .(geometric.mean(NullMetric, na.rm=T)), by = .(Cell_Id55)]
  null_assemblage_data$V1 # null metric value for each assemblage
}

# Same as above, but for the metric based on PCA1-score:
NullValues_PCA1<-foreach(i = 1:1000, .combine = 'cbind', .packages = c("data.table"))  %dopar% {
  trait_data_Num<-trait_data_Num[, c(1:7)]
  trait_data_Num$NullMetric<-sample(trait_data_Num$PCA1, replace=F, prob=(1/(trait_data_Num$n_occ)))
  trait_data_Null<-trait_data_Num[, c(1,8)]
  null_assemblage_data<-merge(assemblage_data, trait_data_Null, by='Binomial', all.x=TRUE, allow.cartesian=TRUE)
  null_assemblage_data<-null_assemblage_data[, .(geometric.mean(NullMetric, na.rm=T)), by = .(Cell_Id55)]
  null_assemblage_data$V1 # null metric value for each assemblage
}
stopCluster(cl)

# Get the observed pattern of verticality metrics for each assemblage:
assemblage_data<-merge(assemblage_data, trait_data_Num, by='Binomial', all.x=TRUE, allow.cartesian=TRUE)
assemblage_data<-assemblage_data[, .(.N, geometric.mean(HabitatNum, na.rm=T),
                                     geometric.mean(TailProp, na.rm=T),
                                     geometric.mean(PCA1, na.rm=T)), by = .(Cell_Id55)] 
names(assemblage_data)<-c('Cell_Id55', 'Richness', 'ObsHabitat', 'ObsTailProp', 'ObsPCA1')

# Get the average value of each verticality metric from the null distribution:
assemblage_data$NullHabitat<-rowSums(NullValues_HabitatNum)/ncol(NullValues_HabitatNum)
assemblage_data$NullTailProp<-rowSums(NullValues_TailProp)/ncol(NullValues_TailProp)
assemblage_data$NullPCA1<-rowSums(NullValues_PCA1)/ncol(NullValues_PCA1)

# Compute the standard deviation of the null values produced for each assemblage:
assemblage_data$NullHabitatSD<-NA
assemblage_data$NullTailPropSD<-NA
assemblage_data$NullPCA1SD<-NA
for(i in 1:nrow(assemblage_data)){
  assemblage_data$NullHabitatSD[i]<-sd(as.data.frame(t(NullValues_HabitatNum))[,i])
  assemblage_data$NullTailPropSD[i]<-sd(as.data.frame(t(NullValues_TailProp))[,i])
  assemblage_data$NullPCA1SD[i]<-sd(as.data.frame(t(NullValues_PCA1))[,i])
} 

# Compute the standardized effect size for each verticality metrics as SES = (observed.value - mean.null.value) / SD.null.values:
assemblage_data$SESHabitat<-((assemblage_data$ObsHabitat - assemblage_data$NullHabitat)/assemblage_data$NullHabitatSD) # compute the z value
assemblage_data$SESTailProp<-((assemblage_data$ObsTailProp - assemblage_data$NullTailProp)/assemblage_data$NullTailPropSD) # compute the z value
assemblage_data$SESPCA1<-((assemblage_data$ObsPCA1 - assemblage_data$NullPCA1)/assemblage_data$NullPCA1SD) # compute the z value

# STEP 3 - EXAMINE THE RELATIONSHIP BETWEEN SNAKE RICHNESS AND VERTICALITY METRICS  
############################################################################################################################
# STEP 3 - EXAMINANOS LA RELACION ENTRE LA RIQUEZA DE ESPECIES Y NUESTRAS METRICAS

# Remove all objects from the workspace, but 'assemblage_data' and 'trait_data_Num':
rm(list=setdiff(ls(), c("assemblage_data", "trait_data_Num")))

# Reorder columns of the assemblage_data object:
names(assemblage_data)
assemblage_data<-assemblage_data[,c(1:5,12:14,6:11)]
assemblage_data<-as.data.frame(assemblage_data)

# Create an empty list to store the scatterplots:
MyPlot<-list()

# Define a vector holding the legend label for the plots:
axis_name<-c('Obs Vertical Habitat', 'Obs Tail Proportion', 'Obs PCA1 Score',
             'SES Vertical Habitat', 'SES Tail Proportion', 'SES PCA1 Score')

# Create the plots in for loop:
for(i in 1:6){
  
  # Select the metric for plotting:
  assemblage_data$Metric<-assemblage_data[,(i+2)]
  
  # Define the plot:
  MyPlot[[i]]<-qplot(data=assemblage_data, x=Richness, y=Metric, alpha=I(1/20)) + 
    geom_point(aes(x=Richness, y=Metric), size=3, alpha=0.5, shape=21, fill="gray50", colour="black") +
    stat_smooth(method="lm", alpha=0.5, colour="red", fill="red") +
    labs(x="", y=axis_name[i]) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=14),
          axis.line = element_line(colour="black"),
          axis.text = element_text(size=12, colour="black"),
          axis.text.x = element_text(angle=0),
          axis.title.x = element_text(margin=margin(t=1, r=0, b=0, l=0), colour="black"),
          axis.title.y = element_text(margin=margin(t=0, r=1, b=0, l=0) ,colour="black"),
          #plot.margin = rep(unit(0,"null"),2),
          panel.spacing = unit(0,"null"),
          plot.background=element_rect(fill="transparent", colour=NA),
          legend.position="none") +
    stat_cor(aes(y=Richness, x=Metric), size=4, method = "pearson", digits=3,
             label.y=(max(assemblage_data$Metric) - (0.025*(max(assemblage_data$Metric) - min(assemblage_data$Metric)))), # max minus 5% of the range
             label.x=(min(assemblage_data$Richness) + (0.025*(max(assemblage_data$Richness) - min(assemblage_data$Richness))))) # min plus 5% of the range
}

# Create a multipanel plot:
Multipanel_plot<-ggpubr::ggarrange(MyPlot[[1]], MyPlot[[2]], MyPlot[[3]],
                                   MyPlot[[4]] + xlab("Snake richness"),
                                   MyPlot[[5]] + xlab("Snake richness"),
                                   MyPlot[[6]] + xlab("Snake richness"),
                                   labels=c("A", "B", "C", "D", "E", "F"), align="hv",
                                   font.label=list(size=12, color = "black"), ncol=3, nrow=2); Multipanel_plot

ggsave("Verticality_vs_Richness_WeightedNullModel.png", plot=Multipanel_plot, width=12, height=7, units="in", bg = "transparent")
ggsave("Verticality_vs_Richness_WeightedNullModel.pdf", plot=Multipanel_plot, width=12, height=7, units="in", bg = "transparent")


######

# STEP 4 - PLOT THE GEOGRAPHICAL PATTERNS OF SNAKE VERTICALITY
############################################################################################################################
# STEP 4 - PLOTEAMOS EL PATRON GEOGRAFICO DE VERTICALIDAD

# Remove all objects from the workspace, but 'assemblage_data' and 'trait_data_Num':
rm(list=setdiff(ls(), c("assemblage_data", "trait_data_Num")))

# Load all shapefiles that will be used for part 1 of Figure 4:
setwd("")
selected_cells<-fread("", stringsAsFactors=TRUE, encoding="UTF-8") #grid used in step 2
terr_cover <- readOGR("") #Chaco boundary or region of interest
grid_cells <- readOGR("") #Grid size in shapefile format

# Set the coordinate reference systems that will be used:
wgs84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # WGS84
equalareaproj<-"+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs" # cylindrical equal area projection

# Verify the coordinate system of each shapefile and convert if necessary:
crs(terr_cover) # wgs84
crs(grid_cells) # equal area
grid_cells<-spTransform(grid_cells, equalareaproj)
terr_cover<-spTransform(terr_cover, equalareaproj)

# Identify those grid cells with at least 5 species (observed):
assemblage_data<-as.data.table(assemblage_data)
cells_to_keep<-assemblage_data[which(assemblage_data$Richness>=5),1] # all cells have more than 5 spp.

# Remove those assemblages with < 5 species:
assemblage_data<-assemblage_data[Cell_Id55 %in% as.character(cells_to_keep$Cell_Id55)]

# Combine assemblage-level estimates with the attribute table of grid_cells object:
grid_cells@data$n_order<-c(1:nrow(grid_cells@data))
grid_cells@data$Cell_Id55<-as.character(grid_cells@data$PageNumber)
assemblage_data$Cell_Id55<-as.character(assemblage_data$Cell_Id55)
grid_cells@data<-merge(grid_cells@data, assemblage_data, by="Cell_Id55", all.x=T)
grid_cells@data<-grid_cells@data[order(grid_cells@data$n_order),]

# Prepare the spatial objects for ggploting: 
grid_cells<-spTransform(grid_cells, equalareaproj)
grid_cells@data$id<-rownames(grid_cells@data)
grid_cell_df<-fortify(grid_cells, region="id")
grid_cell_df<-join(grid_cell_df, grid_cells@data, by="id")
terr_cover@data$id<-rownames(terr_cover@data)
terr_cover_df<-fortify(terr_cover, region="id")


# Create an empty list to store the scatterplots:
MyPlot<-list()

# Define a vector holding the legend label for the plots:
title_list<-c('Snake Richness', 'Obs Vert. Hab.', 'Obs TailProp', 'Obs PCA1 Score',
              'SES Vert. Hab.', 'SES TailProp', 'SES PCA1 Score')


# Create the plots in for loop:
for(i in 1:7){
  
  
  # Get the species richness of one major snake clade and remove assemblages without records:
  grid_cell_df_subset<-grid_cell_df  
  grid_cell_df_subset$Verticality<-grid_cell_df[,(i+16)]
  grid_cell_df_subset$n_order<-c(1:nrow(grid_cell_df_subset))
  grid_cell_df_subset<-droplevels(grid_cell_df_subset)
  grid_cell_df_subset<-grid_cell_df_subset[order(grid_cell_df_subset$n_order),]
  
  # Build the plot:
  MyPlot[[i]]<-ggplot(data=grid_cell_df_subset) +
    geom_polygon(data=terr_cover_df, aes(long, lat, group=group), colour=NA, fill="white", size=0.75) +   # plot landcover boundaries
    geom_polygon(data=grid_cell_df_subset, aes(long, lat, group=group, fill=Verticality)) + 
    geom_polygon(data=terr_cover_df, aes(long, lat, group=group), colour="white", fill=NA, size=0.75) + # replot landcover boundaries but without fill colour
    coord_equal() +
    theme(axis.line=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          legend.justification = c(1, 1),
          legend.position=c(0.05, 0.85), legend.direction="vertical",
          legend.text=element_text(size=rel(0.8), hjust=0),
          legend.background=element_rect(fill="white", size=0.5),
          legend.title=element_text(size=8),
          legend.title.align=0.5,
          legend.spacing.x = unit(0.2, 'cm'),
          legend.spacing.y = unit(0.2, 'cm'),
          panel.background=element_blank(),
          plot.background=element_blank(),
          plot.margin = rep(unit(-0.03,"null"),4),
          panel.spacing = unit(0,"null")) +
    guides(fill=guide_colorbar(title=title_list[i], label=T, nbin=15, barwidth=1, barheight=8, draw.ulim=T, draw.llim=T,
                               frame.colour="black", ticks.linewidth=0.5, ticks=T, ticks.colour="black")) +
    scale_fill_gradientn(colours=viridis_pal(direction=-1, option="magma")(10), na.value = NA)
}

# Create a multipanel plot:
Multipanel_plot<-ggpubr::ggarrange(MyPlot[[2]], MyPlot[[3]], MyPlot[[4]],
                                   MyPlot[[5]], MyPlot[[6]], MyPlot[[7]],
                                   labels=c("A", "B", "C", "D", "E", "F"), align="hv",
                                   font.label=list(size=12, color = "black"), ncol=3, nrow=2); Multipanel_plot

ggsave("GeographicalPattern.png", plot=Multipanel_plot, width=14, height=10, units="in", bg = "transparent")
ggsave("GeographicalPattern.pdf", plot=Multipanel_plot, width=14, height=10, units="in", bg = "transparent")

#####

# STEP 5 - COMPUTE THE SPATIAL CORRELATION BETWEEN SNAKE RICHNESS AND VERTICALITY METRICS
############################################################################################################################
# STEP 5 - CALCULAMOS LA CORRELACION ESPECIAL ENTRE LA RIQUEZA DE ESPECIES Y LAS METRICAS

# Remove all objects from the workspace, but 'grid_cells':
rm(list=setdiff(ls(), "grid_cells"))
# Separate the variables of interest:
assemblage_data<-as.data.frame(grid_cells@data)
assemblage_data<-assemblage_data[,c(1,6,7,8,10:17)]

head(assemblage_data)
head(grid_cells)

# Separate the geographical coordinates to perform correlations with spatially corrected degrees of freedom:
coords<-assemblage_data[,3:4]

# Test the correlation between different variables while accounting for spatial autocorrelation:
SpatialPack::modified.ttest(assemblage_data$Richness, assemblage_data$ObsHabitat, as.data.frame(coords), nclass=10)
SpatialPack::modified.ttest(assemblage_data$Richness, assemblage_data$ObsTailProp, as.data.frame(coords), nclass=10)
SpatialPack::modified.ttest(assemblage_data$Richness, assemblage_data$ObsPCA1, as.data.frame(coords), nclass=10)
SpatialPack::modified.ttest(assemblage_data$Richness, assemblage_data$SESHabitat, as.data.frame(coords), nclass=10)
SpatialPack::modified.ttest(assemblage_data$Richness, assemblage_data$SESTailProp, as.data.frame(coords), nclass=10)
SpatialPack::modified.ttest(assemblage_data$Richness, assemblage_data$SESPCA1, as.data.frame(coords), nclass=10)


# STEP 6 - PERFORM THE ANALYSIS OF DRIVERS OF VERTICAL STRATIFICATION
############################################################################################################################
# STEP 6 - ANALIZAMOS LAS VARIABLES QUE INFLUYEN EN LA ESTRATIFICACION VERTTICAL
#
library(readr)
library(MASS)
library(adespatial) # to compute minimum spanning tree
library(easyGgplot2) # to plot overlapping histograms
library(ggplot2) # to plot geographic variation in species description dates
library(MuMIn) # to perform model selection
library(pgirmess) # to compute spatial corrrelograms
library(raster) # to manipulate spatial data
library(rgdal) # to manipulate spatial data
library(rgeos) # to compute distance to nearest features
library(SoDA) # to convert geodesic coordinates into Euclidean ones
library(SpatialPack) # to compute spatially corrected correlations
library(spdep) # to compute SAR models
library(usdm) # to compute VIF

setwd("")
rm(list=ls())

data_snakes<- read.delim("") #load matrix with metrics and variables information
#####
site_descriptors <- data_snakes
##Perform the SAM for our metrics and variables
# Check for multicollinearity:
site_descriptors[c(12:20)]<-log10(site_descriptors[c(12:20)]+1)
predictors<-as.data.frame(scale(data_snakes[,c(12:21)], center=T, scale=T))
vif(predictors) #sequentially excluded variables with VIF values >10 
pearson_r<-as.data.frame(cor(predictors, use="complete.obs")) # pairwise Pearson correlation among preditors

# Get the geographic coordinates of species assemblages:
coords<-data_snakes[,2:3]
coords<-geoXY(latitude=coords$Lat, longitude=coords$Long, lat0=0, lon0=0, unit = 1)
coords<-coords/1000 # Convert the distance to km 

# Compute the threshold distace to build the neighbourhood connectivity matrix:
max_nbdist<-give.thresh(dist(coords)) # maximum distance from the minimum spanning tree

# Create the spatial weights matrix:
nbdist<-dnearneigh(x=coords, d1=0, d2=max_nbdist) # define connectivity matrix (0/1)
neigh.dist<-nbdists(nbdist, coords, longlat=F) # compute the Euclidean distance between neighbouring sites
inverse<-lapply(neigh.dist, function(x) (1/(x^2))) # compute the inverse distance weigthed matrix
nbdist.W<-nb2listw(neighbours=nbdist, glist=inverse, style="W", zero.policy=FALSE) # coding style W = row standardised

# Scale all continuous variables to get std. coefs.:
site_descriptors<-as.data.frame(scale(data_snakes[c(5:21)])) 

# Create model formulas with all possible combinations of predictor variables:
predictors<-site_descriptors[8:17] # subset of predictors being combined
model_combinations<-list()
for(i in 1:ncol(predictors)){ # create all predictor combinations
  model_combinations[[i]]<-combn(names(predictors), i, simplify=FALSE)
} 

# Merge all lists of combinations into one big list (1023 model combinations here):
model_combinations<-c(model_combinations[[1]], model_combinations[[2]], 
                      model_combinations[[3]], model_combinations[[4]],
                      model_combinations[[5]], model_combinations[[6]],
                      model_combinations[[7]], model_combinations[[8]],
                      model_combinations[[9]], model_combinations[[10]]) 

# Define the response variable, choose between, SES values, or Obs values:
response_variable<-site_descriptors$Richness
#response_variable<-site_descriptors$ObsTailProp
#response_variable<-site_descriptors$SESTailProp
#response_variable<-site_descriptors$ObsPCA1
#response_variable<-site_descriptors$SESPCA1
#response_variable<-site_descriptors$ObsHabitat
#response_variable<-site_descriptors$SESHabitat
# Create and store the model formulas for each predictor combinations:
formlist<-list() 
for(i in 1:length(model_combinations)){ # create model formulas with the combinations of explanatory variables
  formlist[[i]]<-as.formula(paste("response_variable~", paste((model_combinations)[[i]], collapse= "+")))
}
names(predictors)
# Create an empty dataframes to store the fitted values and summary output:
fitted_models<-as.data.frame(matrix(nrow=nrow(data_snakes), ncol=length(formlist))) 
set_of_models<-as.data.frame(matrix(nrow=length(formlist), ncol=(length(predictors)*2)+4)) # store the coefficients and their standard errors in each SAR model 
colnames(set_of_models)<-(c("bio_1", "bio_15","Evenness", "AET", "FRAG","Homogeneity", "NPP", "SAND", "TC","ndvi", 
                            "bio_1_SE","bio_15_SE", "Evenness_SE", "AET_SE", "FRAG_SE","Homogeneity_SE", "NPP_SE", "SAND_SE", "TC_SE","ndvi_SE",
                            "Pseudo_R2", "AICc", "delta_AICc", "wAICc"))


# Extract the pseudo-R², AICc, and standardized coefficents for each SAR model:
for(i in 1:length(formlist)){
  model_tested<-errorsarlm(lm(formlist[[i]], data=predictors), listw=nbdist.W) # build a given SAR model
  fitted_models[,i]<-model_tested[[21]] # store the fitted values from the SAR model
  set_of_models[i,21]<-(cor(response_variable, model_tested[[21]])^2)[[1]] # compute the pseudo-R²
  set_of_models[i,22]<-AICc(model_tested) # compute the AICc
  model_coefs<-as.data.frame(coef(model_tested)) # compute the standardized coefficients
  std_errors<-as.data.frame(model_tested[[5]]) # coefficient standard errors
  row.names(std_errors)<-row.names(as.data.frame(model_tested[[4]])) # rename rows for latter use
  
  # Store the standardized coefficients for each variable in each SAR model:
  if (length(model_coefs[which(row.names(model_coefs)=="bio_1"),1])==1){set_of_models[i,1]<-model_coefs[(row.names(model_coefs)=="bio_1"),1]}
  if (length(model_coefs[which(row.names(model_coefs)=="bio_15"),1])==1){set_of_models[i,2]<-model_coefs[(row.names(model_coefs)=="bio_15"),1]}
  if (length(model_coefs[which(row.names(model_coefs)=="Evenness"),1])==1){set_of_models[i,3]<-model_coefs[(row.names(model_coefs)=="Evenness"),1]}
  if (length(model_coefs[which(row.names(model_coefs)=="AET"),1])==1){set_of_models[i,4]<-model_coefs[(row.names(model_coefs)=="AET"),1]}
  if (length(model_coefs[which(row.names(model_coefs)=="FRAG"),1])==1){set_of_models[i,5]<-model_coefs[(row.names(model_coefs)=="FRAG"),1]}
  if (length(model_coefs[which(row.names(model_coefs)=="Homogeneity"),1])==1){set_of_models[i,6]<-model_coefs[(row.names(model_coefs)=="Homogeneity"),1]}
  if (length(model_coefs[which(row.names(model_coefs)=="NPP"),1])==1){set_of_models[i,7]<-model_coefs[(row.names(model_coefs)=="NPP"),1]}
  if (length(model_coefs[which(row.names(model_coefs)=="SAND"),1])==1){set_of_models[i,8]<-model_coefs[(row.names(model_coefs)=="SAND"),1]}
  if (length(model_coefs[which(row.names(model_coefs)=="TC"),1])==1){set_of_models[i,9]<-model_coefs[(row.names(model_coefs)=="TC"),1]}
  if (length(model_coefs[which(row.names(model_coefs)=="ndvi"),1])==1){set_of_models[i,10]<-model_coefs[(row.names(model_coefs)=="ndvi"),1]}
  
  # Store the standard error for each variable in each SAR model:
  if (length(std_errors[which(row.names(std_errors)=="bio_1"),1])==1){set_of_models[i,11]<-std_errors[(row.names(std_errors)=="bio_1"),1]}
  if (length(std_errors[which(row.names(std_errors)=="bio_15"),1])==1){set_of_models[i,12]<-std_errors[(row.names(std_errors)=="bio_15"),1]}
  if (length(std_errors[which(row.names(std_errors)=="Evenness"),1])==1){set_of_models[i,13]<-std_errors[(row.names(std_errors)=="Evenness"),1]}
  if (length(std_errors[which(row.names(std_errors)=="AET"),1])==1){set_of_models[i,14]<-std_errors[(row.names(std_errors)=="AET"),1]}
  if (length(std_errors[which(row.names(std_errors)=="FRAG"),1])==1){set_of_models[i,15]<-std_errors[(row.names(std_errors)=="FRAG"),1]}
  if (length(std_errors[which(row.names(std_errors)=="Homogeneity"),1])==1){set_of_models[i,16]<-std_errors[(row.names(std_errors)=="Homogeneity"),1]}
  if (length(std_errors[which(row.names(std_errors)=="NPP"),1])==1){set_of_models[i,17]<-std_errors[(row.names(std_errors)=="NPP"),1]}
  if (length(std_errors[which(row.names(std_errors)=="SAND"),1])==1){set_of_models[i,18]<-std_errors[(row.names(std_errors)=="SAND"),1]}
  if (length(std_errors[which(row.names(std_errors)=="TC"),1])==1){set_of_models[i,19]<-std_errors[(row.names(std_errors)=="TC"),1]}
  if (length(std_errors[which(row.names(std_errors)=="ndvi"),1])==1){set_of_models[i,20]<-std_errors[(row.names(std_errors)=="ndvi"),1]}
}  
set_of_models$delta_AICc<-set_of_models$AICc - min(set_of_models$AICc) # compute the delta AICc
set_of_models$wAICc<-Weights(set_of_models$AICc) # compute the wAICc

# Visualize the output table
set_of_models<-set_of_models[order(set_of_models$wAICc, decreasing=T), ] # order models according to wAICc 
row.names(set_of_models)<-1:nrow(set_of_models)
View(set_of_models[,c(1:10, 21,23)])

# Get the averaged model coefficients based on wAICc:
avg_model<-as.data.frame(matrix(nrow=5, ncol=10)) 
colnames(avg_model)<-(c("bio_1","bio_15", "Evenness", "AET", "FRAG","Homogeneity", "NPP", "SAND", "TC", "ndvi"))
row.names(avg_model)<-c("Avg.coef","Std.error","Rel.Import", "Lower_IC", "Upper_IC")

for(i in 1:10){ # averaged coefficients
  candidate_models<-set_of_models[which(set_of_models[,i]!=0),]
  avg_model[3,i]<-sum(candidate_models$wAICc) # compute the relative importance
  candidate_models$wAICc<-Weights(candidate_models$AICc) # rescale the wAICc
  avg_model[1,i]<-sum(candidate_models[,i]*candidate_models$wAICc) # standardized coefficients of the weighted model
}

for(i in 11:20){ # unconditional variance estimator for the averaged model coefficients
  candidate_models<-set_of_models[which(set_of_models[,i]!=0),] # select all models for a given predictor
  weighted_var<-as.data.frame(matrix(nrow=nrow(candidate_models), ncol=1)) 
  
  for(j in 1:nrow(candidate_models)){ # j models, varying from 1 to 512 
    delta_coef<-(candidate_models[j,(i-10)] - avg_model[1,(i-10)])^2 # error due to model selection uncertainty
    var_error<-(candidate_models[j,i])^2 # error in parameter estimation of the model 'j'
    weighted_var[j,1]<-(sqrt(var_error + delta_coef))*candidate_models[j,24]}
  
  avg_model[2,(i-10)]<-sum(weighted_var[,1]) # unconditional standard error
} 

# Get the unconditional confidence intervals:
for(i in 1:10){
  avg_model[4,i]<-avg_model[1,i]-(1.96*avg_model[2,i]) # lower confidence interval
  avg_model[5,i]<-avg_model[1,i]+(1.96*avg_model[2,i]) # upper confidence interval
}
View(t(avg_model))
#transpose and export avg model results to plot
avg_modelRichness <-as.data.frame(t(as.matrix(avg_model)))
library(dplyr)
sarrich_avg_model <- tibble::rownames_to_column(avg_modelRichness, "X")
#changes the names for the predictor variables
write.csv(sarrich_avg_model,"", row.names = FALSE)

# Get the weighted fitted values and the residuals from the average model:
weighted_fitted_values<-as.data.frame(matrix(nrow=nrow(fitted_models), ncol=ncol(fitted_models))) 
for(i in 1:ncol(fitted_models)){
  weighted_fitted_values[,i]<-(set_of_models[i,24]*fitted_models[,i])
}
avg_fitted_values<-rowSums(weighted_fitted_values)
avg_residuals<-response_variable - avg_fitted_values

shapiro.test(avg_residuals)
hist(avg_residuals)

# Get the pseudo-R² of the average model
(cor.test(response_variable,avg_fitted_values)[[4]])^2

##########################################################################################################################
# STEP 7 - Plot the Moran's I correlogram of the average model residuals.

pgi.cor<-as.data.frame(correlog(coords=coords, z=avg_residuals, method="Moran"))
pgi.cor2<-as.data.frame(correlog(coords=coords, z=response_variable, method="Moran"))
plot(coef~dist.class, xlab="Geographic distance", ylab="Moran's I", pch=1, bty="l", cex.lab=1.5, ylim=c(-1,1), data=pgi.cor, cex=1.5)
points(pgi.cor2$dist.class, pgi.cor2$coef, pch=16, cex=1.5) # plot the Moran's I correlogram for the raw values of metrics
abline (h=0, lty=2, col="black", lwd=3)
lines(pgi.cor$dist.class, pgi.cor$coef, type = "l", col="blue", lwd =3) # add connected line segments to Morans'I of dbMEM model residuals
lines(pgi.cor2$dist.class, pgi.cor2$coef, type = "b", col="red", lwd=3) # add connected line segments to Morans'I of the raw value of body size
legend((pgi.cor[14,1])/3.2, 1, c("Richness", "SAR residuals"), pch=c(16,1), lty=c(2,1), bty="n", cex=1,col=c("red", "blue"))

####SENSITIVITY ANALYSIS####
#same analysis, after removing genus of mainly arboreal species
# STEP 7 - PHYLOGENETIC METRICS
##########################################################################################################################
# STEP 7 - METRICAS FILOGENETICAS

library(phyloregion)
library(ape)
library(Matrix)
library(sp)
library(raster)
library(cowplot)
library(rgdal)
library(ggplot2)
library(letsR)
library(picante)
rm(list=ls()); gc()
setwd("")
grid_cells <- readOGR("") #grid in shapefile format

#changing the names to grids
names(grid_cells)
colnames(grid_cells@data)[3] = "grids"
###Add tree
setwd("")
dir()
tree <- read.tree("") #Phylogenetic tree
#add coomunity matrix and convert to sparse community matrix
community <- read.delim("", sep=";", row.names = 1) #presence/absence matrix, species in column
sparse_comm <- dense2sparse(community)
object.size(sparse_comm)

#Weighted endemism
Endemismo<- weighted_endemism(sparse_comm)
m <- merge(grid_cells, data.frame(grids=names(Endemismo), WE=Endemismo), by="grids")
m <- m[!is.na(m@data$WE),]
plot_swatch(m, values = m$WE, leg = 3, border = NA, pos="topleft",
            col = hcl.colors(n=20, palette = "Blue-Red 3", rev=FALSE))

#Phylogenetic diversity

phylo_diversity <- PD(sparse_comm, tree)
mpd <- merge(grid_cells, data.frame(grids=names(phylo_diversity), pd=phylo_diversity), by="grids")
mpd <- mpd[!is.na(mpd@data$pd),]

plot_swatch(mpd, values = mpd$pd, border=NA, leg=4, pos = "topleft",
            col = hcl.colors(n=20, palette = "Blue-Red 3", rev=FALSE))
#You can export the results, if you want to work with GIS programs.

###Phylo endemism
phy_endemism <- phylo_endemism(sparse_comm, tree)
mpe <- merge(grid_cells, data.frame(grids=names(phy_endemism), pe=phy_endemism), by="grids")
mpe <- mpe[!is.na(mpe@data$pe),]
plot_swatch(mpe, values = mpe$pe, border=NA, leg=3, pos="topleft",
            col = hcl.colors(n=20, palette = "Blue-Red 3", rev=FALSE))
#Phyloregion

#phylo beta diversity
pb <- phylobeta(sparse_comm, tree)

##To identify the optimal number of Phyloregions
bd<- beta_diss(sparse_comm)
bdb<- optimal_phyloregion(bd[[1]])
plot(bdb$df$k, bdb$df$ev, ylab = "Explained variances",
     xlab = "Number of clusters", cex.lab=1.2)
lines(bdb$df$k[order(bdb$df$k)], bdb$df$ev[order(bdb$df$k)], pch = 1)
points(bdb$optimal$k, bdb$optimal$ev, pch = 21, bg = "red", cex = 2, lwd=1)
points(bdb$optimal$k, bdb$optimal$ev, pch = 21, bg = "red", type = "h",lwd=4 ,col="red")

#Phyloregion
phy_region <- phyloregion(pb[[1]], shp=grid_cells, k=16)

plot(phy_region, palette="NMDS")
plot_NMDS(y, cex=5, main="stress value = 0.07", cex.main=0.8)
text_NMDS(y, col="white", cex=2)

##You can export the results, if you want to work with GIS programs.
##Phylogenetic species variability
comunity<- read.delim("", sep=";", row.names = 1) #presence/absence matrix, species in column

Macth<- match.phylo.comm(tree, comunity)
PSVsp<- psv(comunity, tree, compute.var = T)


