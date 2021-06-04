library(devtools)
library(usethis)
BiocManager::install("multtest")
library(multtest)
library(metap)
library(mutoss)
library(mvtnorm)
remotes::install_version("SDMTools", "1.1-221")
library(SDMTools)
remotes::install_version("Seurat", version = "3.1.0")
library(Seurat)
library(BiocManager)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(RColorBrewer) 
library(ggraph)
library(clustree)
#install harmony
devtools::install_github("immunogenomics/harmony")
library(Rcpp)
library(harmony)
library(viridisLite)
library(viridis)
library(sctransform)
library(wesanderson)
#install fgsea
BiocManager::install("fgsea")
library(fgsea)
library(msigdbr)
library(ggrepel)
library(cowplot)
#install seurat-wrappers
#register on github
usethis::create_github_token()
install.packages("gitcreds")
#check out git config file
edit_git_config()
#type credentials and connect to github
library(gitcreds)
gitcreds_set()
#connect project to github
library(usethis)
use_git()
use_github()

devtools::install_github('satijalab/seurat-wrappers')
library(remotes)
library(SeuratWrappers)#
library(conos)
install.packages("devtools", type = "win.binary")
devtools::install_github('satijalab/SeuratData')
gh_install_packages('wmacnair/psupertime')
install.packages(SeuratData)
library(SeuratData)#
install.packages("githubinstall")
library(githubinstall)
install_github('wmacnair/psupertime', build_vignettes=TRUE)
library('psupertime')#
gh_install_packages("psupertime")
#install singleCellExperiment
BiocManager::install("SingleCellExperiment")
library('SingleCellExperiment')
library(magrittr)
library(data.table)
library(stringr)
library(Matrix)
library(tidyr)

# import function
data_to_sparse_matrix <- function(data.st_file_path) {
  # read in file with cell index - gene name - values
  # import for one cartridge, one sample
  input <-read.table(data.st_file_path, header = T)
  # transform to matrix (data.frame actually)
  # we take as default values from the column "RSEC_Adjusted_Molecules" (= error corrected UMIs)
  mat <- input %>% pivot_wider(id_cols = Gene, 
                               values_from = RSEC_Adjusted_Molecules, 
                               names_from = Cell_Index, values_fill = 0)  %>% 
    tibble::column_to_rownames("Gene")
  # convert to sparse matrix (~ dgMatrix)
  sparse_mat = Matrix(as.matrix(mat),sparse=TRUE)
  return(sparse_mat)
}


# the dgMatrix is a valid input to create the Seurat object
crypto_infected <- data_to_sparse_matrix("/media/Coco/Collaborations/Crypto/Data analysis/Expression_data/crypto_infected.st")
mock_infected <- data_to_sparse_matrix("/media/Coco/Collaborations/Crypto/Data analysis/Expression_data/mock_infected.st")

crypto <- CreateSeuratObject(crypto_infected,project = "crypto_infected")
mock <- CreateSeuratObject(mock_infected, project = "mock_infected")


####MERGE ALL SEURAT OBJECT####
experiment_1_72h <- merge(crypto,mock, add.cell.ids = c("crypto", "mock"))
length(experiment_1_72h@active.ident)

RidgePlot(experiment_1_72h, features="nFeature_RNA", group.by = "orig.ident")

#check mito
experiment_1_72h[["percent.mt"]] <- PercentageFeatureSet(experiment_1_72h, pattern = "^mt.")
VlnPlot(experiment_1_72h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
RidgePlot(experiment_1, features="percent.mt", group.by = "orig.ident")


#####NORMALIZATION AND VARIABLE GENES####
experiment_1_72h <- subset(experiment_1_72h, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
experiment_1_72h <- NormalizeData(experiment_1_72h, normalization.method = "LogNormalize", scale.factor = 10000)
experiment_1_72h <- FindVariableFeatures(experiment_1_72h, selection.method = "vst", nfeatures = 2000)
top50 <- head(VariableFeatures(experiment_1_72h), 50)
plot1 <- VariableFeaturePlot(experiment_1_72h)
plot2 <- LabelPoints(plot = plot1, points = top50, repel = TRUE)
plot2

####SCALING AND DIMENSIONALITY REDUCTION###
all.genes <- rownames(experiment_1_72h)
experiment_1_72h <- ScaleData(experiment_1_72h, features = all.genes, vars.to.regress = c("nFeature_RNA", "percent.mt")) #should percent.crypto also be regressed out?
experiment_1_72h <- RunPCA(experiment_1_72h, features = VariableFeatures(object = experiment_1_72h))
DimPlot(experiment_1_72h, reduction = "pca", cols=col_vector, label = T)
VizDimLoadings(experiment_1_72h, dims = 3:4, reduction = "pca")
ElbowPlot(experiment_1_72h)
experiment_1_72h <- FindNeighbors(experiment_1_72h, dims = 1:50)
experiment_1_72h <- FindClusters(experiment_1_72h, resolution = 0.4)
experiment_1_72h <- RunUMAP(experiment_1_72h, dims = 1:50)
DimPlot(experiment_1_72h, order=T, group.by = "orig.ident", pt.size = 0.1, label=T, cols = col_vector)
DimPlot(experiment_1_72h, order=T, group.by = "seurat_clusters", pt.size = 0.1, label=T, cols = col_vector)

experiment_1_72h[["percent.crypto"]] <- PercentageFeatureSet(experiment_1_72h, pattern = "^CPATCC-")
FeaturePlot(experiment_1_72h, feature="percent.crypto", order=T, pt.size = 0.1, label = T)
plot <- VlnPlot(experiment_1_72h, feature="percent.crypto", pt.size = F)
plot + 
  theme(axis.text.x = element_text(angle = 45, face="bold", hjust=1), axis.text.y = element_text(face="italic")) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ theme(legend.position="bottom")+
  coord_flip()+ theme(legend.position = 'none')
FeaturePlot(experiment_1_72h, label = T, feature = "percent.crypto")
FeaturePlot(experiment_1_72h, label = T, feature = "Ptprc")
FeaturePlot(experiment_1_72h, label = T, features = c("Lgr5", "Olfm4", "Mki67"))
FeaturePlot(experiment_1_72h, label = T, features = c("Muc2"))


#CLUSTER MARKERS AND ANNOTATION#######
experiment_1_72h_markers <- FindAllMarkers(object = experiment_1_72h, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(experiment_1_72h_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))

Infected_markers <- FindAllMarkers(object = Infected, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
View(Infected_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC))

markers_5_vs_1 <-FindMarkers(Infected, ident.1 = 1, ident.2 = 5, min.pct = 0.25, logfc.threshold = 0.25)
View(markers_5_vs_1 %>% top_n(n = 300, wt = avg_logFC))

#general immune markers
markers.to.plot <-  c("Col1a1", "Vim", "Acta2",
                      "C1qb", "H2-Aa", "Cd74", "Ptprc", #4 = immune cells 
                      "Chga", "Chgb", "Cd24a", #9 = tuft cell markers
                      "Gpx1", "Car4",  #6 = immature 
                      "Gpx2", "Lgals4", "Lypd8", "Minos1", "Reg3b", "Reg3g", "Reg3a", "Nlrp6", "Ccl25",
                      "Apoa4", "Apoc3", "Apoa1", "Npc1l1", "Lgals3", "Egfr", "Klf4", "Junb","Mep1a", "Cre3l3", "Ephx2",
                      "Olfm4", "Lgr5", "Mki67", "Ccnd1", #stem cells /TA$
                      "Mptx2", "Ang4","Lyz1", #5 = Paneth cells  
                      "Atoh1", "Muc2", "Tff3", "percent.crypto") #5=goblet cells)

plot <- DotPlot(experiment_1_72h, features = markers.to.plot)
plot + 
  theme(axis.text.x = element_text(angle = 45, face="bold", hjust=1), axis.text.y = element_text(face="italic")) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))+ theme(legend.position="right")

bottom.landmark.genes <- c("Fabp1", "Gpx2", "Gsta1", "Lgals4", "Lypd8", "Minos1", "Reg3b", "Reg3g", "Reg3a", "Nlrp6", "Il18", "Ccl25")
DotPlot(Infected, features = bottom.landmark.genes) #> cluster 2,3 

top.landmark.genes <- c("Ada", "Apoa4", "Apoc3", "Apoa1", "Npc1l1", "Lgals3", "Egfr", "Klf4", "Junb")
DotPlot(Infected, features = top.landmark.genes)

V1.2.landmark.genes <- c("Slc7a9", "Slc7a8", "Slc7a7")
DotPlot(Infected, features = V1.2.landmark.genes) 

V2.3landmark.genes <- c("Slc5a1", "Slc2a5", "Slc2a1")
DotPlot(Infected, features = V2.3landmark.genes) 

V4.5.landmark.genes <- c("Slc15a1")
DotPlot(Infected, features = V4.5.landmark.genes) 


#cluster annotation
Infected <- RenameIdents(Infected, 
                                 `0` = "crypto-infected enterocytes 1", 
                                 `1` = "enterocytes (mature, villus tip)", 
                                 `2` = "goblet cells", 
                                 `3` = "SC and early precursors", 
                                 `4` = "crypto-infected enterocytes 2", 
                                 `5` = "crypto-infected enterocytes 3", 
                                 `6` = "enteroendocrine cells", 
                                 `7` = "immune cells", 
                                 `8` = "crypto-infected enterocytes 4", 
                                 `9` = "mesenchymal cells")

#cluster 0 > mature enterocytes > Ace, Ace2, apob 
#cluster 1 > crypto infected enterocytes (they still express Ace and other enterocyte markers) > top 50 markers are crypto genes
#cluster 2 > goblet cells > Muc2, Tff3
#cluster 3 > early precursor and SC, high in ribosomal proteins > Ly6a, Krt19, Cd74  
#cluster 4 > crypto infected immature enterocytes (Gpx1, Car4) 
#cluster 5 > crypto infected enterocytes
#cluster 7 > Gzm, conplement, MHCII > immune cells
#cluster 8 > crypto infected enterocytes (they still express Ace and other enterocyte markers) > top 50 markers are crypto genes
#cluster 9 > mesenchymal cells or monocytes > Acta2, collagens, Vim
#cluster 10 > crypto infected enterocytes (they still express Ace and other enterocyte markers) > top 50 markers are crypto genes
#cluster 11 > Enteroendocrine
#cluster 12 > Enteroendocrine & Tuft (Cd24a)
#cluster 13 > crypto infected enterocytes (they still express Ace and other enterocyte markers) > top 50 markers are crypto genes

DimPlot(eosinophils_mito, label = TRUE, pt.size = 1)
