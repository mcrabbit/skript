if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("progeny")
install.packages("pheatmap")


## To install the new version until it is submitted to Bioconductor use:
devtools::install_github("saezlab/progeny", force = TRUE)

library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)

## We create a data frame with the specification of the cells that belong to 
## each cluster to match with the Progeny scores.
Idents(experiment_1_72h)<-"orig.ident"
CellsClusters <- data.frame(Cell = names(Idents(experiment_1_72h)),
                            CellType = as.character(Idents(experiment_1_72h)),
                            stringsAsFactors = FALSE)

## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 

experiment_1_72h <- progeny(experiment_1_72h, scale=FALSE, organism="Mouse", top=500, perm=1, return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
experiment_1_72h <- Seurat::ScaleData(experiment_1_72h, assay = "progeny")

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <-
  as.data.frame(t(GetAssayData(experiment_1_72h, slot = "scale.data",
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white", "red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0,
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength,
                      max(summarized_progeny_scores_df),
                      length.out=floor(paletteLength/2)))

progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14,
                        fontsize_row = 10, cluster_rows=T, cluster_cols = F,
                        color=myColor, breaks = progenyBreaks,
                        main = "PROGENy (conditions)", angle_col = 45,
                        treeheight_col = 0,  border_color = "black")


###PROGENY ON THE CLUSTERS####

## We create a data frame with the specification of the cells that belong to 
## each cluster to match with the Progeny scores.
Idents(experiment_1_72h) <- "seurat_clusters"
CellsClusters <- data.frame(Cell = names(Idents(experiment_1_72h)),
                            CellType = as.character(Idents(experiment_1_72h)),
                            stringsAsFactors = FALSE)

## We compute the Progeny activity scores and add them to our Seurat object
## as a new assay called Progeny. 

experiment_1_72h <- progeny(experiment_1_72h, scale=FALSE, organism="Mouse", top=500, perm=1, return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
experiment_1_72h <- Seurat::ScaleData(experiment_1_72h, assay = "progeny")

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <-
  as.data.frame(t(GetAssayData(experiment_1_72h, slot = "scale.data",
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white", "red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0,
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength,
                      max(summarized_progeny_scores_df),
                      length.out=floor(paletteLength/2)))

progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14,
                        fontsize_row = 10, cluster_rows=T, cluster_cols = T,
                        color=myColor, breaks = progenyBreaks,
                        main = "PROGENy (clusters)", angle_col = 45,
                        treeheight_col = 50,  border_color = "black")

