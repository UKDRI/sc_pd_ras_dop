#!/usr/bin/env Rscript
## [[file:../org/initial_exploration_and_extraction.org::*Data info][Data info:1]]
library(dplyr)
library(Seurat)
library(patchwork)

data_location <- "/Users/samneaves/Documents/test_data/forSam/SNatlas_DaNs_seurat.RData"
load(data_location)
## Data info:1 ends here

## [[file:../org/initial_exploration_and_extraction.org::*Data info][Data info:2]]
head(sn_atlas_dans$CellSubType)
## Data info:2 ends here

## [[file:../org/initial_exploration_and_extraction.org::*Get the normalised expression values][Get the normalised expression values:1]]
genes_norm<-as.matrix(GetAssayData(sn_atlas_dans, slot = "data"))
## Get the normalised expression values:1 ends here

## [[file:../org/initial_exploration_and_extraction.org::*Get the count expression values][Get the count expression values:1]]
genes_counts<-as.matrix(GetAssayData(sn_atlas_dans, slot = "counts"))
## Get the count expression values:1 ends here

## [[file:../org/initial_exploration_and_extraction.org::*Finding the overlap][Finding the overlap:1]]
network_data <- read.table("../../data/pybravo_output/expanded_reg_md10-unified.sif", quote = "", header = FALSE, stringsAsFactors = FALSE, sep = "\t")
head(network_data)
nodes <- unique(c(network_data$V1, network_data$V3))
## Finding the overlap:1 ends here

## [[file:../org/initial_exploration_and_extraction.org::*Finding the overlap][Finding the overlap:2]]
length(nodes)
## Finding the overlap:2 ends here

## [[file:../org/initial_exploration_and_extraction.org::*Finding the overlap][Finding the overlap:3]]
# Convert the rownames of genes_counts into a vector
genes_names <- rownames(genes_counts)

# Get the intersection between the 'name' column in node_table and the rownames of genes_counts_subset_df
common_names <- intersect(nodes, genes_names)

# View the result
length(common_names)
## Finding the overlap:3 ends here

## [[file:../org/initial_exploration_and_extraction.org::*Complexes][Complexes:1]]
# Identify rows in the 'name' column that contain '/'
complexes <- grep("/", nodes)

# Count the number of complexes
num_complexes <- length(complexes)

# View the count
num_complexes
## Complexes:1 ends here

## [[file:../org/initial_exploration_and_extraction.org::*Edge types][Edge types:1]]
# Count occurrences of each edge type in the second column (V2)
edge_counts <- table(network_data$V2)

# Display the counts of each edge type
print(edge_counts)
## Edge types:1 ends here

## [[file:../org/initial_exploration_and_extraction.org::*Edge types][Edge types:2]]
str(sn_atlas_dans)
summary(sn_atlas_dans)
## Edge types:2 ends here

## [[file:../org/initial_exploration_and_extraction.org::*Select cells and genes][Select cells and genes:1]]
sn_atlas_dans_ctrl_vs_pd <- subset(sn_atlas_dans, subset = Disease %in% c("CTR", "PD_B5-6", "PD_B3-4"))
sn_atlas_dans_ctrl_vs_pd_pathway <- sn_atlas_dans_ctrl_vs_pd[common_names,]
dim(sn_atlas_dans_ctrl_vs_pd_pathway)
## Select cells and genes:1 ends here

## [[file:../org/initial_exploration_and_extraction.org::*Select cells and genes][Select cells and genes:2]]
# Extract the expression matrix for the subset
expression_matrix <- as.data.frame(t(as.matrix(sn_atlas_dans_ctrl_vs_pd_pathway@assays$RNA@data)))

# Add Disease column from metadata
expression_matrix$Disease <- sn_atlas_dans_ctrl_vs_pd_pathway@meta.data$Disease

# Save to CSV
write.csv(expression_matrix, file = "../../data/cell_expression_with_disease.csv", row.names = TRUE)
## Select cells and genes:2 ends here

## [[file:../org/initial_exploration_and_extraction.org::*Select cells and genes][Select cells and genes:3]]
# Extract the count matrix for the subset and transpose it for cells as rows
count_matrix <- as.data.frame(t(as.matrix(sn_atlas_dans_ctrl_vs_pd_pathway@assays$RNA@counts)))

# Add Disease column from metadata
count_matrix$Disease <- sn_atlas_dans_ctrl_vs_pd_pathway@meta.data$Disease

# Save to CSV
write.csv(count_matrix, file = "../../data/cell_count_with_disease.csv", row.names = TRUE)
## Select cells and genes:3 ends here
