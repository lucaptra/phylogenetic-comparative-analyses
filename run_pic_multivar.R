# Purpose: Compute phylogenetic independent contrasts (PICs) of multivariate traits
# Usage: Rscript run_pic_multivar.R [working directory] [name of file containing trait data] [phylogeny file in newick format]

# Preliminary steps to run in the Terminal
# # Compile pic_multivar.c
# R64 CMD SHLIB pic_multivar.c # 32 or 64-bit depending on your architecture

# Load libraries
library(ape)
library(caper)

# Set your working directory (which should contain the pic_multivar.R script, your data files, etc.)
setwd(args[1])

# Import the function to calculate the multivariate PICs
source("pic_multivar.R")

# Pass arguments from the command line
args = commandArgs(trailingOnly = TRUE)

# Import your multivariate trait data (e.g., phylogenetic principal components)
# Warnings:
#	The data file should contain informative column names.
#	Column 1 must contain the names of species.
#	Column 2 must contain the taxonomic grouping information (e.g., families, orders).
df <- read.table(file = args[2], header = TRUE, sep = '\t', quote = '', stringsAsFactors = FALSE)

# Create group factors
groups <- factor(df[, 2])

# Import your phylogenetic tree
tree <- read.tree(args[3])

# For each group of interest
for (g in levels(groups)) {
	
	# Extract trait values, excluding column 2 which contains the grouping information
	group.df <- df[df[, 2] == g, c(1, 3:dim(df)[2])]
	
	# Set species as the row names
	row.names(group.df) <- df[df[, 2] == g, 1]
	
	# Create a comparative data set combining your phylogeny with the multivariate trait data
	# Notes:
	#	Data for species not found in the phylogeny are removed.
	#	Tips in the phylogeny with no corresponding trait data are also dropped.
	pics <- comparative.data(tree, group.df, names(group.df)[1], vcv = TRUE, vcv.dim = 3, na.omit = TRUE, warn.dropped = TRUE)
	
	# Create a matrix of the trait data from the comparative data set
	mat <- as.matrix(sapply(pics$data, as.numeric))
	row.names(mat) <- row.names(pics$data)
	
	# Extract the phylogeny from the comparative data set
	phy <- pics$phy
	
	# Calculate the PICs
	pic_multivar(mat, phy, scaled = TRUE, var.contrasts = FALSE, rescaled.tree = FALSE)
	
}
