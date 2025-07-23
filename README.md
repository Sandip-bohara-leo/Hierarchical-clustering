# Hierarchical-clustering
This repository presents various hierarchical clustering methods, with dendrograms visualized in multiple styles, 
including rectangular, horizontal, circular, and phylogenetic layouts. It also includes calculation of cluster means and inter-cluster distances, 
with enhanced visual customization using the dendextend and circlize libraries.

library(readxl)
library(dendextend)
library(circlize)
library(dplyr)

# Load data
data <- WHEAT

# Convert to data frame and set row names
data <- as.data.frame(data)
rownames(data) <- data$`Accession name`

# Select clustering columns
data_for_clustering <- data[, c("DTH","DTF","DOM","PH","SPL","NOT","SPPS","SPSP","SPS","SL","THSW","TY"
)]

# Scale data
scaled_data <- scale(data_for_clustering)

# Compute distance matrix
d <- dist(scaled_data, method = "euclidean")

# Perform hierarchical clustering
hc <- as.dendrogram(hclust(d, method = "average"))

##Cluster dendrogram
library(factoextra)

fviz_dend(hc,cex=0.3,lwd = 0.7,k=3,
          k_colors = c("red","green","magenta"))
fviz_dend(hc,cex=0.5,lwd = 0.7,k=3,
          k_colors = "jco")
fviz_dend(hc,cex=0.5,lwd = 0.7,k=3,
          rect = TRUE,
          rect_border = "gray",
          rect_fill = FALSE)
fviz_dend(hc,cex=0.5,lwd = 0.7,k=3,
          rect = TRUE,
          k_colors = "jco",
          rect_border="jco",
          rect_fill=TRUE,
          ggtheme=theme_bw())
fviz_dend(hc,cex=0.5,lwd = 0.7,k=3,
          rect = TRUE,
          k_colors = c("red","green3","magenta"),
          rect_border= c("red","green3","magenta"),
          rect_fill=TRUE)
fviz_dend(hc,cex=0.5,lwd = 0.7,k=3,
          rect = TRUE,
          k_colors = c("red","green3","magenta"),
          rect_border= c("red","green3","magenta"),
          rect_fill=TRUE,
          horiz=TRUE)
fviz_dend(hc,cex=0.5,lwd = 0.7,k=3,
          rect = TRUE,
          k_colors = "jco",
          rect_border="jco",
          rect_fill=TRUE,
          type = "circular")
fviz_dend(hc,cex=0.8,lwd = 0.7,k=3,
          rect = TRUE,
          k_colors = "jco",
          rect_border="jco",
          rect_fill=TRUE,
          type = "phylogenic",
          repel = TRUE)
#layout for phylogenetci dendrogram
fviz_dend(hc,cex=0.8,lwd = 0.7,k=3,
          rect = TRUE,
          k_colors = "jco",
          rect_border="jco",
          rect_fill=TRUE,
          type = "phylogenic",
          repel = TRUE,
          phylo_layout = "layout.gem")
fviz_dend(hc,cex=0.8,lwd = 0.7,k=3,
          rect = TRUE,
          k_colors = "jco",
          rect_border="jco",
          rect_fill=TRUE,
          type = "phylogenic",
          repel = TRUE,
          phylo_layout = "layout.mds")
fviz_dend(hc,cex=0.8,lwd = 0.7,k=3,
          rect = TRUE,
          k_colors = "jco",
          rect_border="jco",
          rect_fill=TRUE,
          type = "phylogenic",
          repel = TRUE,
          phylo_layout = "layout_as_tree")


# Define 8 custom colors (one for each cluster)
my_colors <- c(
  "#E41A1C",  # Red
  "#377EB8",  # Blue
  "#4DAF4A",  # Green
  "#984EA3",  # Purple
  "#FF7F00",  # Orange
  "#FFFF33",  # Yellow
  "#A65628",  # Brown
  "#F781BF"   # Pink
)

# Apply colors to branches and labels
hc <- hc %>%
  color_branches(k = 3, col = my_colors) %>%  # Color branches for 8 clusters
  color_labels(k = 3, col = my_colors) %>%    # Color labels for 8 clusters
  set("branches_lwd", 3) %>%                  # Line thickness
  set("branches_lty", 1)                      # Solid lines

# ================================================
# CUSTOM COLOR SECTION ENDS HERE
# ================================================

# Create circular dendrogram
circlize_dendrogram(
  hc,
  labels_track_height = 0.3,
  dend_track_height = 0.5
)

# Cluster analysis for k = 8 clusters
k <- 2
clusters <- cutree(hc, k = k)  # Cut the dendrogram into 8 clusters
print("Cluster assignments:")
print(clusters)

# Compute cluster means
# Combine cluster assignments with the data
data_with_clusters <- data_for_clustering
data_with_clusters$Cluster <- clusters

# Calculate mean of each feature per cluster
cluster_means <- data_with_clusters %>%
  group_by(Cluster) %>%
  summarise_all(mean, na.rm = TRUE)

print("Cluster means:")
print(cluster_means)

# Compute distances between cluster centroids
# Extract centroids (means) for each cluster, excluding the Cluster column
centroids <- as.matrix(cluster_means[, -1])  # Remove Cluster column

# Compute Euclidean distance between centroids
cluster_distances <- dist(centroids, method = "euclidean")

# Create a symmetric distance matrix with zeros on the diagonal
n_clusters <- k
dist_matrix <- matrix(0, nrow = n_clusters, ncol = n_clusters)
dist_matrix[lower.tri(dist_matrix)] <- as.vector(cluster_distances)
dist_matrix[upper.tri(dist_matrix)] <- t(dist_matrix)[upper.tri(dist_matrix)]  # Mirror to make symmetric

# Set row and column names as Cluster1, Cluster2, ..., Cluster8
cluster_names <- paste0("Cluster", 1:n_clusters)
rownames(dist_matrix) <- cluster_names
colnames(dist_matrix) <- cluster_names

# Print the formatted distance matrix
print("Cluster distance matrix:")
print(dist_matrix)

####Results
> library(readxl)
> WHEAT <- read_excel("D:/khem data analysis/WHEAT.xlsx", 
+     sheet = "Cluster ")
> View(WHEAT)
> library(readxl)
> library(dendextend)

---------------------
Welcome to dendextend version 1.19.0
Type citation('dendextend') for how to cite the package.

Type browseVignettes(package = 'dendextend') for the package vignette.
The github page is: https://github.com/talgalili/dendextend/

Suggestions and bug-reports can be submitted at: https://github.com/talgalili/dendextend/issues
You may ask questions at stackoverflow, use the r and dendextend tags: 
	 https://stackoverflow.com/questions/tagged/dendextend

	To suppress this message use:  suppressPackageStartupMessages(library(dendextend))
---------------------


Attaching package: ‘dendextend’

The following object is masked from ‘package:stats’:

    cutree

Warning message:
package ‘dendextend’ was built under R version 4.3.3 
> library(circlize)
========================================
circlize version 0.4.16
CRAN page: https://cran.r-project.org/package=circlize
Github page: https://github.com/jokergoo/circlize
Documentation: https://jokergoo.github.io/circlize_book/book/

If you use it in published research, please cite:
Gu, Z. circlize implements and enhances circular visualization
  in R. Bioinformatics 2014.

This message can be suppressed by:
  suppressPackageStartupMessages(library(circlize))
========================================

> library(dplyr)

Attaching package: ‘dplyr’

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

> # Load data
> data <- WHEAT
> # Convert to data frame and set row names
> data <- as.data.frame(data)
> rownames(data) <- data$`Accession name`
> # Select clustering columns
> data_for_clustering <- data[, c("DTH","DTF","DOM","PH","SPL","NOT","SPPS","SPSP","SPS","SL","THSW","TY"
+ )]
> # Scale data
> scaled_data <- scale(data_for_clustering)
> # Compute distance matrix
> d <- dist(scaled_data, method = "euclidean")
> # Perform hierarchical clustering
> hc <- as.dendrogram(hclust(d, method = "average"))
> ##Cluster dendrogram
> library(factoextra)
Loading required package: ggplot2
Welcome! Want to learn more? See two factoextra-related books at https://goo.gl/ve3WBa
Warning messages:
1: package ‘factoextra’ was built under R version 4.3.3 
2: package ‘ggplot2’ was built under R version 4.3.3 
> fviz_dend(hc,cex=0.3,lwd = 0.7,k=3,
+           k_colors = c("red","green","magenta"))
Warning message:
The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as of ggplot2
3.3.4.
ℹ The deprecated feature was likely used in the factoextra package.
  Please report the issue at <https://github.com/kassambara/factoextra/issues>.
This warning is displayed once every 8 hours.
Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated. 
> fviz_dend(hc,cex=0.5,lwd = 0.7,k=3,
+           k_colors = "jco")
> fviz_dend(hc,cex=0.5,lwd = 0.7,k=3,
+           rect = TRUE,
+           rect_border = "gray",
+           rect_fill = FALSE)
> fviz_dend(hc,cex=0.5,lwd = 0.7,k=3,
+           rect = TRUE,
+           k_colors = "jco",
+           rect_border="jco",
+           rect_fill=TRUE,
+           ggtheme=theme_bw())
> fviz_dend(hc,cex=0.5,lwd = 0.7,k=3,
+           rect = TRUE,
+           k_colors = c("red","green3","magenta"),
+           rect_border= c("red","green3","magenta"),
+           rect_fill=TRUE)
Error in if (color == "cluster") color <- "default" : 
  the condition has length > 1
> # Apply colors to branches and labels
> hc <- hc %>%
+   color_branches(k = 3, col = my_colors) %>%  # Color branches for 8 clusters
+   color_labels(k = 3, col = my_colors) %>%    # Color labels for 8 clusters
+   set("branches_lwd", 3) %>%                  # Line thickness
+   set("branches_lty", 1)                      # Solid lines
Error: object 'my_colors' not found
> # Create circular dendrogram
> circlize_dendrogram(
+   hc,
+   labels_track_height = 0.3,
+   dend_track_height = 0.5
+ )
> # Cluster analysis for k = 8 clusters
> k <- 2
> clusters <- cutree(hc, k = k)  # Cut the dendrogram into 8 clusters
> print("Cluster assignments:")
[1] "Cluster assignments:"
> print(clusters)
NGRCO 5102  NGRCO5103  NGRCO5104  NGRCO5188  NGRCO5189  NGRCO5190  NGRCO5191  NGRCO5192 
         1          2          1          1          1          1          1          1 
 NGRCO5193  NGRCO5194  NGRCO5195  NGRCO5196  NGRCO5197  NGRCO5198  NGRCO5199  NGRCO5200 
         1          1          1          1          1          1          1          1 
 NGRCO5201  NGRCO6199  NGRCO6200  NGRCO6201  NGRCO6230  NGRCO6231  NGRCO6232  NGRCO6233 
         1          1          1          1          1          1          1          1 
 NGRCO6234  NGRCO6235  NGRCO6236  NGRCO6237  NGRCO6238  NGRCO6239  NGRCO6240  NGRCO6241 
         1          1          1          1          1          1          1          1 
 NGRCO6242  NGRCO6243  NGRCO6244  NGRCO6245  NGRCO6246  NGRCO6247  NGRCO6248  NGRCO6249 
         1          1          1          1          1          1          1          1 
 NGRCO6250  NGRCO6251  NGRCO6252  NGRCO6253  NGRCO6254  NGRCO6255  NGRCO6256  NGRCO6257 
         1          1          1          1          1          1          1          1 
 NGRCO6258  NGRCO6259  NGRCO6260  NGRCO6261  NGRCO6262  NGRCO6263  NGRCO6264  NGRCO6265 
         1          1          1          1          1          1          1          1 
 NGRCO6266  NGRCO6267  NGRCO6268  NGRCO6269  NGRCO6270  NGRCO6448 NGRCO10285 NGRCO10286 
         1          1          1          1          1          1          1          1 
NGRCO10287 NGRCO10288 NGRCO10289 NGRCO10290 NGRCO10291 NGRCO10294 NGRCO10295 NGRCO10296 
         1          1          1          1          1          1          1          1 
NGRCO10299 NGRCO10300 NGRCO10301 NGRCO10302 NGRCO10303 NGRCO10284 NGRCO10292 NGRCO10304 
         1          1          1          1          1          1          1          1 
NGRCO10305 NGRCO10306 NGRCO10307 NGRCO10308 NGRCO10309 NGRCO10310 NGRCO10311 NGRCO10393 
         1          1          1          1          1          1          1          1 
   CO15402    CO15555 NGRCO00021 NGRCO00108 NGRCO00128 NGRCO00129 NGRCO00131 NGRCO00132 
         1          1          1          1          1          1          1          1 
NGRCO00134 NGRCO00135 NGRCO00136 NGRCO00137 NGRCO00138 NGRCO00145 NGRCO00146 NGRCO00147 
         1          1          1          1          1          1          1          1 
NGRCO00148 NGRCO00157 NGRCO00158 NGRCO00159 NGRCO00166 NGRCO00168 NGRCO00169 NGRCO00195 
         1          2          1          1          1          1          1          1 
NGRCO00208 NGRCO06348 NGRCO06349 NGRCO06350 
         1          1          1          1 
> # Compute cluster means
> # Combine cluster assignments with the data
> data_with_clusters <- data_for_clustering
> data_with_clusters$Cluster <- clusters
> # Calculate mean of each feature per cluster
> cluster_means <- data_with_clusters %>%
+   group_by(Cluster) %>%
+   summarise_all(mean, na.rm = TRUE)
> print("Cluster means:")
[1] "Cluster means:"
> print(cluster_means)
# A tibble: 2 × 13
  Cluster   DTH   DTF   DOM    PH   SPL   NOT  SPPS  SPSP   SPS    SL  THSW    TY
    <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
1       1  108.  113.  144.  93.0  9.14  8.27  9.47  2.92  41.2  5.96  36.6  1.42
2       2  128   136.  158   56.8  6.96  9     7.9   2.1   24.9  6.31  25.0  1.87
> # Compute distances between cluster centroids
> # Extract centroids (means) for each cluster, excluding the Cluster column
> centroids <- as.matrix(cluster_means[, -1])  # Remove Cluster column
> # Compute Euclidean distance between centroids
> cluster_distances <- dist(centroids, method = "euclidean")
> # Create a symmetric distance matrix with zeros on the diagonal
> n_clusters <- k
> dist_matrix <- matrix(0, nrow = n_clusters, ncol = n_clusters)
> dist_matrix[lower.tri(dist_matrix)] <- as.vector(cluster_distances)
> dist_matrix[upper.tri(dist_matrix)] <- t(dist_matrix)[upper.tri(dist_matrix)]  # Mirror to make symmetric
> # Set row and column names as Cluster1, Cluster2, ..., Cluster8
> cluster_names <- paste0("Cluster", 1:n_clusters)
> rownames(dist_matrix) <- cluster_names
> colnames(dist_matrix) <- cluster_names
> # Print the formatted distance matrix
> print("Cluster distance matrix:")
[1] "Cluster distance matrix:"
> print(dist_matrix)
         Cluster1 Cluster2
Cluster1  0.00000 53.50554
Cluster2 53.50554  0.00000
