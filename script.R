dtf <- read.csv("2026-09-09_data_PROT_RNA_TP.csv", 
                stringsAsFactors = TRUE)
plot(dtf[, 7:10])
summary(dtf[, 7:10])
dtf_rEV <- dtf[, c(2, 7:9)]
summary(dtf_rEV)

# PCA on rEVij
dtf_rEV <- dtf[, c(2, 7:9)]
summary(dtf_rEV)
library(FactoMineR)
pca1 <- PCA(dtf_rEV, scale.unit = FALSE,
            graph = FALSE, quali.sup = 1)
plot(pca1, choix = "var")
cor(dtf_rEV[, -1])
plot(pca1, choix = "ind", label = "none", 
     habillage = 1)

# PCA on rEVij, Num_Atoms and Rg
Rg <- sqrt(dtf$Eigval.1 +
             dtf$Eigval.2 +
             dtf$Eigval.3)
dtf_all <- data.frame(dtf_rEV, 
                      size = dtf$Num_Atoms, 
                      Rg = Rg)
summary(dtf_all)
pca2 <- PCA(dtf_all, scale.unit = TRUE,
            graph = FALSE, quali.sup = 1)
summary(pca2)
plot(pca2, choix = "var", axes = c(1, 3))
plot(pca2, choix = "ind", habillage = 1, 
     label = "none", axes = c(1, 3))
cor(dtf_all[, -1])
corrplot::corrplot(cor(dtf_all[, -1]))

# hierarchical clustering on rEVij
dist_rEV <- dist(dtf_rEV[, -1])
par(mfrow = c(2, 2))
hclust_rEV_complet <- hclust(dist_rEV, method = "complet")
plot(hclust_rEV_complet, hang = -1,
     labels = FALSE)
hclust_rEV_ward <- hclust(dist_rEV, method = "ward.D2")
plot(hclust_rEV_ward, hang = -1,
     labels = FALSE)
hclust_rEV_average <- hclust(dist_rEV, method = "average")
plot(hclust_rEV_average, hang = -1,
     labels = FALSE)
par(mfrow = c(1, 1))
cut_complet <- cutree(hclust_rEV_complet, k = 6)
cut_ward <- cutree(hclust_rEV_ward, k = 5)
cut_average <- cutree(hclust_rEV_average, k = 6)
par(mfrow = c(1, 1))

## how many clusters
library(cluster)
par(mfrow = c(2, 2))
sil_complet <- silhouette(cut_complet, dist = dist_rEV)
plot(sil_complet, border = NA)
par(mfrow = c(1, 1))

library(NbClust)
nb_k_ward <- NbClust(dtf_rEV[, -1], method = "ward.D2")
nb_k_complet <- NbClust(dtf_rEV[, -1], method = "complete")

cut_ward4 <- cutree(hclust_rEV_ward, k = 4)
table(cut_ward4)
plot(hclust_rEV_ward, hang = -1, labels = FALSE)
rect.hclust(hclust_rEV_ward, k = 4)
sil_ward <- silhouette(x = cut_ward4, 
                          dist = dist_rEV)
plot(sil_ward, border = NA)

## stability of the clustering = bootstrap
set.seed(2026)

cut_ward <- matrix(nrow = 100, ncol = nrow(dtf_rEV))
for (b in 1:100) {
  # create new datas from dtf_rEV
  ind_boot <- sample(1:nrow(dtf_rEV), 
                     size = nrow(dtf_rEV),
                     replace = TRUE)
  dtf_temp <- dtf_rEV[ind_boot,]
  dist_temp <- dist(dtf_temp[, -1])
  hclust_temp_ward <- hclust(dist_temp, 
                            method = "ward.D2")
  cut_ward[b, ] <- cutree(hclust_temp_ward, k = 4)
}

jaccard <- function(vec1, vec2) {
  length(intersect(vec1, vec2)) / length(union(vec1, vec2))
}

### on est dans une impasse de code
### on utilise donc une fonction toute faite.

## avec une fonction de boostrap toute faite.
set.seed(2026)
k <- 3
cboot_k <- fpc::clusterboot(dist_rEV, B = 100, 
                        distances = TRUE,
                        clustermethod = disthclustCBI,
                        method = "ward.D2", k = k,
                        count = FALSE)

cboot_k$bootmean  # Jaccard moyen par cluster -- equivalent a stabilite_k4


## silhouette and Jaccard
k <- 2
### silhouette
cut_ward <- cutree(hclust_rEV_ward, k = k)
sil_ward <- silhouette(x = cut_ward, 
                       dist = dist_rEV)
plot(sil_ward, border = NA)
### jaccard
set.seed(2026)
cboot_k <- fpc::clusterboot(dist_rEV, B = 100, 
                            distances = TRUE,
                            clustermethod = disthclustCBI,
                            method = "ward.D2", k = k,
                            count = FALSE)

cboot_k$bootmean  # Jaccard moyen par cluster -- equivalent a stabilite_k4

# PCA and clustering
# hierarchical clustering on rEVij
dist_all <- dist(scale(dtf_all[, -1]))
hclust_all_ward <- hclust(dist_all, method = "ward.D2")
cut_ward_5d <- cutree(hclust_all_ward, k = 4)

attributes(pca1$ind)
head(pca1$ind$coord)
## 3 and 5 descr on pca1
par(mfrow = c(1,2))
plot(pca1$ind$coord[, 1], pca1$ind$coord[, 2],
     pch = 20, col = cut_ward4)
abline(h = 0, v = 0, lty = 2)
plot(pca1$ind$coord[, 1], pca1$ind$coord[, 2],
     pch = 20, col = cut_ward_5d)
abline(h = 0, v = 0, lty = 2)
par(mfrow = c(1,1))

## 5 descr on pca2
plot(pca2$ind$coord[, 1], pca2$ind$coord[, 2],
     pch = 20, col = cut_ward_5d)
abline(h = 0, v = 0, lty = 2)

factoextra::fviz_pca_biplot(pca1)

# heatmap
library(pheatmap)
pheatmap(dtf_rEV[, -1])
pheatmap(scale(dtf_all[, -1]))

# describe the clusters
aggregate(dtf_all[, -1], by = list(cluster = factor(cut_ward_5d)), FUN = mean)
aggregate(dtf_all[, -1], by = list(cluster = factor(cut_ward_5d)), FUN = sd)
