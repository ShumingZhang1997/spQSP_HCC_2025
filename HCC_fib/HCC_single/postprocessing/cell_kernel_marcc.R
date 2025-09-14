###############################################
#  do_Cellular_Influence_(CI)_quantification  #
###############################################
library(spatstat)
library(dplyr)
library(reshape2)
#library(ggpubr)

do_CI_quantification<- function(expr = df_rm, # should contain sample_id, X, Y position, celltype in the df
                                kernels ="gaussian",
                                clusterlevels = clusterlevels,
                                sigma = 10){
  require(spatstat);
  
  spatwt<- c()
  expr_k<- expr
  coord <- expr_k[, c("X_position", "Y_position")]
  
  rownames(coord)<- paste0(rownames(expr_k), ":", expr_k$celltype)
  
  gaussian<- data.frame(matrix(0, nrow = nrow(coord), ncol = length(clusterlevels)))
  colnames(gaussian) <- clusterlevels
  
  coord <- cbind(coord, gaussian)
  
  sapply(1:nrow(coord), function(r)
    coord[r, as.character(lapply(strsplit(as.character(rownames(coord)[r]), split = ":"),"[[",2))] <<- 1)
  
  # calculation of spatial weight 
  mywin <- owin(xrange = range(coord$X_position), yrange = range(coord$Y_position))
  spObj <- ppp(x=coord$X_position,y=coord$Y_position,window = mywin, marks = coord[,3:length(coord)])
  
  # adjust sigma 
  # leaveoneout = T -> removing the self
  K <- Smooth(spObj,kernel = kernels,at = "points",leaveoneout = T, sigma=sigma) 
  rownames(K)<-rownames(coord)
  K<- as.data.frame(K)
  #K$sample_id <- i 
  
  spatwt<- rbind(spatwt, K)
  
  
  return (spatwt) 
}

#preprocessing

# Initialize an empty list to store flattened matrices
print(getwd())
flattened_list <- list()
file_path <- "./snapShots/cell_0.csv"

if (!file.exists(file_path)) {
  message("Skipping sample: file not found.")
  next  # This is the equivalent of "continue"
}

HCCdata <- read.csv(file_path)

clusterlevels_all <- c("Tumor", "CD8T", "CD4T", "Tregs", "MDSCs", "Macrophages", "Fibroblasts", "Endothelials")

interaction_z_list <- list()

for (j in seq(0, 49)) {
  temp_HCC <- HCCdata %>% filter(z == as.character(j)) %>%
    rename(X_position = x, Y_position = y) %>%
    filter(!is.na(X_position) & !is.na(Y_position)) %>%
    mutate(
      X_position = (X_position + 1) * 20 + rnorm(n(), mean = 0, sd = 2),
      Y_position = (Y_position + 1) * 20 + rnorm(n(), mean = 0, sd = 2),
      celltype = case_when(
        Type == 1 ~ "Tumor",
        Type == 2 ~ "CD8T",
        Type == 3 & State == 9 ~ "CD4T",
        Type == 3 & State == 10 ~ "Tregs",
        Type == 4 ~ "MDSCs",
        Type == 5 ~ "Macrophages",
        Type == 6 ~ "Fibroblasts",
        Type == 8 ~ "Endothelials",
        TRUE ~ "Other"
      )
    )
  
  if (nrow(temp_HCC) == 0 || all(temp_HCC$celltype == "Tumor")) next
  
  spatial_weights <- do_CI_quantification(temp_HCC, clusterlevels = clusterlevels_all)
  spatial_weights$celltype <- sapply(strsplit(rownames(spatial_weights), ":"), `[`, 2)
  
  interaction_matrix <- spatial_weights %>%
    group_by(celltype) %>%
    summarise(across(all_of(clusterlevels_all), ~ mean(.x, na.rm = TRUE))) %>%
    as.data.frame()
  
  rownames(interaction_matrix) <- interaction_matrix$celltype
  interaction_matrix <- interaction_matrix[, -1, drop = FALSE]
  
  # Add missing rows
  missing_rows <- setdiff(clusterlevels_all, rownames(interaction_matrix))
  if (length(missing_rows) > 0) {
    add_matrix <- matrix(0, nrow = length(missing_rows), ncol = length(clusterlevels_all))
    rownames(add_matrix) <- missing_rows
    colnames(add_matrix) <- clusterlevels_all
    interaction_matrix <- rbind(interaction_matrix, add_matrix)
  }
  
  # Add missing columns
  missing_cols <- setdiff(clusterlevels_all, colnames(interaction_matrix))
  if (length(missing_cols) > 0) {
    for (mc in missing_cols) {
      interaction_matrix[[mc]] <- 0
    }
  }
  
  # Reorder both rows and columns
  interaction_matrix <- interaction_matrix[clusterlevels_all, clusterlevels_all]
  
  interaction_z_list[[as.character(j)]] <- interaction_matrix
}

# Average across all z-slices
if (length(interaction_z_list) > 0) {
  avg_matrix <- Reduce("+", interaction_z_list) / length(interaction_z_list)
  flat_matrix <- as.data.frame(as.table(as.matrix(avg_matrix)))
  flat_matrix$Source_Neighbor <- paste0("Reference_", flat_matrix$Var1, "_Weight_", flat_matrix$Var2)
  flat_matrix <- flat_matrix[, c("Source_Neighbor", "Freq")]
  colnames(flat_matrix)[2] <- paste("Core", 1)
  flattened_list[[1]] <- flat_matrix
}


flattened_list <- flattened_list[!sapply(flattened_list, is.null)]
spatial_merged_mat <- Reduce(function(x, y) merge(x, y, by = "Source_Neighbor", all = TRUE), flattened_list)
spatial_merged_mat[is.na(spatial_merged_mat)] <- 0
spatial_merged_df <- as.data.frame(t(spatial_merged_mat))
rownames_temp <- rownames(spatial_merged_df)
colnames(spatial_merged_df) <- spatial_merged_df[1, ]
spatial_merged_df <- spatial_merged_df[-1, ]
rownames(spatial_merged_df) <- rownames_temp[-1]
spatial_merged_df <- spatial_merged_df %>% mutate_all(~ as.numeric(as.character(.)))

#print(paste("spatial_merged_df$Reference_CD4T_Weight_Tumor:", spatial_merged_df$Reference_CD4T_Weight_Tumor))
#print(paste("spatial_merged_df dimensions: ", dim(spatial_merged_df)))
#print(paste("Column names: ", colnames(spatial_merged_df)))
