data_path <- "your_data_path"
#install.packages("readxl")
library(readxl)

all_files <- list.files(data_path, pattern = "^data_05.*\\.xlsx$", full.names = TRUE)

#make a list
data_list <- list()

# for loop 
for (i in seq_along(all_files)) {
  data_list[[i]] <- read_excel(all_files[i], sheet = 1)
}

names(data_list) <- paste0("data", 1:length(all_files))

####### FOR Bar plot figures | Fig 1 ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

for (i in seq_along(all_files)) {
  
  # Read file
  data <- read_excel(all_files[i], sheet = 1)
  
  # Remove NA 
  cleaned_data <- data %>%
    filter(!is.na(Ratio_PSI))

  #Filtration non-numeric values
  cleaned_data <- cleaned_data %>%
    mutate(Ratio_PSI = as.numeric(Ratio_PSI)) %>%
    filter(!is.na(Ratio_PSI))
  
  #Computituonal AS events
  as_events <- cleaned_data %>%
    group_by(gp1, gp2) %>%
    summarise(avg_ratio_psi = mean(Ratio_PSI))
  
  # Bar Plot
  p <- ggplot(as_events, aes(x = interaction(gp1, gp2), y = avg_ratio_psi, fill = gp2)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = paste0("Average Alternative Splicing Events by Stage: ", basename(all_files[i])), 
         x = "Stage Comparison", y = "Average Ratio PSI") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save a plot
  ggsave(filename = paste0("AS_Barplot_", i, ".png"), plot = p, width = 10, height = 6)
  
  # List of results
  data_list[[i]] <- list(file_name = basename(all_files[i]), plot = p)
}

####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 
# Figure 2 | Box plots 
data$FDR <- p.adjust(data$p_value, method = "BH")


selected_genes <- c("AMD1", "ARMC8", "ATP6V1C2", "TDP1", "THOC2", "UBAP2L") # genes are identified from pre-analysis results
data <- data{i}
filtered_data <- data %>% filter(gene_name %in% selected_genes)

data_tidy <- filtered_data %>% 
  gather(key = "Group", value = "PSI", PSI_gp1, PSI_gp2)

# Grup info is merged
data_tidy <- data_tidy %>% 
  mutate(Group_Label = paste(gp1, "vs", gp2))

# Box plot  - with facet_wrap for all gp1 ve gp2 ccombination:
ggplot(data_tidy, aes(x = gene_name, y = PSI, fill = Group)) +
  geom_boxplot() +
  labs(title = "PSI Changes for Selected Genes", 
       x = "Genes", y = "PSI Values") +
  theme_minimal() +
  facet_wrap(~ Group_Label, scales = "free_y", ncol = 2)  # Multi box plot 


####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 
# Figure 3 | Heatmap

library(ggplot2)
library(pheatmap)
library(dendextend)
library(reshape2)
library(dplyr)

selected_genes <- c("AMD1", "ARMC8", "ATP6V1C2", "TDP1", "THOC2", "UBAP2L")


merged_data <- rbind(data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11)

scats_data <- merged_data
# change form
scats_long <- melt(scats_data, id.vars = c("gp1", "gp2", "gene_name", "AS_exons"), 
                   variable.name = "metric", value.name = "value")

# Selected gene_name
scats_long_selected <- scats_long[scats_long$gene_name %in% selected_genes, ]

#define stage of GP1
gp1_stage <- "Eblastocyst"  # TODO

# Filtration of GP1
subset_data <- scats_long_selected[scats_long_selected$gp1 == gp1_stage, ]

# define stage of GP2
gp2_stages <- unique(subset_data$gp2)

# Heatmap function
plot_heatmap <- function(data_matrix, title) {
  pheatmap(data_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
           main = title,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           show_rownames = TRUE, show_colnames = TRUE)
}

# make a matrix 
if (nrow(subset_data) > 0) {
  # Data matrix: rows are genes, columns are GP2 stages
  data_matrix <- acast(subset_data, gene_name ~ gp2, value.var = "value")
  
  # make a heatmap
  title <- paste("Heatmap for GP1 =", gp1_stage)
  plot_heatmap(data_matrix, title)
  
} else {
  print(paste("Yetersiz veri - GP1 Türü:", gp1_stage))
}







