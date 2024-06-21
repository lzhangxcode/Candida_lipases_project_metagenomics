# this script is used for visualization of the diversity 
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(metagenomeSeq)
library(pheatmap)
library(ggVennDiagram)
library(viridis)

input_dir <- "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output/diversity/"
output_dir <- "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output/diversity/plots/"

####################################################################

#                     boxplot for alpha diversity                  # 

####################################################################


# 1. inputs 
input_alpha <- read_tsv(str_c(input_dir, "/alpha_diversity/1_alpha_diversity.tsv")) %>% 
  select(samples, chao1, diversity_gini_simpson, diversity_shannon) %>% 
  rename(simpson = diversity_gini_simpson, shannon = diversity_shannon)
input_paried <- read_tsv(str_c(input_dir, "/alpha_diversity//1_alpha_diversity_paired_metatable.tsv"))
input_alpha_statistic <- read_tsv(str_c(input_dir, "/alpha_diversity/1_alpha_diversity_wilcoxon_paired.tsv"))

# 2. visualization 

vec3 <- c("ABX", "No_ABX")
vec4 <- c("shannon", "simpson", "chao1")

# creating a data frame of vector1 and vector2 combinations
compar_alpha_diversity_table <- expand.grid("abx" = vec3, "alpha" = vec4) %>% add_column(output = str_c(output_dir, "/alpha_diversity_plots/"))

draw_boxplot_time_point <- function(input_alpha_number = input_alpha, metatable = input_paried, abx = NA, index = NA, output_dir = output){
  
  print(abx)
  print(index)
  
  input_boxplot <- input_alpha_number %>% 
    left_join(., metatable, by = c('samples' = "sample")) %>% 
    filter(antibiotic == abx) 
  
  input_boxplot <- input_boxplot[, c(1, which(colnames(input_boxplot)==index), 7, 8, 9)]
  
  
  label <- paste(abx, index, sep = "_")
  colnames(input_boxplot)[2] <- "index"
  input_boxplot$time <- as.character(input_boxplot$time)
  input_boxplot <- input_boxplot %>% arrange(time, group, diet)
  stat.test <- input_boxplot %>% group_by(group, diet) %>% rstatix::wilcox_test(., index ~ time, paired = T, p.adjust.method = "fdr") %>% rstatix::add_xy_position(x = "time")
  
  p_b <- ggpaired(
    input_boxplot %>% mutate(time_group = str_c("Day",time, "_", group)), 
    x = "time", y = "index", color = "time_group", line.color = "gray", line.size = 0.4,
    facet.by = c("group", "diet"), scales = "free", add = c("jitter")
  ) +
    stat_pvalue_manual(stat.test, label = "p", hide.ns = TRUE) + 
    ggtitle(label) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab(index) +
    scale_color_manual(values=c("#ff006645", "#00000045", "#ff00668c", "#0000008c"))
  
  # 009392 # E9E29C
  
  write_tsv(stat.test, paste0(output_dir, "/", label, ".tsv"))
  ggsave(paste0(output_dir, "/", label, ".pdf"), p_b, width = 7, height = 8)
  
  return(p_b)
  
  
}

p_b <- compar_alpha_diversity_table %>% rowwise() %>% 
  mutate(plot = list(draw_boxplot_time_point(abx = abx, index = alpha, output_dir = output)))


####################################################################

#        pcoa plot for beta diversity - normal day14 vs day0       # 

####################################################################

# input 
beta_statistics_genera <- readRDS(str_c(input_dir, "/beta_diversity/2_beta_statistic_genera.rds"))
beta_statistic_result_unnest <- read_tsv(str_c(input_dir, "/beta_diversity/2_beta_statistic_tibble.tsv"))

# input-fraction 
beta_statistics_genera_fraction <- readRDS(str_c(input_dir, "/beta_diversity/2_beta_statistic_genera_fraction.rds"))
beta_statistic_result_unnest_fraction <- read_tsv(str_c(input_dir, "/beta_diversity/2_beta_statistic_tibble_fraction.tsv"))


# input-fraction - robust aitchison 
beta_statistics_genera_fraction_v2 <- readRDS(str_c(input_dir, "/beta_diversity/2_beta_statistic_genera_fraction_v2.rds"))
beta_statistic_result_unnest_fraction_v2 <- read_tsv(str_c(input_dir, "/beta_diversity/2_beta_statistic_tibble_fraction_v2.tsv"))


# function for draw the pcoa plot 
PcoA_plot_time <- function(points_df = NA, eig = NA, p_value = NA, rsqaure = NA, tag = NA){
  
  
  tag <- str_replace_all(tag, "AND", "_") %>% str_remove_all(., "_all")
  if(is.na(str_extract(tag, "Mutant")) ){
    tag <- str_replace_all(tag, "WT", "WT_")
    p <- ggpubr::ggscatter(points_df, x = "x", y = "y",
                           color = "time_md", shape = "time_md", fill = "time_md", ellipse = TRUE, palette = c("#00000045", "#0000008c")) + 
      labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
      annotate(geom="text", x = max(points_df$x), y = (max(points_df$y) + 1.5*sd(points_df$y)), label = list(bquote(P == .(as.vector(p_value))))) + # x = 0.2, y = 0.3, 
      annotate(geom="text", x = max(points_df$x), y = (max(points_df$y) + sd(points_df$y)), label = list(bquote(R^2 == .(as.vector(rsqaure))))) + # x = 0.2, y = 0.3, 
      ggtitle(tag) + 
      theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
    
  }else{
    tag <- str_replace_all(tag, "Mutant", "Mutant_")
    
    p <- ggpubr::ggscatter(points_df, x = "x", y = "y",
                           color = "time_md", shape = "time_md", fill = "time_md", ellipse = TRUE, palette = c("#ff006645", "#ff00668c")) + 
      labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
      annotate(geom="text", x = max(points_df$x), y = (max(points_df$y) + 1.5*sd(points_df$y)), label = list(bquote(P == .(as.vector(p_value))))) + # x = 0.2, y = 0.3, 
      annotate(geom="text", x = max(points_df$x), y = (max(points_df$y) + sd(points_df$y)), label = list(bquote(R^2 == .(as.vector(rsqaure))))) + # x = 0.2, y = 0.3, 
      ggtitle(tag) + 
      theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
    
  }
  
  
  return(p)
}



# function for calculate the points for pcoa plot 
PCoA_plots_all <- function(distance.bray = NA, permanova_result = NA, output_dir = NA, metadat_sel = NA, filetag = NA){
  
  set.seed(1234)
  pcoa <- cmdscale(distance.bray, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
  points <- as.data.frame(pcoa$points) # get coordinate string, format to dataframme
  colnames(points) <- c("x", "y", "z")
  eig <- pcoa$eig
  
  metadat_sel <- input_paried %>% filter(sample %in% rownames(points)) %>% 
    column_to_rownames(., var = "sample")
  metadat_sel <- metadat_sel[rownames(points),]
  metadat_sel$time_md <- as.character(metadat_sel$time)
  
  points <- cbind(points, metadat_sel) 
  
  group_pcoa <- PcoA_plot_time(points_df = points, eig = eig,
                               p_value = permanova_result$`Pr(>F)`[1],
                               rsqaure = round(permanova_result$R2[1], 3),
                               tag = filetag)
  
  
  
  system(str_c("mkdir -p ", output_dir, "/beta_diversity_plots/tmp/"))
  saveRDS(points, str_c(output_dir, "/beta_diversity_plots/tmp/pcoa_points_", filetag, ".rds"))
  saveRDS(eig, str_c(output_dir, "/beta_diversity_plots/tmp/pcoa_eig_", filetag, ".rds"))
  
  
  return(group_pcoa)    
  
}

# function for exceute function: PCoA_plots_all
produce_pcoa_plots <- function(distance_result = NA, output_dir = NA, metadat_sel = NA, category = NA){
  library(ggpubr)
  
  pcoa_result_bray <- list()
  pcoa_result_bray <- PCoA_plots_all(distance.bray = distance_result[["bray_distance"]], 
                                     permanova_result = distance_result[["bray_adonis_time"]], 
                                     output_dir = output_dir, 
                                     metadat_sel = metadat_sel, 
                                     filetag = paste(category, "bray_curtis", sep = ""))
  pcoa_result_aitchison <- list()
  pcoa_result_aitchison <- PCoA_plots_all(distance.bray = distance_result[["aitchison"]], 
                                          permanova_result = distance_result[["aitchison_adonis_time"]], 
                                          output_dir = output_dir, 
                                          metadat_sel = metadat_sel,
                                          filetag = paste(category, "aitchison", sep = ""))
  pcoa_result_all <- list()
  pcoa_result_all[["bray"]] <- pcoa_result_bray
  pcoa_result_all[["aitchison"]] <- pcoa_result_aitchison
  
  
  return(pcoa_result_all)
  
}


# function for combine the pcoa plots 

produce_combined_plots <- function(input_list_figures = NA){
  new_plot_list <- list()
  for (each_plot in names(input_list_figures)){
    new_plot_list[[each_plot]] <- input_list_figures[[each_plot]][["bray"]]
  }
  
  bray_combined <- ggarrange(plotlist = new_plot_list, nrow = 2, ncol = 4, widths = 15)
  
  new_plot_list <- list()
  for (each_plot in names(input_list_figures)){
    new_plot_list[[each_plot]] <- input_list_figures[[each_plot]][["aitchison"]]
  }
  aitchison_combined <- ggarrange(plotlist = new_plot_list, nrow = 2, ncol = 4, widths = 15)
  
  combined <- list(bray_combined, aitchison_combined)
  return(combined)
  
}


# fraction - robust aitchison 
beta_plot_time_fraction_v2 <- list()
beta_plot_time_combined_fraction_v2 <- list()
for (g in names(beta_statistics_genera_fraction_v2)){
  label <- g
  
  beta_plot_time_fraction_v2[[g]] <- produce_pcoa_plots(distance_result = beta_statistics_genera_fraction_v2[[g]], 
                                                     output_dir = output_dir,
                                                     metadat_sel = input_paried, 
                                                     category = label)
}

beta_plot_time_combined_fraction_v2 <- produce_combined_plots(input_list_figures = beta_plot_time_fraction_v2)
ggsave(paste(output_dir, "/beta_diversity_plots_fraction_v2/beta_diversity_bray_fraction.pdf", sep = "/"), beta_plot_time_combined_fraction_v2[[1]], height = 8, width = 15)
ggsave(paste(output_dir, "/beta_diversity_plots_fraction_v2/beta_diversity_ait_fraction.pdf", sep = "/"), beta_plot_time_combined_fraction_v2[[2]], height = 8, width = 15)


save.image(str_c(output_dir, "genera_vis.RData"))



