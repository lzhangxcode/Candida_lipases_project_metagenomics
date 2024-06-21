# this script is for doing alpha diversity analysis & beta diversity analysis for the genera level 

#####################################################

# alpha diversity

#####################################################

load("/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output/taxa_ditribution/kingdom_distribution_all_taxa.RData")
library(tidyverse)
library(phyloseq)
library(metagenomeSeq)
library(ggVennDiagram)
library(ggpubr)

physeq_object_bac <- subset_taxa(physeq_object, Kingdom == "k__Bacteria")
physeq_object_bac_rm_cf <- physeq_object_bac

combined_genera <- function(input_phyloseq = NA){
  
  input_phyloseq_genus <- tax_glom(input_phyloseq, "Genus")
  
  number_of_genera <- tax_table(input_phyloseq) %>% as.data.frame(.) %>% .$Genus %>% unique()
  
  number_of_genera_phyloseq <- tax_table(input_phyloseq_genus) %>% as.data.frame(.) %>% .$Genus %>% unique() 
  
  print(setdiff(number_of_genera, number_of_genera_phyloseq))
  print(length(number_of_genera))
  print(length(number_of_genera_phyloseq))
  
  
  taxa_names(input_phyloseq_genus) <- taxa_names(input_phyloseq_genus) %>% str_extract(., pattern = "^k__.*g__.*\\|s__") %>% str_remove(., "\\|s__") 
  
  return(input_phyloseq_genus)
  
}

physeq_object_bac_rm_cf_list <- list("kneaddata_custom" = physeq_object_bac_rm_cf)
phyloseq_objects_genera_combined_from_species <- lapply(physeq_object_bac_rm_cf_list, combined_genera)
output_dir <- "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output/diversity/"
saveRDS(phyloseq_objects_genera_combined_from_species, str_c(output_dir, "phyloseq_objects_genera_combined_from_species.rds"))
rm(list = ls())

output_dir <- "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output/diversity/"
phyloseq_objects_genera_combined_from_species <- readRDS(str_c(output_dir,"phyloseq_objects_genera_combined_from_species.rds"))
metatable <- read_tsv("~/Documents/202210_mouse_project/04_diversity/metatable_paired.tsv")

pathway_alpha_diversity <- microbiome::alpha(phyloseq_objects_genera_combined_from_species$kneaddata_custom, index = "all")


input_alpha <- pathway_alpha_diversity
input_paired <- metatable


make_compar_table <- function(input_paired = NA){
  
  compar_table <- input_paired %>% 
    select(diet, antibiotic, group) %>% 
    unique() 
  
  return(compar_table)
}


compare_alpha_diversity_wilcox_paried <- function(input_alpha = NA, input_paired = NA, diet_choose = NA, antibiotic_choose = NA, group_choose = NA){
  
  # first extract alpha diversity for paired samples 
  input_alpha_sel <- input_alpha %>% 
    rownames_to_column(., var = "sample") %>%
    select(sample, chao1, diversity_gini_simpson, diversity_shannon) %>% 
    filter(sample %in% input_paired$sample) %>% 
    left_join(., input_paired, by = c("sample" = "sample")) %>% 
    filter(diet == diet_choose, antibiotic == antibiotic_choose, group == group_choose) %>% 
    arrange(time, sample)
  
  
  # second do wilcox paired test 
  
  wilcox_paired_table1 <- rstatix::wilcox_test(chao1 ~ time, data = input_alpha_sel, paired = T, p.adjust.method = "fdr")
  wilcox_paired_table2 <- rstatix::wilcox_test(diversity_gini_simpson ~ time, data = input_alpha_sel, paired = T, p.adjust.method = "fdr")
  wilcox_paired_table3 <- rstatix::wilcox_test(diversity_shannon ~ time, data = input_alpha_sel, paired = T, p.adjust.method = "fdr")
  
  
  wilcox_paired_result <- rbind(wilcox_paired_table1, wilcox_paired_table2, wilcox_paired_table3)
  return(wilcox_paired_result)
  
  
}



compar_table <- make_compar_table(input_paired)
alpha_diversity_wilcox_paried <-  compar_table %>% 
  rowwise() %>% 
  mutate(table = list(compare_alpha_diversity_wilcox_paried(input_alpha = input_alpha, input_paired = input_paired, diet_choose = diet, antibiotic_choose = antibiotic, group_choose = group)))


# write the result into tsv file 
write_tsv(alpha_diversity_wilcox_paried_table, str_c(output_dir, "/alpha_diversity/1_alpha_diversity_wilcoxon_paired.tsv"))
write_tsv(input_alpha %>% rownames_to_column(., var = "samples"), str_c(output_dir, "/alpha_diversity/1_alpha_diversity.tsv"))
write_tsv(input_paired, str_c(output_dir, "/alpha_diversity/1_alpha_diversity_paired_metatable.tsv"))



############################################
#     beta diversity  - day14 vs day0      # 
############################################

# beta diversity - original comparison day 14 vs day 0 

## produce input for calculating distance 

prepare_beta_diversity_input <- function(input = NA, metadat = NA){
  
  beta_diversity_input <- input %>% 
    dplyr::select(., c("name", metadat$sample)) %>% 
    column_to_rownames(., var = "name") #
  
  beta_diversity_input <- beta_diversity_input[rowSums(abs(beta_diversity_input)) != 0,]
  
  
  beta_diversity_input_t <- beta_diversity_input %>% t(.)
  
  
  input_list <- list()
  input_list[[1]] <- metadat
  input_list[[2]] <- beta_diversity_input_t
  
  return(input_list)
}



condition_separate_also_time <- function(input = NA, table = NA){
  beta_diversity_result <- list()
  for (abx in unique(table$antibiotic)){
    print(abx)
    for (diet_target in unique(table$diet)){
      print(diet_target)
      for (t in unique(table$group)){
        print(t)
        table_temp <- table %>% dplyr::filter(antibiotic == abx & diet == diet_target & group == t)
        input_temp <- input %>% dplyr::select(name, table_temp$sample)
        label <- paste(abx, diet_target, t, sep = "AND")
        
        beta_diversity_result[[label]] <- prepare_beta_diversity_input(input = input_temp, metadat = table_temp)
      }
      
    }
  }
  
  
  return(beta_diversity_result)
}


genera_profile_beta_diversity <- phyloseq_objects_genera_combined_from_species$kneaddata_custom %>% otu_table() %>% as.data.frame() %>% rownames_to_column(., var = "name")
genera_separate_profile <- list()
genera_separate_profile <- condition_separate_also_time(input = genera_profile_beta_diversity, table = input_paired)


# convert the count into fraction table 
phy_genus_relative <- microbiome::transform(phyloseq_objects_genera_combined_from_species$kneaddata_custom, trans = "compositional")
genera_profile_beta_diversity_fraction <- phy_genus_relative %>% otu_table() %>% as.data.frame() %>% rownames_to_column(., var = "name")
genera_separate_profile_fraction <- list()
genera_separate_profile_fraction <- condition_separate_also_time(input = genera_profile_beta_diversity_fraction, table = input_paired)



## calculate distance 


distance_produce <- function(beta_diversity = NA){ # 1st list metadata; 2nd list ec_profile
  
  
  metadat_sel <- beta_diversity[[1]]
  
  profile <- beta_diversity[[2]]
  profile <- profile[metadat_sel$sample, ]
  
  
  # distance calculation 
  set.seed(1234)
  distance.bray <- profile %>%
    vegan::vegdist(.,method = 'bray') %>%
    as.matrix() 
  
  set.seed(1234)
  aitchson_input <- profile
  min_tmp <- min(aitchson_input[aitchson_input>0])
  aitchson_input[aitchson_input == 0] <- min_tmp/2
  distance.aitchison <- as.matrix(robCompositions::aDist(aitchson_input))
  
  
  # separate adonis for each factor - bray curtis 
  metadat_sel$time <- as.character(metadat_sel$time)
  print(metadat_sel$time)
  permanova_result_bray_time <- vegan::adonis2(distance.bray ~ time, data = metadat_sel)

  
  # separate adonis for each factor - aitchison 
  permanova_result_aitchison_time <- vegan::adonis2(distance.aitchison ~ time, data = metadat_sel)

  
  # beta dispersion 
  bray_betadispersion <- vegan::betadisper(as.dist(distance.bray), metadat_sel$time)
  ait_betadispersion <- vegan::betadisper(as.dist(distance.aitchison), metadat_sel$time)
  
  
  # statistical result 
  statistic_result <- list()
  statistic_result[["bray_distance"]] <- distance.bray
  statistic_result[["aitchison"]] <- distance.aitchison
 
  statistic_result[["bray_adonis_time"]] <- permanova_result_bray_time

  
  statistic_result[["aitchison_adonis_time"]] <- permanova_result_aitchison_time

  statistic_result[["bray_beta_dispersion_time"]] <- bray_betadispersion
  statistic_result[["aitchison_beta_dispersion_time"]] <- ait_betadispersion
  statistic_result[["meta"]] <- metadat_sel
  
  return(statistic_result)
  
}




batch_beta_diversity_statistic <- function(beta_diversity_separate = NA){
  
  beta_statistics_time <- list()
  for (p in names(beta_diversity_separate)){
    print(p)
    beta_statistics_time[[p]] <- distance_produce(beta_diversity = beta_diversity_separate[[p]])
  }
  return(beta_statistics_time)
}


beta_statistics_genera <- batch_beta_diversity_statistic(beta_diversity_separate = genera_separate_profile)
beta_statistics_genera_fraction <- batch_beta_diversity_statistic(beta_diversity_separate = genera_separate_profile_fraction)


extract_r_square <- function(input = NA){
  
  x <- input
  r2 <- x[["bray_adonis_time"]]$R2[1]
  p <- x[["bray_adonis_time"]]$`Pr(>F)`[1]
  
  aitchison_r2 <- x[["aitchison_adonis_time"]]$R2[1]
  aitchison_p <- x[["aitchison_adonis_time"]]$`Pr(>F)`[1]
  
  dispersion_p_tmp <- vegan::permutest(x[["bray_beta_dispersion_time"]])
  dispersion_p <- dispersion_p_tmp$tab$`Pr(>F)`[1]
  
  dispersion_p_aitchison_tmp <- vegan::permutest(x[["aitchison_beta_dispersion_time"]])
  dispersion_p_aitchison <- dispersion_p_aitchison_tmp$tab$`Pr(>F)`[1]
  
  summary_tb <- tibble(bray_r2 = r2, bray_p = p, 
                       aitchison_r2 = aitchison_r2, aitchison_p = aitchison_p,
                       bray_beta_dispersion = dispersion_p, 
                       aitchison_beta_dispersion = dispersion_p_aitchison)
  
  return(summary_tb)
  
}

beta_statistic_result <- tibble(list_name = names(beta_statistics_genera)) %>% rowwise() %>% 
  mutate(r2_p = list(extract_r_square(input = beta_statistics_genera[[list_name]])))

beta_statistic_result_unnest <- beta_statistic_result %>% unnest(r2_p)
saveRDS(beta_statistics_genera, str_c(output_dir, "/beta_diversity/2_beta_statistic_genera.rds"))
write_tsv(beta_statistic_result_unnest, str_c(output_dir, "/beta_diversity/2_beta_statistic_tibble.tsv"))


# fraction 
beta_statistic_result_fraction <- tibble(list_name = names(beta_statistics_genera_fraction)) %>% rowwise() %>% 
  mutate(r2_p = list(extract_r_square(input = beta_statistics_genera_fraction[[list_name]])))

beta_statistic_result_unnest_fraction <- beta_statistic_result_fraction %>% unnest(r2_p)


# fraction 
saveRDS(beta_statistics_genera_fraction, str_c(output_dir, "/beta_diversity/2_beta_statistic_genera_fraction.rds"))
write_tsv(beta_statistic_result_unnest_fraction, str_c(output_dir, "/beta_diversity/2_beta_statistic_tibble_fraction.tsv"))



############################################

# calculate robust aitchison 

############################################

distance_produce_robust_aitchison <- function(beta_diversity = NA){ # 1st list metadata; 2nd list ec_profile
  
  
  metadat_sel <- beta_diversity[[1]]
  
  profile <- beta_diversity[[2]]
  profile <- profile[metadat_sel$sample, ]
  
  
  # distance calculation 
  set.seed(1234)
  distance.bray <- profile %>%
    vegan::vegdist(.,method = 'bray') %>%
    as.matrix() 
  
  set.seed(1234)
  aitchson_input <- profile
  
  
  distance.aitchison <- profile %>%
    vegan::vegdist(.,method = 'robust.aitchison') %>%
    as.matrix() 
  
  
  
  
  # separate adonis for each factor - bray curtis 
  metadat_sel$time <- as.character(metadat_sel$time)
  print(metadat_sel$time)
  permanova_result_bray_time <- vegan::adonis2(distance.bray ~ time, data = metadat_sel)
 
  
  
  # separate adonis for each factor - aitchison 
  permanova_result_aitchison_time <- vegan::adonis2(distance.aitchison ~ time, data = metadat_sel)
 
  
  
  # beta dispersion 
  bray_betadispersion <- vegan::betadisper(as.dist(distance.bray), metadat_sel$time)
  ait_betadispersion <- vegan::betadisper(as.dist(distance.aitchison), metadat_sel$time)
  
  
  # statistical result 
  statistic_result <- list()
  statistic_result[["bray_distance"]] <- distance.bray
  statistic_result[["aitchison"]] <- distance.aitchison
 
  statistic_result[["bray_adonis_time"]] <- permanova_result_bray_time

  
  statistic_result[["aitchison_adonis_time"]] <- permanova_result_aitchison_time

  statistic_result[["bray_beta_dispersion_time"]] <- bray_betadispersion
  statistic_result[["aitchison_beta_dispersion_time"]] <- ait_betadispersion
  statistic_result[["meta"]] <- metadat_sel
  
  return(statistic_result)
  
}

batch_beta_diversity_statistic_v2 <- function(beta_diversity_separate = NA){
  
  beta_statistics_time <- list()
  for (p in names(beta_diversity_separate)){
    print(p)
    beta_statistics_time[[p]] <- distance_produce_robust_aitchison(beta_diversity = beta_diversity_separate[[p]])
  }
  return(beta_statistics_time)
}


beta_statistics_genera_fraction_v2 <- batch_beta_diversity_statistic_v2(beta_diversity_separate = genera_separate_profile_fraction)

beta_statistic_result_fraction_v2 <- tibble(list_name = names(beta_statistics_genera_fraction_v2)) %>% rowwise() %>% 
  mutate(r2_p = list(extract_r_square(input = beta_statistics_genera_fraction_v2[[list_name]])))

beta_statistic_result_unnest_fraction_v2 <- beta_statistic_result_fraction_v2 %>% unnest(r2_p)


# fraction 
saveRDS(beta_statistics_genera_fraction_v2, str_c(output_dir, "/beta_diversity/2_beta_statistic_genera_fraction_v2.rds"))
write_tsv(beta_statistic_result_unnest_fraction_v2, str_c(output_dir, "/beta_diversity/2_beta_statistic_tibble_fraction_v2.tsv"))

save.image(str_c(output_dir, "diversity_genera.RData"))


