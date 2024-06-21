library(microbiome)
library(phyloseq)
library(tidyverse)
library(rstatix)
library(vegan)



output_dir <- "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/functional_analysis/output/"
############################################
#                 input                    # 
############################################

metatable_paired <- read_tsv("~/Documents/202210_mouse_project/04_diversity/metatable_paired.tsv")

input_dir <- "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/functional_analysis/"
pathway_profile <- read_tsv(paste(input_dir, "KOstat.cpm.txt", sep = "/"))

# D0_3S8 is not successful, so I have to remove this sample D0_3S8 & D14_3S8
metatable_paired_abx <- metatable_paired %>% filter(antibiotic == "ABX") %>% 
  filter(!sample %in% c("D0_3S8", "D14_3S8"))

############################################
#           build a phyloseq object        # 
############################################

build_phyloseq_object <- function(input_pathway = NA, meta = NA){
  
  input_pathway <- input_pathway %>% 
    column_to_rownames(., var = "KO") %>% 
    rename_with(., ~ gsub(".m8", "", .x, fixed = TRUE)) %>% 
    select(all_of(meta$sample))
  
  # remove all zero samples 
  
  rm_samples <- input_pathway %>% colSums() %>% as.data.frame() %>% add_column(sumvalue = .[,1]) %>% filter(sumvalue == 0)
  input_pathway <- input_pathway %>% select(-all_of(rownames(rm_samples)))
  keep_pathway <- rowSums(input_pathway) %>% as.data.frame() %>% add_column(rowsum = .[,1]) %>% filter(rowsum > 0)
  print(nrow(input_pathway))
  print(nrow(keep_pathway))
  
  fake_taxa <- data.frame(taxon = rownames(input_pathway)) %>% as.matrix()
  rownames(fake_taxa) <- rownames(input_pathway)
  
  meta <- meta %>% 
    filter(!sample %in% all_of(rownames(rm_samples))) %>% 
    column_to_rownames(., var = "sample") 
  
  
  ps_pathways <- phyloseq(otu_table(input_pathway, taxa_are_rows = T), tax_table(fake_taxa), sample_data(meta))
  
  return(ps_pathways)
  
}


phyloseq_pathway <- build_phyloseq_object(input_pathway = pathway_profile, meta = metatable_paired_abx)



############################################
#     alpha diversity - comparison         # 
############################################


pathway_alpha_diversity <- microbiome::alpha(phyloseq_pathway, index = "all")
input_alpha <- pathway_alpha_diversity
input_paired <- metatable_paired_abx


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
    filter(diet == diet_choose, antibiotic == antibiotic_choose, group_choose == group) %>% 
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


alpha_diversity_wilcox_paried_table <- alpha_diversity_wilcox_paried %>% unnest(table)

# write the result into tsv file 
write_tsv(alpha_diversity_wilcox_paried_table, str_c(output_dir, "/1_alpha_diversity_wilcoxon_paired.tsv"))
write_tsv(input_alpha %>% rownames_to_column(., var = "samples"), str_c(output_dir, "/1_alpha_diversity.tsv"))
write_tsv(input_paired, str_c(output_dir, "/1_alpha_diversity_paired_metatable.tsv"))


##########################################################################
#     sig.pathways - 10% prevalence filter + abundance filter 0.1%       #
##########################################################################

# filter 10% prevalence + 0.1% abundance filtering 

abundance_filter <- function(input_laststep = NA){
  
  relab_phyloseq <- transform_sample_counts(input_laststep, 
                                            function(x) x / sum(x) ) 
  
  relab_otu_table <- relab_phyloseq %>% otu_table() %>% as.data.frame()
  
  relab_otu_table$max <- -1
  for (i in 1:nrow(relab_otu_table)){
    
    max_value <- max(relab_otu_table[i,])
    relab_otu_table$max[i] <- max_value
  }
  relab_otu_table_0.1 <- relab_otu_table %>% filter(max >= 0.001) %>% select(-max)  
  
  sig_taxon <- prune_taxa(rownames(relab_otu_table_0.1), input_laststep)
  
  
  return(sig_taxon)
  
  
}


ten_percent_prevalence_filtering_group <- function(phyloseq_objects_input = NA){
  
  phyloseq_objects_target <- list()
  for (all_object in names(phyloseq_objects_input)){
    phyloseq_objects_target[[all_object]] <- list()
    #print(all_object)
    abx <- c("ABX")
    diet <- c("Standard_diet", "High_fat")
    group <- c("Mutant", "WT")
    #add new column
    sample_data(phyloseq_objects_input[[all_object]])$anti_diet <-  paste(sample_data(phyloseq_objects_input[[all_object]])$anti, sample_data(phyloseq_objects_input[[all_object]])$diet, sample_data(phyloseq_objects_input[[all_object]])$group,
                                                                          sep = "_")
    for (anti in unique(abx)){
      for (diet_target in unique(diet)){
        for (group_target in unique(group)){
          target_label <- paste(anti, diet_target, group_target, sep = "_")
          #print(target_label)
          
          #filter by the new column : code based on here https://github.com/joey711/phyloseq/issues/1522
          phyloseq_objects_antibiotic_sd <- do.call("subset_samples", list(quote(phyloseq_objects_input[[all_object]]), substitute(anti_diet == target_label))) 
          
          
          nsamples <- nsamples(phyloseq_objects_antibiotic_sd)
          phyloseq_objects_subset <- phyloseq::filter_taxa(phyloseq_objects_antibiotic_sd, function(x){(sum(x > 0) > nsamples*0.1)},  prune=TRUE) # paste from here: https://github.com/joey711/phyloseq/issues/1555
          
          
          phyloseq_objects_subset_v2 <- abundance_filter(phyloseq_objects_subset)
          
          phyloseq_objects_target[[all_object]][[target_label]] <- phyloseq_objects_subset_v2
        }
      }
    }
  }
  return(phyloseq_objects_target)
}


sig_taxon_list <- list("diamond" = phyloseq_pathway)
phyloseq_objects_ten_prevalence_group <- ten_percent_prevalence_filtering_group(phyloseq_objects_input = sig_taxon_list)

saveRDS(phyloseq_objects_ten_prevalence_group, str_c(output_dir, "/3_phyloseq_objects_ten_prevalence_group.rds"))


############################################
#     sig.pathways - metagenomeSeq         # 
############################################
# code checked 


library(metagenomeSeq)
separate_dataset_for_metagenomesq_time <- function(phyloseq_target = NA){
  
  phyloseq_sel_wt <- subset_samples(phyloseq_target, time == 0) 
  phyloseq_sel_mutant <- subset_samples(phyloseq_target, time == 14) 
  metadata_sel <- phyloseq_target %>% sample_data() 
  dataset_1 <- as(otu_table(phyloseq_sel_wt), "matrix")
  dataset_2 <- as(otu_table(phyloseq_sel_mutant), "matrix") # tips from here: https://github.com/joey711/phyloseq/issues/613
  if(taxa_are_rows(phyloseq_sel_wt)){dataset_1 <- t(dataset_1)}
  if(taxa_are_rows(phyloseq_sel_mutant)){dataset_2 <- t(dataset_2)}
  
  # Coerce to data.frame
  dataset_1_df = as.data.frame(dataset_1)
  dataset_2_df = as.data.frame(dataset_2)
  
  result <- list(dataset_1_df, dataset_2_df, metadata_sel)
  return(result)
  
  
}


metagenomeseq_difference_abundance <- function(Dataset_1 = NA, Dataset_2 = NA, prevCutoff = 0.1, Combined_metadata = NA) {#this is not used prevCutoff
  
  
  dataset_name_1 = deparse(substitute(Dataset_1))
  dataset_name_2 = deparse(substitute(Dataset_2))
  # 
  ### Setting the sets
  
  full_features = rbind(Dataset_1, Dataset_2)
  
  ## setting the metadata samples that exist in the current comparison
 
  full_metadata  <- Combined_metadata  %>% .[which(rownames(.) %in% rownames(full_features)),]
  
  
  
  full_metadata$Group <- factor(full_metadata$group) 
  
  ### Order them 
  full_metadata = full_metadata[order(rownames(full_metadata)),] #%>% .[complete.cases(.), ]
  full_features = full_features[order(rownames(full_features)),] %>% .[which(rownames(.) %in% rownames(full_metadata)),]
  
  
  ### transoposed data frame of features
  tran = t(full_features) %>% as.data.frame(.)
  
  ### Taxonomy dataframe , item needed for MRExperiment object
  taxaData = data.frame(colnames(full_features))
  rownames(taxaData) = names(full_features)
  colnames(taxaData) = "OTU"
  taxaData <- AnnotatedDataFrame(taxaData)
  
  ### Metadata in annotated daa frame format with teh same order and rownames as the full features
  phenoData <- AnnotatedDataFrame(full_metadata)
  
  ##3 Creating MR experiment Object
  prof_mr <- newMRexperiment(counts = tran, phenoData = phenoData, featureData = taxaData)
  
  
  prevalence = ncol(prof_mr)*0.1

  prof_mr <- filterData(prof_mr, depth = 1000, present = prevalence)
  
  ## Cumulative Sum Scaling Normalisation
  p <- cumNormStatFast(prof_mr)
  tmpMR <- cumNorm(prof_mr, p)
  

  mod <-  model.matrix(~ time, data = pData(tmpMR))
  
  
  ## FitZig =  zero inflated model ideal for metaphlan data
  results_zeller <-  fitZig(tmpMR, mod)
  
  ## Results
  result = MRtable(results_zeller, number = nrow(results_zeller))
  #result = result[result$adjPvalues < 0.05,]
  result = result[result$pvalues < 0.05,]
  return(result)
}


batch_separate_data_for_metagenomeseq_time <- function(input_phyloseq = NA){
  
  # loop different qc & kneaddata method: 
  metagenome_result <- list()
  for (all_object in names(input_phyloseq)){
    abx <- c("ABX")
    diet <- c("Standard_diet", "High_fat")
    group <- c("Mutant", "WT")
    for (anti in unique(abx)){
      for (diet_target in unique(diet)){
        for (group_target in unique(group)){
          #print(group_target)
          
          target_label <- paste(anti, diet_target, group_target, sep = "_")
          print(target_label)
          phyloseq_target <- input_phyloseq[[all_object]][[target_label]]
          #print(target_label)
          # convert time point into vector 
          sample_data(phyloseq_target)$time <- as.character(sample_data(phyloseq_target)$time)
          
          dataset_list <- list()
          dataset_list <- separate_dataset_for_metagenomesq_time(phyloseq_target = phyloseq_target)
          metagenome_result[[all_object]][[target_label]] <- metagenomeseq_difference_abundance(Dataset_1 = dataset_list[[1]], Dataset_2 = dataset_list[[2]], Combined_metadata = dataset_list[[3]]) %>% 
            mutate(taxon = rownames(.))
        }
      }
    }
  }
  return(metagenome_result)
  
}



metagenomseq_result_list_time <- batch_separate_data_for_metagenomeseq_time(input_phyloseq = phyloseq_objects_ten_prevalence_group)

# no threshold 

metagenomeseq_difference_abundance_full <- function(Dataset_1 = NA, Dataset_2 = NA, prevCutoff = 0.1, Combined_metadata = NA) {#this is not used prevCutoff
  
 
  dataset_name_1 = deparse(substitute(Dataset_1))
  dataset_name_2 = deparse(substitute(Dataset_2))
  
  full_features = rbind(Dataset_1, Dataset_2)
  

  full_metadata  <- Combined_metadata  %>% .[which(rownames(.) %in% rownames(full_features)),]
  
  
  
  full_metadata$Group <- factor(full_metadata$group) # Levels: Mutant WT
  
  ### Order them 
  full_metadata = full_metadata[order(rownames(full_metadata)),] #%>% .[complete.cases(.), ]
  full_features = full_features[order(rownames(full_features)),] %>% .[which(rownames(.) %in% rownames(full_metadata)),]
  
  
  ### transoposed data frame of features
  tran = t(full_features) %>% as.data.frame(.)
  
  ### Taxonomy dataframe , item needed for MRExperiment object
  taxaData = data.frame(colnames(full_features))
  rownames(taxaData) = names(full_features)
  colnames(taxaData) = "OTU"
  taxaData <- AnnotatedDataFrame(taxaData)
  
  ### Metadata in annotated daa frame format with teh same order and rownames as the full features
  phenoData <- AnnotatedDataFrame(full_metadata)
  
  ##3 Creating MR experiment Object
  prof_mr <- newMRexperiment(counts = tran, phenoData = phenoData, featureData = taxaData)
  
  
  prevalence = ncol(prof_mr)*0.1
 
  prof_mr <- filterData(prof_mr, depth = 1000, present = prevalence)
  
  ## Cumulative Sum Scaling Normalisation
  p <- cumNormStatFast(prof_mr)
  tmpMR <- cumNorm(prof_mr, p)
  
  #NORMFACTOR <- normFactors(prof_mr) # lU ADDING 
  ## Creating the model 
  mod <-  model.matrix(~ time, data = pData(tmpMR))
  
  
  ## FitZig =  zero inflated model ideal for metaphlan data
  results_zeller <-  fitZig(tmpMR, mod)
  
  ## Results
  result = MRtable(results_zeller, number = nrow(results_zeller))

  return(result)
}

batch_separate_data_for_metagenomeseq_time_full <- function(input_phyloseq = NA){
  
  # loop different qc & kneaddata method: 
  metagenome_result <- list()
  for (all_object in names(input_phyloseq)){
    abx <- c("ABX")
    diet <- c("Standard_diet", "High_fat")
    group <- c("Mutant", "WT")
    for (anti in unique(abx)){
      for (diet_target in unique(diet)){
        for (group_target in unique(group)){
          #print(group_target)
          
          target_label <- paste(anti, diet_target, group_target, sep = "_")
          print(target_label)
          phyloseq_target <- input_phyloseq[[all_object]][[target_label]]
          #print(target_label)
          # convert time point into vector 
          sample_data(phyloseq_target)$time <- as.character(sample_data(phyloseq_target)$time)
          
          dataset_list <- list()
          dataset_list <- separate_dataset_for_metagenomesq_time(phyloseq_target = phyloseq_target)
          metagenome_result[[all_object]][[target_label]] <- metagenomeseq_difference_abundance_full(Dataset_1 = dataset_list[[1]], Dataset_2 = dataset_list[[2]], Combined_metadata = dataset_list[[3]]) %>% 
            mutate(taxon = rownames(.))
        }
      }
    }
  }
  return(metagenome_result)
  
}


metagenomseq_result_list_time_no_threshold <- batch_separate_data_for_metagenomeseq_time_full(input_phyloseq = phyloseq_objects_ten_prevalence_group)

############################################
#     sig.pathways - metagenomeSe-write tsv# 
############################################
# code checked 

# write the tsv 

write_tsv_batch <- function(input_list = NA, output_dir = NA, tag = NA){
  
  for (i in names(input_list)){
    for (j in names(input_list[[i]])){
      write_tsv(input_list[[i]][[j]], paste(output_dir, paste0(i, "_", j, "_", tag, ".tsv"), sep = ""))
    }
  }
}


write_tsv_batch(input_list = metagenomseq_result_list_time, output_dir = str_c(output_dir, "/metagenomeSeq/"), tag = "metagenomeseq")


############################################
#     sig.pathways - metagenomeSeq - qvalue# 
############################################

# code checked 

metagenomeseq_q_value_filter_0.1 <- function(metagenomeseq_result_as_input = NA, q = 0.1){
  
  metagenomeseq_result_as_result <- list()
  
  for (i in names(metagenomeseq_result_as_input)){
    
    metagenomeseq_result_q_value <- metagenomeseq_result_as_input[[i]] %>% 
      filter(adjPvalues < q)
    
    metagenomeseq_result_as_result[[i]] <- metagenomeseq_result_q_value
  }
  
  
  return(metagenomeseq_result_as_result)
  
}


metagenomeseq_q_value_filter_0.2 <- function(metagenomeseq_result_as_input = NA, q = 0.2){
  
  metagenomeseq_result_as_result <- list()
  
  for (i in names(metagenomeseq_result_as_input)){
    
    metagenomeseq_result_q_value <- metagenomeseq_result_as_input[[i]] %>% 
      filter(adjPvalues < q)
    
    metagenomeseq_result_as_result[[i]] <- metagenomeseq_result_q_value
  }
  
  
  return(metagenomeseq_result_as_result)
  
}



metagenomeseq_q_value_filter_0.05 <- function(metagenomeseq_result_as_input = NA, q = 0.05){
  
  metagenomeseq_result_as_result <- list()
  
  for (i in names(metagenomeseq_result_as_input)){
    
    metagenomeseq_result_q_value <- metagenomeseq_result_as_input[[i]] %>% 
      filter(adjPvalues < q)
    
    metagenomeseq_result_as_result[[i]] <- metagenomeseq_result_q_value
  }
  
  
  return(metagenomeseq_result_as_result)
  
}


metagenomeseq_q_value_filter_0.01 <- function(metagenomeseq_result_as_input = NA, q = 0.01){
  
  metagenomeseq_result_as_result <- list()
  
  for (i in names(metagenomeseq_result_as_input)){
    
    metagenomeseq_result_q_value <- metagenomeseq_result_as_input[[i]] %>% 
      filter(adjPvalues < q)
    
    metagenomeseq_result_as_result[[i]] <- metagenomeseq_result_q_value
  }
  
  
  return(metagenomeseq_result_as_result)
  
}


metagenomseq_result_list_time_q_0.05 <- lapply(metagenomseq_result_list_time, metagenomeseq_q_value_filter_0.05)
metagenomseq_result_list_time_q_0.1 <- lapply(metagenomseq_result_list_time, metagenomeseq_q_value_filter_0.1)
metagenomseq_result_list_time_q_0.2 <- lapply(metagenomseq_result_list_time, metagenomeseq_q_value_filter_0.2)
metagenomseq_result_list_time_q_0.01 <- lapply(metagenomseq_result_list_time, metagenomeseq_q_value_filter_0.01)


saveRDS(metagenomseq_result_list_time_q_0.05, str_c(output_dir, "/metagenomeSeq/metagenomeseq_result_list_time_q_0.05.rds"))
saveRDS(metagenomseq_result_list_time_q_0.01, str_c(output_dir, "/metagenomeSeq/metagenomeseq_result_list_time_q_0.01.rds"))
saveRDS(metagenomseq_result_list_time_q_0.1, str_c(output_dir, "/metagenomeSeq/metagenomeseq_result_list_time_q_0.1.rds"))
saveRDS(metagenomseq_result_list_time_q_0.2, str_c(output_dir, "/metagenomeSeq/metagenomeseq_result_list_time_q_0.2.rds"))
saveRDS(metagenomseq_result_list_time, str_c(output_dir, "/metagenomeSeq/metagenomeseq_result_list_time.rds"))
saveRDS(metagenomseq_result_list_time_no_threshold, str_c(output_dir, "/metagenomeSeq/metagenomeseq_result_list_time_full.rds"))

############################################
#     sig.pathways - metagenomeSeq - table # 
############################################

# code checked 

compare_the_overlap_metagenomeseq_table <- function(input_file = NA){
  
  compar_tibble <- tibble(group1 = names(input_file) %>% str_remove_all("_WT") %>% str_remove_all("_Mutant") %>% unique()) %>% 
    mutate(group2 = group1) %>% 
    mutate(group = group1) %>% 
    mutate(group1 = paste(group1, "_WT", sep = "")) %>% 
    mutate(group2 = paste(group2, "_Mutant", sep = "")) 
  
  compar_tibble
  
  table_final <- compar_tibble %>% 
    rowwise() %>% 
    mutate(nums = list(compare_DA_directions(group1, group2, group, input_file))) %>% 
    pull(nums) %>%
    bind_rows()
  
  
  return(table_final)
  
  
  
}

#input_file <- metagenomseq_result_list_time_q_0.2[[1]]
compare_DA_directions <- function(group1 = NA, group2 = NA, compa = NA, input_file = NA){
  
  wt_up <- input_file[[group1]] %>% filter(time14 > 0)
  wt_down <- input_file[[group1]] %>% filter(time14 < 0)
  mutant_up <- input_file[[group2]] %>% filter(time14 > 0)
  mutant_down <- input_file[[group2]] %>% filter(time14 < 0)
  
  both_up <- intersect(wt_up$taxon, mutant_up$taxon)
  both_down <- intersect(wt_down$taxon, mutant_down$taxon)
  
  summary_tibble <- tibble(taxon = rep("Pathway", 6), group = rep(c("WT", "MUT", "Both"), 2), direction = rep(c("Up", "Down"), each = 3),
                           num = c(nrow(wt_up), nrow(mutant_up), length(both_up), nrow(wt_down), nrow(mutant_down), length(both_down)),
                           comparison = compa) %>% 
    mutate(antibiotic = str_remove_all(comparison, "_Standard_diet") %>% str_remove_all(., "_High_fat"), .before = taxon) %>% 
    mutate(diet = str_remove_all(comparison, "No_ABX_") %>% str_remove_all(., "ABX_"), .before = taxon) %>% 
    dplyr::select(-comparison) %>% 
    add_column(method = "metagenomeseq")
  
  return(summary_tibble)
  
  
}


metagenomseq_result_list_time_q_0.1_table <- lapply(metagenomseq_result_list_time_q_0.1, compare_the_overlap_metagenomeseq_table)
metagenomseq_result_list_time_q_0.2_table <- lapply(metagenomseq_result_list_time_q_0.2, compare_the_overlap_metagenomeseq_table)
metagenomseq_result_list_time_q_0.05_table <- lapply(metagenomseq_result_list_time_q_0.05, compare_the_overlap_metagenomeseq_table)
metagenomseq_result_list_time_q_0.01_table <- lapply(metagenomseq_result_list_time_q_0.01, compare_the_overlap_metagenomeseq_table)

write_summary_table_batch <- function(input_table = NA, tag = NA, output_dir = "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/functional_analysis/output/metagenomeSeq/sig_pathway_table/"){
  
  for (i in names(input_table)){
    
    write_tsv(input_table[[i]], paste0(output_dir, i, "_", tag, ".tsv"))
    
  }
  
  
}
write_summary_table_batch(metagenomseq_result_list_time_q_0.1_table, tag = "0.1")
write_summary_table_batch(metagenomseq_result_list_time_q_0.2_table, tag = "0.2")
write_summary_table_batch(metagenomseq_result_list_time_q_0.05_table, tag = "0.05")
write_summary_table_batch(metagenomseq_result_list_time_q_0.01_table, tag = "0.01")


save.image(str_c(output_dir, "Diamond_pathway.RData"))








