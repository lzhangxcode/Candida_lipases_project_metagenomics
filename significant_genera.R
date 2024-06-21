# this script is used to produce the genera level analysis for 20% prevalence at least 0.1% abundance 
# author: lu.zhang
# date: dec.14th,2023 + 15,jan.2024

##########################################################################

# combine the genera phyloseq object 

##########################################################################
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

###########################################################################################

# 20% prevalence filtering 

###########################################################################################

ten_percent_prevalence_filtering_group <- function(phyloseq_objects_input = NA){
  
  phyloseq_objects_target <- list()
  for (all_object in names(phyloseq_objects_input)){
    phyloseq_objects_target[[all_object]] <- list()
    #print(all_object)
    abx <- c("ABX", "No_ABX")
    diet <- c("Standard_diet", "High_fat")
    group <- c("Mutant", "WT")
    #add new column
    sample_data(phyloseq_objects_input[[all_object]])$anti_diet <-  paste(sample_data(phyloseq_objects_input[[all_object]])$anti, sample_data(phyloseq_objects_input[[all_object]])$diet, sample_data(phyloseq_objects_input[[all_object]])$group,
                                                                          sep = "_")
    for (anti in unique(abx)){
      for (diet_target in unique(diet)){
        for (group_target in unique(group)){
          target_label <- paste(anti, diet_target, group_target, sep = "_")
          
          phyloseq_objects_antibiotic_sd <- do.call("subset_samples", list(quote(phyloseq_objects_input[[all_object]]), substitute(anti_diet == target_label))) 
          
          
          
          
          nsamples <- nsamples(phyloseq_objects_antibiotic_sd)
          phyloseq_objects_subset <- phyloseq::filter_taxa(phyloseq_objects_antibiotic_sd, function(x){(sum(x > 0) > nsamples*0.2)},  prune=TRUE) 
          
          phyloseq_objects_target[[all_object]][[target_label]] <- phyloseq_objects_subset
        }
      }
    }
  }
  return(phyloseq_objects_target)
}

phyloseq_objects_ten_prevalence_group <- ten_percent_prevalence_filtering_group(phyloseq_objects_input = phyloseq_objects_genera_combined_from_species)
phyloseq_objects_ten_prevalence_group_ori <- phyloseq_objects_ten_prevalence_group


####################################################################################################

#     prevalence and taxonomic abundance filter 

####################################################################################################

for (k in names(phyloseq_objects_ten_prevalence_group[[1]])){
  print(k)
  
  temp_phyloseq_object <- phyloseq_objects_ten_prevalence_group[[1]][[k]]
  temp_phyloseq_relab <- microbiome::transform(temp_phyloseq_object, trans = "compositional") %>% 
    psmelt()
  total_n <- temp_phyloseq_relab$Sample %>% unique() %>% length()
  temp_keep_species <- temp_phyloseq_relab %>% 
    filter(Abundance > 0.001) %>% 
    group_by(OTU) %>% 
    summarise(n = n()) %>%
    ungroup() %>% 
    mutate(total_sample = total_n) %>% 
    mutate(prevalene = n/total_n) %>% 
    filter(prevalene >= 0.2) 
  
  phyloseq_objects_ten_prevalence_group[[1]][[k]] <- prune_taxa(temp_keep_species$OTU, temp_phyloseq_object)
  
}


phyloseq_objects_ten_prevalence_group

####################################################################################################

#     MetagenomeSeq analysis between day14 vs day0 

####################################################################################################


### metage day0 vs day14 : code checked 
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

  full_features = rbind(Dataset_1, Dataset_2)
  
 
  full_metadata  <- Combined_metadata  %>% .[which(rownames(.) %in% rownames(full_features)),]
  
  
  full_metadata$Group <- factor(full_metadata$group) # Levels: Mutant WT: this code doesnt help 
  
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
  
  ## Creating MR experiment Object
  prof_mr <- newMRexperiment(counts = tran, phenoData = phenoData, featureData = taxaData)
  
  
  prevalence = ncol(prof_mr)*0.1
  
  prof_mr <- filterData(prof_mr, depth = 1000, present = prevalence)
  
  ## Cumulative Sum Scaling Normalisation
  p <- cumNormStatFast(prof_mr)
  tmpMR <- cumNorm(prof_mr, p)
  
  ## Creating the model 
  mod <-  model.matrix(~ time, data = pData(tmpMR))
  
  
  ## FitZig =  zero inflated model ideal for metaphlan data
  results_zeller <-  fitZig(tmpMR, mod)
  
  ## Results
  result = MRtable(results_zeller, number = nrow(results_zeller))
  result = result[result$pvalues < 0.05,]
  return(result)
}



batch_separate_data_for_metagenomeseq_time <- function(input_phyloseq = NA){
  
  # loop different qc & kneaddata method: 
  metagenome_result <- list()
  for (all_object in names(input_phyloseq)){
    abx <- c("ABX", "No_ABX")
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


####################################################################################################

#     MetagenomeSeq analysis between day14 vs day0  no threshold 
 
####################################################################################################

batch_separate_data_for_metagenomeseq_time_full <- function(input_phyloseq = NA){
  
  # loop different qc & kneaddata method: 
  metagenome_result <- list()
  for (all_object in names(input_phyloseq)){
    abx <- c("ABX", "No_ABX")
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

metagenomeseq_difference_abundance_full <- function(Dataset_1 = NA, Dataset_2 = NA, prevCutoff = 0.1, Combined_metadata = NA) {#this is not used prevCutoff
  
  #Dataset_1 = otu_hydrac_crab 
  #Dataset_2 = otu_control_with_crab
  #Combined_metadata =  metatable_sel_biorep
  dataset_name_1 = deparse(substitute(Dataset_1))
  dataset_name_2 = deparse(substitute(Dataset_2))
  # 
  ### Setting the sets
  #full_features = bind_rows(Dataset_1[,2:ncol(Dataset_1)], Dataset_2[,2:ncol(Dataset_2)])
  #full_features[is.na(full_features)] <- 0
  
  full_features = rbind(Dataset_1, Dataset_2)
  
  ## setting the metadata samples that exist in the current comparison
  #Combined_metadata <- Combined_metadata %>%
  #  column_to_rownames(., var = "uniquebiorep")
  full_metadata  <- Combined_metadata  %>% .[which(rownames(.) %in% rownames(full_features)),]
  
  
  #full_metadata$Group = factor(c(rep(dataset_name_1,nrow(Dataset_1)),
  #                                 rep(dataset_name_2,nrow(Dataset_2))))
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
  ## Filtering for 1) Depth > 0.001 2) 10% prevalence
  #prof_mr <- filterData(prof_mr, depth = 0.001, present = prevalence)
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
  #result = result[result$adjPvalues < 0.05,]
  #result = result[result$pvalues < 0.05,]
  return(result)
}

metagenomseq_result_list_time_no_threshold <- batch_separate_data_for_metagenomeseq_time_full(input_phyloseq = phyloseq_objects_ten_prevalence_group)

####################################################################################################

#  write the result to tsv 

####################################################################################################


# write the tsv 

write_tsv_batch <- function(input_list = NA, output_dir = NA, tag = NA){
  
  for (i in names(input_list)){
    for (j in names(input_list[[i]])){
      write_tsv(input_list[[i]][[j]], paste(output_dir, paste0(i, "_", j, "_", tag, ".tsv"), sep = ""))
    }
  }
}


write_tsv_batch(input_list = metagenomseq_result_list_time, output_dir = "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output_strict/sig_genera.timepoint/metagenomeseq/", tag = "metagenomeseq")


####################################################################################################

#  qvalue filtering 

####################################################################################################


#### metagenomeseq: qvalue filtering 


metagenomeseq_result_as_input <- metagenomseq_result_list_time[[1]]
# group0 is timepoint 0, and group1 is timepoint 14. 
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


####################################################################################################

#  overlap between WT and Mutant 

####################################################################################################

input_file <- metagenomseq_result_list_time_q_0.1[[1]]

venn_plot_function <- function(a = NA, b = NA, file = NA, input_file = NA){
  
  a_list <- input_file[[a]]$taxon
  b_list <- input_file[[b]]$taxon
  
  list_venn <- list('WT' = a_list, 'Mutant' = b_list)
  p <- ggvenn::ggvenn(list_venn, c("WT", "Mutant")) + ggtitle(file %>% str_replace_all("_", " "))

  return(p)
  
}



compare_the_overlap_metagenomeseq <- function(input_file = NA){
  
  compar_tibble <- tibble(group1 = names(input_file) %>% str_remove_all("_WT") %>% str_remove_all("_Mutant") %>% unique()) %>% 
    mutate(group2 = group1) %>% 
    mutate(group = group1) %>% 
    mutate(group1 = paste(group1, "_WT", sep = "")) %>% 
    mutate(group2 = paste(group2, "_Mutant", sep = "")) 
  
  compar_tibble
  
  p_final <- compar_tibble %>% 
    rowwise() %>% 
    mutate(plots = list(venn_plot_function(group1, group2, group, input_file)))
  
  
  return(p_final)
  
  
  
}


#venn_list <- metagenomseq_result_list_time_q_0.1_venn 

save_venn_plot <- function(venn_list = NA, qtag = NA, output_dir=NA){
  
  p_combine <- list()
  for (i in names(venn_list)){
    
    venn_list[[i]]$tag <- i 
    p_combine[[i]] <- ggarrange(plotlist = venn_list[[i]]$plots)
    ggsave(paste0(output_dir, i, "_", qtag, ".pdf"), p_combine[[i]], width = 9)
    
  }
  
  
}

metagenomseq_result_list_time_q_0.1_venn <- lapply(metagenomseq_result_list_time_q_0.1, compare_the_overlap_metagenomeseq)
save_venn_plot(metagenomseq_result_list_time_q_0.1_venn, qtag = "0.1", output_dir="/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output_strict/sig_genera.timepoint/vennplot/metag/")

metagenomseq_result_list_time_q_0.2_venn <- lapply(metagenomseq_result_list_time_q_0.2, compare_the_overlap_metagenomeseq)
save_venn_plot(metagenomseq_result_list_time_q_0.2_venn, qtag = "0.2", output_dir="/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output_strict/sig_genera.timepoint/vennplot/metag/")

metagenomseq_result_list_time_q_0.05_venn <- lapply(metagenomseq_result_list_time_q_0.05, compare_the_overlap_metagenomeseq)
save_venn_plot(metagenomseq_result_list_time_q_0.05_venn, qtag = "0.05", output_dir="/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output_strict/sig_genera.timepoint/vennplot/metag/")


metagenomseq_result_list_time_q_0.01_venn <- lapply(metagenomseq_result_list_time_q_0.01, compare_the_overlap_metagenomeseq)
save_venn_plot(metagenomseq_result_list_time_q_0.01_venn, qtag = "0.01", output_dir="/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output_strict/sig_genera.timepoint/vennplot/metag/")



save.image("/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output_strict/sig_genera.timepoint/sig_species.RData")# even though the RData called species.. but it is reality genus, it is a typo here .


