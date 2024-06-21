# this script is used to visualized the pathway output 

## loaded package 

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(metagenomeSeq)
library(pheatmap)
library(ggVennDiagram)
library(viridis)

input_dir <- "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/functional_analysis/output/"
output_dir <- "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/functional_analysis/output/plots/"
input_paried <- read_tsv(str_c(input_dir, "/1_alpha_diversity_paired_metatable.tsv"))

####################################################################

#                 heatmap for sig.pathways                         # 

####################################################################


# input 

phyloseq_objects_ten_prevalence_group <- readRDS(str_c(input_dir, 
                                                       "/3_phyloseq_objects_ten_prevalence_group.rds"))

sig_pathways_q_0.2 <- readRDS(str_c(input_dir, "metagenomeSeq/metagenomeseq_result_list_time_q_0.2.rds"))


# function1: normalized by metagenomeSeq 

Normalization_metagenome <- function(OTUfile = NA, comp = NA, otudir = NA){
  library(data.table)
  print(comp)
  # otu_table(phyloseq_objects_ten_prevalence_group[['diamond']][['No_ABX_High_fat_Mutant']])
  otu_table <- phyloseq::otu_table(OTUfile[[comp]]) %>% as.data.frame()
  
  normfile <-  otu_table %>%
    metagenomeSeq::newMRexperiment(.) %>%       ## creating MR object                      
    metagenomeSeq::cumNorm(., p = metagenomeSeq::cumNormStatFast(.)) %>% ## normalization  
    metagenomeSeq::MRcounts (., norm = TRUE, log = TRUE) %>% # produce counts 
    as_tibble(.) %>%
    add_column(OTUID = rownames(otu_table(OTUfile[[comp]]))) %>%
    as.data.table() 
  write_csv(normfile, paste(otudir, comp, "_normalized_otu_table.csv", sep = ""))
  return(normfile)
}

# function2: batch run 

Run_batch <- function(phyloseq_object = NA){
  
  compar_tibble <- tibble(compar = names(phyloseq_object))
  
  compar_tibble <- compar_tibble %>% 
    rowwise() %>% 
    mutate(normalized_table = list(Normalization_metagenome(OTUfile = phyloseq_objects_ten_prevalence_group[["diamond"]],
                                                            comp = compar, 
                                                            otudir = str_c(output_dir, "/metagenomeSeq_heatmap/"))))
  
  return(compar_tibble)
}


normalized_css_pathway <- Run_batch(phyloseq_object = phyloseq_objects_ten_prevalence_group[["diamond"]])



# function3: only select the sig.pathways
# function4: save the sig.pathway pvlaue and effect size indicates the direction 


extract_sig_pathways <- function(normalized_pathways = NA, sig_pathways = NA, comp = NA){
  
  norm_all <- normalized_pathways$normalized_table[[which(normalized_pathways$compar == comp)]]
  
  sig_path <- sig_pathways[[comp]]
  print(nrow(sig_path))
  norm_sel <- norm_all %>% 
    filter(OTUID %in% rownames(sig_path))
  
  print(nrow(norm_sel))
  
  return(norm_sel)
  
}

extract_pvalue <- function(sig_pathways = NA, comp = NA){
  
  sig_path <- sig_pathways[[comp]]
  
  return(sig_path)
  
}


batch_extract <- function(normalized_pathways = NA, sig_pathways = NA){
  
  
  normalized_path_sel <- normalized_pathways %>% 
    rowwise() %>% 
    mutate(normalized_table_sel = list(extract_sig_pathways(normalized_pathways = normalized_pathways,
                                                            sig_pathways = sig_pathways_q_0.2[["diamond"]],
                                                            comp = compar))) %>%
    mutate(sig_pathways = list(extract_pvalue(sig_pathways = sig_pathways_q_0.2[["diamond"]],
                                              comp = compar)))
  
  return(normalized_path_sel)
  
}



normalized_css_pathway_sel <- batch_extract(normalized_pathways = normalized_css_pathway,
                                            sig_pathways = sig_pathways_q_0.2[["diamond"]])




# function5: draw plots based on the normalized metagenomSeq 

pheatmap_f <- function(otu = NA, meta = NA, row_anno = NA){
  
  
  if (nrow(otu) == 0){
    p <- NULL
  }else{
    
    
    # meta 
    
    meta_sel <-  meta%>% 
      filter(sample %in% colnames(otu)) %>% 
      select(sample, time) %>% 
      column_to_rownames(., var = "sample") %>% 
      mutate(time = as.character(time)) %>% 
      select(time)
    
    print(meta_sel)
    
    #newCols <- colorRampPalette(grDevices::rainbow(length(unique(meta_sel$time))))
    #annoCol <- newCols(length(unique(meta_sel$time)))
    
    # "increase" = "#ef3647", "decrease" = "#298dc4"
    # 
    
    # rowanno 
    
    annotation_row <- row_anno %>% 
      mutate(Direction = ifelse(time14 > 0, "Increase","Decrease"),
             Pvalues = cut(pvalues, breaks = c(0.05, 0.01, 0)),
             AdjustedP = cut(adjPvalues, breaks = c(0.2, 0.1, 0.05, 0.01, 0))) %>% 
      select(Direction, Pvalues, AdjustedP)
    
    
    # draw plots 
    p <- pheatmap::pheatmap(otu %>% column_to_rownames(., var = "OTUID"), 
                            annotation_col = meta_sel,
                            annotation_colors = list(time = c("0" = "#fbb4ae", "14" = "#b3cde3"),
                                                     Direction = c("Increase" = "#ef3647", "Decrease" = "#298dc4"),
                                                     Pvalues = c("(0,0.01]" = "#fec44f", "(0.01,0.05]" = "#fff7bc"),
                                                     AdjustedP = c("(0,0.01]" = "#54278f", "(0.01,0.05]" = "#756bb1", "(0.05,0.1]" = "#9e9ac8", "(0.1,0.2]" = "#cbc9e2")),
                            annotation_row = annotation_row,
                            show_rownames = TRUE,
                            cluster_cols = FALSE,
                            scale = "row",
                            fontsize_col = 5,
                            fontsize_row = 8,
                            cellheight = 12, cellwidth = 6)#,
    #labels_row = as.expression(newnames))
    
    return(p)
  }
  
  
}


batch_pheatmap <- function(normalized_sel = NA){
  
  normalized_sel_p <- normalized_sel %>% 
    rowwise() %>% 
    mutate(plots = list(pheatmap_f(otu = normalized_table_sel, 
                                   meta = input_paried, 
                                   row_anno = sig_pathways)))
  
  
  
}


normalized_css_pathway_sel_plots <- batch_pheatmap(normalized_sel = normalized_css_pathway_sel)
normalized_css_pathway_sel_plots$path <- str_c(output_dir, "metagenomeSeq_heatmap/")

saveRDS(normalized_css_pathway_sel_plots, str_c(output_dir, "/metagenomeSeq_heatmap/normalized_css_pathway_sel_plots.rds"))

batch_savepdf <- function(input = NA, output_dir = NA, filename = NA, otu = NA){
  
  if (nrow(otu) > 150){
    pdf(str_c(output_dir, filename, ".pdf"), width = 20, height = 40)
    print(input)
    dev.off()
  }else{
    pdf(str_c(output_dir, filename, ".pdf"), width = 20*0.6, height = 40*0.6)
    print(input)
    dev.off()
  }
  
  
}

normalized_css_pathway_sel_plots %>% 
  rowwise() %>% 
  mutate(list(batch_savepdf(input = plots,
                            output_dir = path,
                            filename = compar,
                            otu = normalized_table_sel)))


####################################################################

#         KEGG lipid related pathways                              # 

####################################################################

# 1. I can check one by one whether the lipid related pathways are in the sig.pathways list

# 2. a barplot showing the pathways in WT or MUTANT 

# 1.3 Lipid metabolism
# 00061 MFatty acid biosynthesis
# 00062 MFatty acid elongation
# 00071 M NFatty acid degradation
# 00073 Cutin, suberine and wax biosynthesis
# 00100 M NSteroid biosynthesis
# 00120 M NPrimary bile acid biosynthesis
# 00121 Secondary bile acid biosynthesis
# 00140 M NSteroid hormone biosynthesis
# 00561 M Glycerolipid metabolism
# 00564 M Glycerophospholipid metabolism
# 00565 M Ether lipid metabolism
# 00600 M N Sphingolipid metabolism
# 00590 Arachidonic acid metabolism
# 00591 Linoleic acid metabolism
# 00592 Malpha-Linolenic acid metabolism
# 01040 M N Biosynthesis of unsaturated fatty acids

# I can use relative abundance to compare the lipid related pathways between WT and MUTANT 
lipid_related_pathways <- tribble(~ ko, ~pathway,
                                  "ko00061", "Fatty acid biosynthesis",
                                  "ko00062", "Fatty acid elongation",
                                  "ko00071", "Fatty acid degradation",
                                  "ko00073", "Cutin, suberine and wax biosynthesis",
                                  "ko00100", "Steroid biosynthesis",
                                  "ko00120", "Primary bile acid biosynthesis",
                                  "ko00121", "Secondary bile acid biosynthesis",
                                  "ko00140", "Steroid hormone biosynthesis",
                                  "ko00561", "Glycerolipid metabolism",
                                  "ko00564", "Glycerophospholipid metabolism",
                                  "ko00565", "Ether lipid metabolism",
                                  "ko00600", "Sphingolipid metabolism",
                                  "ko00590", "Arachidonic acid metabolism",
                                  "ko00591", "Linoleic acid metabolism",
                                  "ko00592", "Malpha-Linolenic acid metabolism",
                                  "ko01040", "Biosynthesis of unsaturated fatty acids")


# 3. pick up all the lipid-related pathways 


# code checked 
metagenomseq_result_list_time <- readRDS(str_c(input_dir, "/metagenomeSeq/metagenomeseq_result_list_time_full.rds"))


select_lipid_related_path <- function(input = NA, target_lipid = NA){
  
  input_sel <- input %>% 
    filter(taxon %in% target_lipid) 
  return(input_sel)
  
}


select_lipid_related_path <- function(input_list = NA, target_lipid) {
  combined_table <- data.frame()  # Create an empty data frame to store the combined tables
  
  # Loop through each table in the input list
  for (i in seq_along(input_list)) {
    table <- input_list[[i]]  # Get the table from the list
    
    # Create a new column for the list name and assign the corresponding list name to all rows in the table
    table$list_name <- names(input_list)[i]
    
    # Filter the table based on the target lipid
    table_sel <- table %>% 
      dplyr::filter(taxon %in% target_lipid)
    
    # Combine the filtered table with the combined table using rbind()
    combined_table <- rbind(combined_table, table_sel)
  }
  
  return(combined_table)
}



# 4. combine all comparison lipid-related pathways together 

lipid_path <- select_lipid_related_path(input_list = metagenomseq_result_list_time[[1]], 
                                        target_lipid = str_c(lipid_related_pathways$pathway, lipid_related_pathways$ko, sep = ":"))


# 5. ABX highfat as one group et al 


lipid_path_long <- lipid_path %>% 
  select(taxon, time14, pvalues, adjPvalues, list_name)


produce_combined_table_for_each_group <- function(input_long = lipid_path_long, comp = NA){
  
  
  # select the target group 
  lipid_path_long_sel <- input_long %>% 
    mutate(group = str_remove_all(list_name, "_WT|_Mutant")) %>% 
    filter(group == comp) %>% 
    rename(effect_size = time14) %>% 
    mutate(mark = ifelse(adjPvalues < 0.2, "#", "")) %>% 
    mutate(list_name = ifelse(grepl("Mutant", list_name), "Mutant", "WT")) %>% 
    complete(taxon, list_name) %>% 
    arrange(taxon)
  
  return(lipid_path_long_sel)
  
}

produce_bar_plots <- function(long_input = NA, comp = NA, out = NA){
  
  comp_title <- str_replace_all(comp, "_", " ")
  long_input <- long_input %>% rename(Group = list_name)
  
  lipid_plots <- ggplot(long_input, aes(x = taxon, y = effect_size, group = Group, fill = Group, label = mark)) +
    geom_bar(stat = "identity", position = "dodge", show.legend = TRUE) + # position dodge for >= 2 groups x pos
    geom_text(aes(x = taxon, y = sign(effect_size) * 0.1, group = Group), position = position_dodge(width = 1), size = 5, color = "black", fontface = "bold") + 
    #geom_text(aes(hjust = 0.5 - sign(effect_size)/1.5), size=7) +
    coord_flip() +
    theme_classic2() + 
    theme(
      #plot.title = element_text(size = 36, hjust = 0.5),
      #text = element_text(size = 30),
      axis.title.x = element_text(colour='black', size=15),
      axis.title.y = element_text(colour='black', size=20),
      axis.text.y = element_text(colour='black', size=12),
      legend.text = element_text(size=8),
      legend.title = element_text(size=12),
      legend.position = 'right') +
    scale_fill_manual(values = c("#ff006645", "#00000045")) +
    labs(title = comp_title) + xlab(str_c("Lipid related pathways")) + ylab("Effect size")
  
  ggsave(str_c(out, "/", comp, "_lipid_barplot.pdf"), lipid_plots)
  
  return(lipid_plots)
  
  
}

compar_tibble <- tibble(comp = names(sig_pathways_q_0.2[[1]]))


compar_tibble_lipid <- compar_tibble %>% mutate(comp = str_remove_all(comp, "_WT|_Mutant")) %>% 
  distinct() %>% rowwise() %>% 
  mutate(combined_table = list(produce_combined_table_for_each_group(comp = comp))) %>% 
  mutate(output = str_c(output_dir, "/metagenomeSeq_lipid_barplot/")) %>% 
  mutate(combined_plot = list(produce_bar_plots(long_input = combined_table, comp = comp, out = output)))



save.image(file = str_c(output_dir, "/pathway_vis.RData"))

####################################################################

#         KEGG amino acids related pathways                        # 

####################################################################

output_dir <- "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/functional_analysis/output/plots/"
load(file = str_c(output_dir, "/pathway_vis.RData"))

# exceute again the previous barplot code, to change it to black color 

#https://www.genome.jp/brite/query=01212&htext=br08901.keg&option=-a&node_proc=br08901_org&proc_enabled=map&panel=collapse
#00250 Alanine, aspartate and glutamate metabolism
#00260 Glycine, serine and threonine metabolism
#00270 Cysteine and methionine metabolism
#00280 Valine, leucine and isoleucine degradation
#00290 Valine, leucine and isoleucine biosynthesis
#00300 Lysine biosynthesis
#00310 Lysine degradation
#00220 Arginine biosynthesis
#00330 Arginine and proline metabolism
#00340 Histidine metabolism
#00350 Tyrosine metabolism
#00360 Phenylalanine metabolism
#00380 Tryptophan metabolism
#00400 Phenylalanine, tyrosine and tryptophan biosynthesis

amino_acids_pathways <- tribble(~ ko, ~pathway,
                                "ko00250", "Alanine, aspartate and glutamate metabolism",
                                "ko00260", "Glycine, serine and threonine metabolism",
                                "ko00270", "Cysteine and methionine metabolism",
                                "ko00280", "Valine, leucine and isoleucine degradation",
                                "ko00290", "Valine, leucine and isoleucine biosynthesis",
                                "ko00300", "Lysine biosynthesis",
                                "ko00310", "Lysine degradation",
                                "ko00220", "Arginine biosynthesis",
                                "ko00330", "Arginine and proline metabolism",
                                "ko00340", "Histidine metabolism",
                                "ko00350", "Tyrosine metabolism",
                                "ko00360", "Phenylalanine metabolism",
                                "ko00380", "Tryptophan metabolism",
                                "ko00400", "Phenylalanine, tyrosine and tryptophan biosynthesis")

aminoacid_path <- select_lipid_related_path(input_list = metagenomseq_result_list_time[[1]], 
                                        target_lipid = str_c(amino_acids_pathways$pathway, amino_acids_pathways$ko, sep = ":"))


aminoacid_path_long <- aminoacid_path %>% 
  select(taxon, time14, pvalues, adjPvalues, list_name)


save.image(file = str_c(output_dir, "/pathway_vis_add_amino_acid.RData"))


####################################################################

#                 heatmap for sig.pathways                         # 
#                     revised version 
####################################################################

# remove the false positive pathways: human/disease related or non-bacteria pathways and also one antibiotic produce pathways 

output_dir <- "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/functional_analysis/output/plots/"
load(file = str_c(output_dir, "/pathway_vis_add_amino_acid.RData"))

### remove false positive pathways manually 
### sig_pathways_q_0.2[["diamond"]]
sig_pathways_q_0.2_manual <- sig_pathways_q_0.2

## ABX_standard_diet_mutant 
nrow(sig_pathways_q_0.2[["diamond"]]$ABX_Standard_diet_Mutant) #48 
abx_std_mutant_removed <- c("Cell cycle - yeast:ko04111", "Non-alcoholic fatty liver disease (NAFLD):ko04932", "Cardiac muscle contraction:ko04260", "Type I diabetes mellitus:ko04940", "Hepatocellular carcinoma:ko05225", "Human T-cell leukemia virus 1 infection:ko05166", "Parkinson disease:ko05012", "Meiosis - yeast:ko04113", "Pathways in cancer:ko05200", "Chemical carcinogenesis:ko05204", "Alzheimer disease:ko05010", "Huntington disease:ko05016", "Protein processing in endoplasmic reticulum:ko04141", "Central carbon metabolism in cancer:ko05230", "Mitophagy - yeast:ko04139", "MAPK signaling pathway - yeast:ko04011", "Insulin resistance:ko04931", "Glucagon signaling pathway:ko04922")

sig_pathways_q_0.2_manual[["diamond"]]$ABX_Standard_diet_Mutant <- sig_pathways_q_0.2[["diamond"]]$ABX_Standard_diet_Mutant %>% filter(!taxon %in% abx_std_mutant_removed)

View(sig_pathways_q_0.2_manual[["diamond"]]$ABX_Standard_diet_Mutant)
nrow(sig_pathways_q_0.2_manual[["diamond"]]$ABX_Standard_diet_Mutant) #32 

## removed pathways are cancer/disease related or eukaryotes specific 
## "ko04922" Glucagon signaling pathway - Homo sapiens (human) 
## endoplasmic reticulum is in the eukarytoes. 


## ABX_Standard_diet_WT
nrow(sig_pathways_q_0.2[["diamond"]]$ABX_Standard_diet_WT) #0 


## ABX_High_fat_Mutant
nrow(sig_pathways_q_0.2[["diamond"]]$ABX_High_fat_Mutant) #54 

abx_hf_mutant_removed <- c("Leukocyte transendothelial migration:ko04670", "Endometrial cancer:ko05213", "Adherens junction:ko04520", "Hippo signaling pathway:ko04390", "Gastric cancer:ko05226", "Parkinson disease:ko05012", "Retrograde endocannabinoid signaling:ko04723", "Taste transduction:ko04742", "Retrograde endocannabinoid signaling:ko04723", "Arrhythmogenic right ventricular cardiomyopathy (ARVC):ko05412", "Cardiac muscle contraction:ko04260", "Non-alcoholic fatty liver disease (NAFLD):ko04932", "Huntington disease:ko05016", "Meiosis - yeast:ko04113", "Thermogenesis:ko04714", "Alzheimer disease:ko05010", "Central carbon metabolism in cancer:ko05230", "Peroxisome:ko04146", "Streptomycin biosynthesis:ko00521", "Herpes simplex virus 1 infection:ko05168", "Herpes simplex virus 1 infection:ko05168")
abx_hf_mutant_removed <- c(abx_hf_mutant_removed, abx_std_mutant_removed)

## https://www.genome.jp/kegg-bin/show_pathway?ko04520: multiple cells, not bacteria 
## Hippo signaling is an evolutionarily conserved signaling pathway that controls organ size from flies to humans.
## https://www.genome.jp/dbget-bin/www_bget?ko04742 Taste transduction: Five basic tastes are recognized by humans and most other animals
## https://www.genome.jp/dbget-bin/www_bget?ko05168 Herpes simplex virus 1(HSV-1) is a common human pathogen, which initially infects orofacial mucosal surfaces. The virus replicates in epithelial cells at these sites. i don't think bacteria has orofacial mucosal surfaces. 
## https://www.genome.jp/dbget-bin/www_bget?ko04723 Endogenous cannabinoids (endocannabinoids) serve as retrograde messengers at synapses in various regions of the brain. I don't think bacteria has brain.
## Parkinson disease:ko05012 
## Metabolism; Xenobiotics biodegradation and metabolism: https://www.genome.jp/dbget-bin/www_bget?ko00642 : not very strong evidence it is not bacteria 
## Arrhythmogenic right ventricular cardiomyopathy (ARVC):ko05412 human disease 
## Thermogenesis is essential for warm-blooded animals, https://www.genome.jp/dbget-bin/www_bget?ko04714 
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8326754/ they include ko04068	FoxO signaling pathway 
## P450 enzymes have been identified in all kingdoms of life: animals, plants, fungi, protists, bacteria, and archaea, as well as in viruses.[4] 
## bacteria doesn't have nucleus https://www.genome.jp/dbget-bin/www_bget?ko03013   RNA transport from the nucleus to the cytoplasm is fundamental for gene expression. 
## https://www.genome.jp/dbget-bin/www_bget?ko00860 it can be in bacteria : search bacteria 
## The phosphoenolpyruvate (PEP)-dependent phosphotransferase system (PTS) is a major mechanism used by bacteria for uptake of carbohydrates https://www.genome.jp/dbget-bin/www_bget?ko02060
## A peroxisomes a membrane-bound organelle, a type of microbody, found in the cytoplasm of virtually all eukaryotic cells.[2][3] 
## Streptomycin is an aminocyclitol-aminoglycoside antibiotic produced by Streptomyces griseus. Streptomycin consists of aminocyclitol (streptidine), 6-deoxyhexose (streptose), and N-methyl-L-glucosamine moieties, which are formed by independent biosynthetic pathways: excluded because of antibiotic
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7244492/ in this paper it is mentioned in microbiota Synthesis and degradation of ketone bodies:ko00072
## Certain antibiotics, such as tetracyclines, prevent the aminoacyl-tRNA from binding to the ribosomal subunit in prokaryotes. 


sig_pathways_q_0.2_manual[["diamond"]]$ABX_Standard_diet_Mutant <- sig_pathways_q_0.2[["diamond"]]$ABX_Standard_diet_Mutant %>% filter(!taxon %in% abx_hf_mutant_removed)
nrow(sig_pathways_q_0.2_manual[["diamond"]]$ABX_Standard_diet_Mutant) #31 
sig_pathways_q_0.2_manual[["diamond"]]$ABX_High_fat_Mutant <- sig_pathways_q_0.2[["diamond"]]$ABX_High_fat_Mutant %>% filter(!taxon %in% abx_hf_mutant_removed)
nrow(sig_pathways_q_0.2_manual[["diamond"]]$ABX_High_fat_Mutant) #35 


## ABX_High_fat_WT
# Endocytosis is a fundamental process of membrane-trafficking in eukaryotes, but has not been known to occur in bacteria or archaea.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3038067/ Endocytosis:ko04144
# Shear stress represents the frictional force that the flow of blood exerts at the endothelial surface of the vessel wall and plays a central role in vascular biology and contributes to the progress of atherosclerosis. https://www.genome.jp/entry/ko05418
# d-Alanine (D-Ala) metabolism is essential for bacterial peptidoglycan biosynthesis., it should be included D-Alanine metabolism:ko00473
# https://en.wikipedia.org/wiki/Homologous_recombination#In_bacteria it can be in bacteria 

abx_high_fat_wt <- c("Endocytosis:ko04144", "Autophagy - yeast:ko04138", "Fluid shear stress and atherosclerosis:ko05418")
abx_high_fat_wt <- c(abx_high_fat_wt, abx_hf_mutant_removed)


nrow(sig_pathways_q_0.2[["diamond"]]$ABX_High_fat_WT) #31 
sig_pathways_q_0.2_manual[["diamond"]]$ABX_High_fat_WT <- sig_pathways_q_0.2[["diamond"]]$ABX_High_fat_WT %>% filter(!taxon %in% c(abx_high_fat_wt))

nrow(sig_pathways_q_0.2_manual[["diamond"]]$ABX_High_fat_WT) #12 


sig_pathways_q_0.2_manual[["diamond"]]$ABX_High_fat_Mutant <- sig_pathways_q_0.2[["diamond"]]$ABX_High_fat_Mutant %>% filter(!taxon %in% abx_high_fat_wt)
nrow(sig_pathways_q_0.2_manual[["diamond"]]$ABX_High_fat_Mutant) #35 

sig_pathways_q_0.2_manual[["diamond"]]$ABX_Standard_diet_Mutant <- sig_pathways_q_0.2[["diamond"]]$ABX_Standard_diet_Mutant %>% filter(!taxon %in% abx_high_fat_wt)
nrow(sig_pathways_q_0.2_manual[["diamond"]]$ABX_Standard_diet_Mutant) #30




batch_extract_change <- function(normalized_pathways = NA, sig_pathways = NA){
  
  
  normalized_path_sel <- normalized_pathways %>% 
    rowwise() %>% 
    mutate(normalized_table_sel = list(extract_sig_pathways(normalized_pathways = normalized_pathways,
                                                            sig_pathways = sig_pathways,
                                                            comp = compar))) %>%
    mutate(sig_pathways = list(extract_pvalue(sig_pathways = sig_pathways,
                                              comp = compar)))
  
  return(normalized_path_sel)
  
}



normalized_css_pathway_sel_manual <- batch_extract_change(normalized_pathways = normalized_css_pathway,
                                            sig_pathways = sig_pathways_q_0.2_manual[["diamond"]])



normalized_css_pathway_sel_plots_manual <- batch_pheatmap(normalized_sel = normalized_css_pathway_sel_manual)
normalized_css_pathway_sel_plots_manual$path <- str_c(output_dir, "metagenomeSeq_heatmap_change/")

saveRDS(normalized_css_pathway_sel_plots_manual, str_c(output_dir, "/metagenomeSeq_heatmap_change/normalized_css_pathway_sel_plots_manual.rds"))

batch_savepdf <- function(input = NA, output_dir = NA, filename = NA, otu = NA){
  
  if (nrow(otu) > 150){
    pdf(str_c(output_dir, filename, ".pdf"), width = 20, height = 40)
    print(input)
    dev.off()
  }else{
    pdf(str_c(output_dir, filename, ".pdf"), width = 20*0.6, height = 40*0.6)
    print(input)
    dev.off()
  }
  
  
}

normalized_css_pathway_sel_plots_manual %>% 
  rowwise() %>% 
  mutate(list(batch_savepdf(input = plots,
                            output_dir = path,
                            filename = compar,
                            otu = normalized_table_sel)))



### venn plot of the sig.diff patways 

metagenomseq_result_list_time_q_0.2_venn_manual <- lapply(sig_pathways_q_0.2_manual, compare_the_overlap_metagenomeseq)
save_venn_plot(metagenomseq_result_list_time_q_0.2_venn_manual, qtag = "0.2", 
               output_dir=str_c(output_dir, "/metagenomeSeq_venn_manual/"))



save.image(file = str_c(output_dir, "/pathway_vis_manual_remove_fp.RData"))


####################################################################

#                 heatmap for sig.pathways                         # 
#                 revised version 
####################################################################

output_dir <- "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/functional_analysis/output/plots/"
load(str_c(output_dir, "/pathway_vis_manual_remove_fp.RData"))

# check the overlap pathways between wt and mutant in high-fat diet 
overlap_pathways <- intersect(normalized_css_pathway_sel_plots_manual$sig_pathways[[3]]$taxon, normalized_css_pathway_sel_plots_manual$sig_pathways[[4]]$taxon) # four pathways 

# check the lipid related pathways 
lipid_pathways_plot <- rownames(lipid_path)

# check the amino acid related pathways 
amino_acid_pathways_plot <- aminoacid_path$taxon

# change the font color for the plot 
## first change the ones for the abx high fat mutant 

normalized_css_pathway_sel_plots_manual_modified <- normalized_css_pathway_sel_plots_manual
library(pheatmap)
library(grid)

labels_for_abx_hf_mutant <- normalized_css_pathway_sel_plots_manual_modified$normalized_table_sel[[3]]$OTUID[normalized_css_pathway_sel_plots_manual_modified$plots[[3]]$tree_row$order]
color_for_abx_hf_mutant <- ifelse(labels_for_abx_hf_mutant %in% overlap_pathways, "#f46d43", "black")
color_for_abx_hf_mutant <- ifelse(labels_for_abx_hf_mutant %in% lipid_pathways_plot, "#fdae61", color_for_abx_hf_mutant) # 
color_for_abx_hf_mutant <- ifelse(labels_for_abx_hf_mutant %in% amino_acid_pathways_plot, "#4575b4", color_for_abx_hf_mutant) 
normalized_css_pathway_sel_plots_manual_modified$plots[[3]]$gtable$grobs[[4]]$gp <- gpar(col=color_for_abx_hf_mutant)
normalized_css_pathway_sel_plots_manual_modified$plots[[3]]




## second change the ones for the abx high fat wt

labels_for_abx_hf_wt <- normalized_css_pathway_sel_plots_manual_modified$normalized_table_sel[[4]]$OTUID[normalized_css_pathway_sel_plots_manual_modified$plots[[4]]$tree_row$order]
color_for_abx_hf_wt <- ifelse(labels_for_abx_hf_wt %in% overlap_pathways, "#f46d43", "black")
color_for_abx_hf_wt <- ifelse(labels_for_abx_hf_wt %in% lipid_pathways_plot, "#fdae61", color_for_abx_hf_wt) # 
color_for_abx_hf_wt <- ifelse(labels_for_abx_hf_wt %in% amino_acid_pathways_plot, "#4575b4", color_for_abx_hf_wt) 
normalized_css_pathway_sel_plots_manual_modified$plots[[4]]$gtable$grobs[[4]]$gp <- gpar(col=color_for_abx_hf_wt)
normalized_css_pathway_sel_plots_manual_modified$plots[[4]]


normalized_css_pathway_sel_plots_manual_modified$path <- str_replace(normalized_css_pathway_sel_plots_manual_modified$path, "metagenomeSeq_heatmap_change", "metagenomeSeq_heatmap_change_modified_plot/")


batch_savepdf_v2 <- function(input = NA, output_dir = NA, filename = NA, otu = NA){
  
  if (nrow(otu) > 150){
    pdf(str_c(output_dir, filename, ".pdf"), width = 20, height = 40)
    print(input)
    dev.off()
  }else{
    pdf(str_c(output_dir, filename, ".pdf"), width = 20*0.8, height = 40*0.8)
    print(input)
    dev.off()
  }
  
  
}

normalized_css_pathway_sel_plots_manual_modified %>% 
  rowwise() %>% 
  mutate(list(batch_savepdf_v2(input = plots,
                            output_dir = path,
                            filename = compar,
                            otu = normalized_table_sel)))



save.image(str_c(output_dir, "/pathway_vis_manual_remove_fp_pheatmpa_change.RData"))

########################################################
# change color of the venn plot 
# 2024.04.24 
########################################################


output_dir <- "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/functional_analysis/output/plots/"
load(file = str_c(output_dir, "/pathway_vis_manual_remove_fp_pheatmpa_change.RData"))



venn_plot_function <- function(a = NA, b = NA, file = NA, input_file = NA){
  
  a_list <- input_file[[a]]$taxon
  b_list <- input_file[[b]]$taxon
  
  list_venn <- list('WT' = a_list, 'Mutant' = b_list)
  p <- ggvenn::ggvenn(list_venn, c("WT", "Mutant"), fill_color = c("#0000008c", "#ff00668c"), fill_alpha = 0.4) + ggtitle(file %>% str_replace_all("_", " "))
  
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


save_venn_plot <- function(venn_list = NA, qtag = NA, output_dir=NA){
  
  p_combine <- list()
  for (i in names(venn_list)){
    
    venn_list[[i]]$tag <- i 
    p_combine[[i]] <- ggarrange(plotlist = venn_list[[i]]$plots)
    ggsave(paste0(output_dir, i, "_", qtag, ".pdf"), p_combine[[i]], width = 9)
    
  }
  
  
}

metagenomseq_result_list_time_q_0.2_venn_manual_v2 <- lapply(sig_pathways_q_0.2_manual, compare_the_overlap_metagenomeseq)
save_venn_plot(metagenomseq_result_list_time_q_0.2_venn_manual_v2, qtag = "0.2", 
               output_dir=str_c(output_dir, "/metagenomeSeq_venn_manual/v2/"))

save.image(file = str_c(str_c(output_dir, "/metagenomeSeq_venn_manual/v2/"), "/venn_plot_v2.RData"))



