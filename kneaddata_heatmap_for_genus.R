# this script is to draw the bubble plot for HF diet in genera level 
# author: lu.zhang 
# date: 15,Jan,2023 


load("/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output_strict/sig_genera.timepoint/sig_species.RData")

###########################################################
## input 
###########################################################

metagenomseq_result_list_time[[1]][["ABX_High_fat_Mutant"]]
metagenomseq_result_list_time[[1]][["ABX_High_fat_WT"]]
metagenomseq_result_list_time_no_threshold[[1]][["ABX_High_fat_Mutant"]]
metagenomseq_result_list_time_no_threshold[[1]][["ABX_High_fat_WT"]]

###########################################################
## union of the sig.sp for both wt and mutant in ABX 
###########################################################

library(tidyverse)
all_sig_sp <- union(metagenomseq_result_list_time[[1]][["ABX_High_fat_Mutant"]]$taxon, metagenomseq_result_list_time[[1]][["ABX_High_fat_WT"]]$taxon) #35 genera 


mutant_abx_hf_sig <- metagenomseq_result_list_time_no_threshold[[1]][["ABX_High_fat_Mutant"]] %>% 
  filter(taxon %in% all_sig_sp)
wt_abx_hf_sig <- metagenomseq_result_list_time_no_threshold[[1]][["ABX_High_fat_WT"]] %>% 
  filter(taxon %in% all_sig_sp)


abx_hf_sig <- bind_rows(mutant_abx_hf_sig %>% select(time14, pvalues, adjPvalues, taxon) %>% mutate(group = "Day14 vs Day0 (Mutant)"),
                        wt_abx_hf_sig %>% select(time14, pvalues, adjPvalues, taxon) %>% mutate(group = "Day14 vs Day0 (WT)")) %>% 
  dplyr::rename(effect_size = time14)

###########################################################
## extract sig.sp for ABX, highfat diet 
###########################################################

bubble_input <- abx_hf_sig %>% 
  mutate(EffectSize = abs(effect_size), direction = (ifelse(effect_size >0, "increase", "decrease"))) %>% 
  mutate(species = str_split(taxon, pattern = "g__", simplify = T)[,2] %>% str_replace_all("_", " ")) %>% arrange(species) %>% 
  mutate(color = "gray")

bubble_input$direction[bubble_input$pvalues > 0.05] <- "No difference"
bubble_input$color[bubble_input$direction == "increase"] <- "#ef3647" # "red"
bubble_input$color[bubble_input$direction == "decrease"] <- "#298dc4" #"dodgerblue"

# annotation bar1 for mutant 
anno_bar1 <- bubble_input %>% select(species) %>% unique() %>% arrange(species) %>% 
  left_join(., bubble_input %>% filter(group == "Day14 vs Day0 (Mutant)") %>% select(color, species))
anno_bar1$color[is.na(anno_bar1$color)] <- "black"
# annotation bar2 for wt 
anno_bar2 <- bubble_input %>% select(species) %>% unique() %>% arrange(species) %>% 
  left_join(., bubble_input %>% filter(group == "Day14 vs Day0 (WT)") %>% select(color, species))
anno_bar2$color[is.na(anno_bar2$color)] <- "black"



bubble_input$adjP_range <- cut(
  bubble_input$adjPvalues, 
  right = FALSE,
  breaks = c(0, 0.01, 0.05, 0.1, 0.2, 1),
  labels = c("[0-0.01)", "[0.01, 0.05)", "[0.05, 0.1)", "[0.1, 0.2)", "[0.2, 1]")
)


text_color <- "black"
###########################################################
# main figure 
###########################################################

library(RColorBrewer)
library(grid)
p3 <- bubble_input %>% mutate(`Effect Size` = EffectSize, `FDR-correct P` = adjP_range) %>%
  ggplot(aes(x=group, y = species, color = `FDR-correct P`, size = `Effect Size`, fill = direction)) +
  geom_point()+
  coord_flip() + 
  theme_minimal() + 
  scale_color_manual(
    values = c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8","#4575B4"), #brewer.pal(6, "RdYlBu"), 
    labels = c("[0-0.01)", "[0.01, 0.05)", "[0.05, 0.1)", "[0.1, 0.2)", "[0.2, 1.0)") 
  ) + 
  scale_y_discrete(guide = guide_axis(angle = 45)) + 
  ylab('Significantly differentially abundant genera') +
  xlab('Group') +
  annotate(
    xmin = 2.75, xmax = 3.00,
    ymin = c(-Inf, seq(1,34) + 0.5), #c(-Inf, 3, 6),
    ymax = c(seq(1,34) + 0.5, Inf), #c(3, 6, Inf),
    geom = "rect",
    fill = anno_bar1$color, #c("dodgerblue", "gray", "red"),
    color = "black"
  ) + 
  annotate(
    xmin = 3.00, xmax = 3.25,
    ymin = c(-Inf, seq(1,34) + 0.5), #c(-Inf, 3, 6),
    ymax = c(seq(1,34) + 0.5, Inf), #c(3, 6, Inf),
    geom = "rect",
    fill = anno_bar2$color,
    color = "black"
  ) + 
  scale_fill_manual(
    values = c("increase" = "#ef3647", "decrease" = "#298dc4", "No difference" = "gray"),
    labels = c("Increase", "Decrease", "No difference"),
    name = "Direction of sig.DA genera"
  ) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1.5), axis.text.x=element_text(size=9, face = "italic", color = text_color), plot.margin = unit(c(2,1,1,1), "cm"), legend.key.size = unit(0.5, 'cm'), legend.key = element_rect(color = "black")) +  # up,right, bottom,  
  ggtitle("Bubble plot of significantly differentially abundant genera")


ggsave("/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output_strict/sig_genera.timepoint/heatmap/bubble_abx_hf_range_genera.pdf", p3, width = 20, height = 5)


# add prevalence 

phyloseq_objects_ten_prevalence_group # I can just select the abx high fat groups 

antibiotic_hf_mutant <- phyloseq_objects_ten_prevalence_group[[1]]$ABX_High_fat_Mutant
antibiotic_hf_wt <- phyloseq_objects_ten_prevalence_group[[1]]$ABX_High_fat_WT

# calculate the prevalence for each species 
species_list <- unique(bubble_input$taxon)

antibiotic_hf_mutant_df <- phyloseq::psmelt(antibiotic_hf_mutant) %>% 
  mutate(existence = case_when(Abundance > 0 ~ 1, Abundance == 0 ~ 0)) %>% 
  filter(OTU %in% species_list) %>% 
  mutate(num = 1) %>% 
  group_by(OTU) %>% 
  summarise(count_mutant = sum(existence), total_sample_mutant = sum(num)) %>% 
  ungroup() %>% 
  mutate(prevalence_mutant = count_mutant/total_sample_mutant) %>% 
  left_join(., mutant_abx_hf_sig, by = c("OTU" = "taxon")) %>% 
  rename(taxon = OTU, effect_size = time14) %>% 
  select(-`(Intercept)`, -scalingFactor) %>% 
  mutate(group = "mutant_hf_abx") %>% 
  rename(count = count_mutant, total_sample = total_sample_mutant, prevalence = prevalence_mutant)

antibiotic_hf_wt_df <- phyloseq::psmelt(antibiotic_hf_wt) %>% 
  mutate(existence = case_when(Abundance > 0 ~ 1, Abundance == 0 ~ 0)) %>% 
  filter(OTU %in% species_list) %>% 
  mutate(num = 1) %>% 
  group_by(OTU) %>% 
  summarise(count_wt = sum(existence), total_sample_wt = sum(num)) %>% 
  ungroup() %>% 
  mutate(prevalence_wt = count_wt/total_sample_wt) %>% 
  left_join(., wt_abx_hf_sig, by = c("OTU" = "taxon")) %>% 
  rename(taxon = OTU, effect_size = time14) %>% 
  select(-`(Intercept)`, -scalingFactor) %>% 
  mutate(group = "wt_hf_abx") %>% 
  rename(count = count_wt, total_sample = total_sample_wt, prevalence = prevalence_wt)


# average relative abundance for each species : added data 
input_phyloseq_genus <- tax_glom(physeq_object_bac_rm_cf, "Genus")

phyloseq_input.rel <- microbiome::transform(input_phyloseq_genus, "compositional") 
phyloseq_input.rel.df <- phyloseq::psmelt(phyloseq_input.rel) %>%
  mutate(OTU_new = OTU %>% str_extract(., pattern = "^k__.*g__.*\\|s__") %>% str_remove(., "\\|s__")) %>% 
  filter(OTU_new %in% species_list) %>% 
  filter(antibiotic == "ABX", diet == "High_fat") %>% 
  select(OTU_new, Sample, Abundance, diet, antibiotic, group)%>% 
  group_by(OTU_new, group) %>% 
  summarise(mean_abundance = mean(Abundance)) %>% 
  ungroup() %>% 
  mutate(group = case_when(group %in% "Mutant" ~ "mutant_hf_abx", group %in% "WT" ~ "wt_hf_abx"))



sig_species_prevalence_abundance <- bind_rows(antibiotic_hf_mutant_df, antibiotic_hf_wt_df) %>% 
  left_join(., phyloseq_input.rel.df, by = c("taxon" = "OTU_new", "group" = "group"))

write_tsv(sig_species_prevalence_abundance, "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output_strict//sig_genera.timepoint/heatmap/bubble_plot_prevalence.tsv")
saveRDS(sig_species_prevalence_abundance, "/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output_strict//sig_genera.timepoint/heatmap/bubble_plot_prevalence.rds")
save.image("/data/202210_mouse_project/Result/kneaddata_second_time_kraken2/output_strict//sig_genera.timepoint/heatmap/bubble_plot_heatmap.RData")

