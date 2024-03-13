## this script is for all analyses for all datasets. write functions


############### loading packages
############################################
require(tidyverse)
require(Maaslin2)
require(dplyr)
require(tidyr)
require(ggplot2)
require(DESeq2)
require(ggsignif)
require(nnet)
require(gtsummary)
require(tidyr)
require(broom)
require(ggvenn)
require(ggVennDiagram)
require(ggsci)
require(ggtext)
require(readxl)
require(flextable)
library(cooccur)
library(visNetwork)
require(rstatix)
require(ggpubr)
require(ggrepel)
require(tidygraph)
require(ggraph)
require(ggside)
require(ggpmisc)
require(irr)
require(psych)
require(enrichplot)
require(MicrobiomeProfiler)

# log transformation from Maaslin2
LOG <- function(x) {
  y <- replace(x, x==0, min(x[x>0]) /2)
  return(log2(y))
}

############################ frequency plots
############################################

chisq_sex <-
  lapply(c("crcahus_biopsi_ASV", "crcahus_faeces_ASV", "norccap_ASV", "norccap_MG", "CRCbiome_meta_MG", "curated_MG"), function(dataset) { 
    if (dataset == "crcahus_biopsi_ASV") tmp_dataset <- crcahus_biopsi_ASV
    if (dataset == "crcahus_faeces_ASV") tmp_dataset <- crcahus_faeces_ASV
    if (dataset == "norccap_ASV") tmp_dataset <- norccap_ASV
    if (dataset == "norccap_MG") tmp_dataset <- norccap_MG
    if (dataset == "CRCbiome_meta_MG") tmp_dataset <- CRCbiome_meta_MG
    if (dataset == "curated_MG") tmp_dataset <- curated_MG
    
    tmp_freq_table <- tmp_dataset %>% 
      select(c(sex, pres_sp, pres_succinatutens, pres_faecium)) %>% 
      tbl_summary(by="sex") %>% 
      add_p(
        pvalue_fun = function(x) style_pvalue(x, digits=3)) %>% 
      add_stat_label() %>% 
      modify_caption("Table 1. Chi2 test") %>% 
      as_gt() %>% 
      .$`_data` %>% 
      mutate(dataset = dataset) %>% 
      as.data.frame()
  }) %>% 
  bind_rows()


## Making 
(barplot_sex <- freq_table_sex %>% 
    select(sex, freq, taxa, dataset) %>% 
    mutate(taxa = recode(taxa, "pres_succinatutens" = "P. succinatutens",
                         "pres_faecium" = "P. faecium",
                         "pres_sp" = "P. sp")) %>% 
    mutate(dataset = recode(dataset,                          
                            "norccap_ASV" = "NORCCAP 16S",
                            "norccap_MG" = "NORCCAP MG",
                            "CRCbiome_meta_MG" = "CRCbiome", 
                            "curated_MG" = "CuratedMG",
                            "crcahus_biopsi_ASV" = "CRCAhus biopsy",
                            "crcahus_faeces_ASV" = "CRCAhus feces")) %>%     
    dplyr::rename(Dataset = dataset) %>% 
    mutate(taxa=factor(taxa, levels=c( "P. faecium", "P. sp","P. succinatutens"))) %>% 
    mutate(Dataset = factor(Dataset, levels=c("CRCAhus biopsy", "CRCAhus feces", "NORCCAP 16S", "NORCCAP MG", "CRCbiome", "CuratedMG"))) %>%
    ggplot(aes(x=taxa, y=freq, fill=sex))+
    geom_bar(stat="identity", width = 0.8, position = position_dodge(),colour="black") +
    scale_fill_manual(values=c("#00A087FF","#F39B7FFF", "#3C5488FF"))+  
    scale_x_discrete(breaks=c("P. succinatutens",
                              "P. faecium",
                              "P. sp"),
                     labels=c("P. succinatutens",
                              "P. faecium",
                              "P. sp")) +
    labs(x=NULL,
         y="Prevalence") +
    facet_wrap(~ Dataset, ncol=1) +
    coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_markdown(size=10),
          axis.text.y = element_markdown(size=10, face="italic"),
          axis.title.x = element_markdown(size=10),
          legend.title = element_markdown(size=12),
          legend.title.align = 0.5,
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size=10),
          legend.key.height = unit(12, "pt")))

############################################

#################### overlap between species
############################################

freq_upset <-
  lapply(c("crcahus_biopsi_ASV", "crcahus_faeces_ASV", "norccap_ASV", "norccap_MG", "CRCbiome_meta_MG", "CRCbiome_mags_MG", "curated_MG"), function(dataset) { 
    if (dataset == "crcahus_biopsi_ASV") tmp_dataset <- crcahus_biopsi_ASV %>% select(sample_id, pres_succinatutens, pres_faecium, pres_sp) %>% mutate(dataset = "crcahus_biopsi_ASV")
    if (dataset == "crcahus_faeces_ASV") tmp_dataset <- crcahus_faeces_ASV %>% select(sample_id, pres_succinatutens, pres_faecium, pres_sp) %>% mutate(dataset = "crcahus_faeces_ASV")
    if (dataset == "norccap_ASV") tmp_dataset <- norccap_ASV %>% select(sample_id, pres_succinatutens, pres_faecium, pres_sp) %>% mutate(dataset = "norccap_ASV")
    if (dataset == "norccap_MG") tmp_dataset <- norccap_MG %>% select(sample_id, pres_succinatutens, pres_faecium, pres_sp) %>% mutate(dataset = "norccap_MG")
    if (dataset == "CRCbiome_meta_MG") tmp_dataset <- CRCbiome_meta_MG %>% select(sample_id, pres_succinatutens, pres_faecium, pres_sp) %>% mutate(dataset = "CRCbiome_meta_MG")
    if (dataset == "curated_MG") tmp_dataset <- curated_MG %>% select(sample_id, pres_succinatutens, pres_faecium, pres_sp) %>% mutate(dataset = "curated_MG")
    
    lapply(c("pres_succinatutens", "pres_faecium", "pres_sp"), function(taxa) {
      
      tmp_freq_table <- tmp_dataset %>% 
        dplyr::rename(tmp = all_of(taxa)) %>% 
        select(tmp) %>% 
        group_by(tmp) %>% 
        summarise(PA = n()) %>% 
        ungroup() %>% 
        pivot_wider(names_from=tmp, values_from=PA) %>% 
        rowwise() %>% 
        mutate(freq = yes/(yes+no)) %>% 
        mutate(dataset = paste0(dataset)) %>% 
        mutate(taxa = taxa) %>% 
        as.data.frame()
    }) %>% 
      bind_rows()
  }) %>% 
  bind_rows()


(upset_plot <- frequency_upset %>% 
    select(freq, taxa, dataset) %>% 
    filter(taxa != "succ_faec_sp") %>% 
    mutate(freq = freq*100) %>%  
    mutate(taxa = factor(taxa, levels=c("pres_faecium", "pres_succinatutens", "pres_sp",  "faec_sp", "succ_faec", "succ_sp"))) %>% 
    mutate(dataset = recode(dataset, "CRCAhus_biopsi" = "CRCAhus biopsy",
                            "CRCAhus_faeces" = "CRCAhus feces", 
                            "norccap_ASV" = "NORCCAP 16S", 
                            "norccap_MG" = "NORCCAP MG",
                            "CRCbiome_meta_MG" = "CRCbiome", 
                            "curated_MG" = "CuratedMG")) %>% 
    dplyr::rename(Dataset = dataset) %>% 
    mutate(Dataset = factor(Dataset, levels=c("CRCAhus biopsy", "CRCAhus feces", "NORCCAP 16S", "NORCCAP MG", "CRCbiome", "CuratedMG"))) %>%
    ggplot(aes(x=taxa, y=freq, fill=Dataset))+
    geom_bar(stat="identity") +
    scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF", "#3C5488FF", "#B09C85FF", "#F39B7FFF", "#8491B4FF"))+
    labs(x=NULL,
         y=NULL) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_markdown(size=14),
          legend.title = element_markdown(size=18),
          legend.title.align = 0.5,
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.text = element_text(size=18),
          legend.key.height = unit(18, "pt")))

############################################

##################### linear model + heatmap
############################################
## performing a linear model with log transformed abundance data as dependent variable and clinical group as independent variable adjusted for sex, age and geography (screening center or study)
lm_adjusted <- lapply(c("CRCbiome_meta_MG", "curated_MG"), function(dataset) { 
  if (dataset == "CRCbiome_meta_MG") tmp_dataset <- CRCbiome_meta_MG
  if (dataset == "curated_MG") tmp_dataset <- curated_MG
  
  lapply(c("log_succinatutens", "log_sp", "log_faecium"), function(taxa) {
    
    tmp <- tmp_dataset %>% 
      dplyr::rename(tmp = all_of(taxa)) %>% 
      lm(tmp ~ cc_status + sex + age + geography, data = ., na.action = na.omit) %>% 
      tidy() %>% 
      mutate(dataset = dataset, taxa = taxa) %>% 
      filter(term != "(Intercept)")
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows()

heatmap_lm <- lm_results %>% 
  filter(dataset == "CRCbiome_meta_MG" | dataset == "curated_MG") %>% 
  mutate(term = recode(term, "sexMale" = "Male vs Female",
                       "cc_statusCRC" = "CRC vs Controls",
                       "cc_statusAdenoma" = "Adenoma vs Controls")) %>% 
  mutate(term = recode(term, "sexmale" = "Male vs Female",
                       "cc_statusadenoma" = "Adenoma vs Controls")) %>% 
  mutate(taxa = recode(taxa, "log_succinatutens" = "P. succinatutens",
                       "log_faecium" = "P. faecium",
                       "log_sp" = "P. sp")) %>% 
  mutate(pstar = case_when(
    p.value < 0.05 ~ "*",
    TRUE ~ "")) %>% 
  mutate(dataset = recode(dataset, 
                          "CRCbiome_meta_MG" = "CRCbiome",
                          "curated_MG" = "CuratedMG")) %>% 
  mutate(dataset=factor(dataset, levels=c("CRCbiome", "CuratedMG"))) %>% 
  mutate(taxa=factor(taxa, levels=c("P. succinatutens", "P. sp", "P. faecium"))) %>% 
  ggplot(aes(x=taxa, y=term, fill=estimate)) +
  geom_tile() +
  geom_text(aes(label = pstar), size=1.5/0.35) + 
  scale_x_discrete(breaks=c("P. succinatutens",
                            "P. faecium",
                            "P. sp"),
                   labels=c("P.<br>succinatutens",
                            "P.<br>faecium",
                            "P.<br>sp"),
                   expand = c(0,0),
                   position = "top") +
  scale_y_discrete(breaks=c("Male vs Female",
                            "CRC vs Controls",
                            "Adenoma vs Controls"),
                   labels=c("Male vs<br>Female",
                            "CRC vs<br>Controls",
                            "Adenoma vs<br>Controls")) +
  scale_fill_gradient2(name = "Log2FC",
                       low = "#3C5488FF", mid = "#FFFFFF", high = "#F39B7FFF",
                       expand = c(0,0)) +
  labs(x=NULL,
       y=NULL) +
  theme_bw() +
  theme(axis.text.x.top = element_markdown(vjust=0.5, size=10, face="italic"),
        axis.text.y = element_markdown(size = 10),
        legend.title = element_markdown(size=10),
        legend.title.align = 0.3,
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size=10),
        legend.key.height = unit(10, "pt")
  ) +
  coord_fixed(ratio = 0.9) +
  facet_wrap("dataset", ncol = 2, strip.position = "bottom") +
  theme(strip.text = element_text(size = 10))

############################################

####################### correlation analyses
############################################
list_of_dfs <- list(crcahus_biopsi_ASV = crcahus_biopsi_ASV, crcahus_faeces_ASV = crcahus_faeces_ASV, norccap_maaslin_ASV = norccap_maaslin_ASV, norccap_maaslin_MG = norccap_maaslin_MG, CRCbiome_meta_maaslin_MG = CRCbiome_meta_maaslin_MG, curated_maaslin_MG = curated_maaslin_MG)
lapply(list_of_dfs, dim)
names(list_of_dfs)

filter_bacteria <- function(df) {
  # Calculate the proportion of non-zero entries for each column
  proportion_nonzero <- df %>%
    select(-1) %>%  # Exclude the first column
    summarise(across(everything(), ~ mean(.x != 0)))
  
  # Identify columns to keep (including the first column)
  columns_to_remove <- names(proportion_nonzero)[proportion_nonzero < 0.05]
  
  # Filter the dataframe to keep only the identified columns
  df %>% select(-all_of(columns_to_remove))
}

filtered_dfs <- map(list_of_dfs, filter_bacteria)
lapply(filtered_dfs, dim)
names(filtered_dfs)

## performing correlation
correlation_filtered <- lapply(names(filtered_dfs), function(dataset_name) {
  tmp_dataset <- filtered_dfs[[dataset_name]]
  
  tmp_dataset <- tmp_dataset %>% 
    column_to_rownames("sample_id") %>% 
    LOG() 
  
  vars <- c("Phascolarctobacterium_faecium","Phascolarctobacterium_sp_CAG_266","Phascolarctobacterium_succinatutens")
  vars2 <- names(tmp_dataset)
  
  tmp <- cor_test(tmp_dataset, method = "spearman", conf.level = 0.95, use = "na.or.complete", vars=all_of(vars), vars2=all_of(vars2)) %>% 
    mutate(dataset=dataset_name)
}) %>% 
  bind_rows()

## making correlation table for import into cytoscape
# aslo want to set phylum level 

taxonomy <- read_tsv("metaphlan_table.tsv", skip = 1, col_names=T) %>% 
  select(clade_name) %>%
  mutate(clade_name = gsub("[|]", "_", clade_name)) %>% 
  mutate(clade_name = gsub("k__", "", clade_name)) %>% 
  filter(grepl("s__", clade_name)) %>% 
  separate(clade_name, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "__", 
           remove = FALSE, 
           extra = "merge", 
           fill = "right") %>% 
  mutate(across(Kingdom:Genus, ~sub("_.$", "", .))) %>% 
  select(-c(clade_name, Kingdom))

pair_counts <- correlation_filtered %>% 
  select(-c(statistic, method, p)) %>% 
  filter(p_adjust < 0.05) %>% 
  group_by(var1, var2) %>%
  summarise(count = n_distinct(dataset), .groups = "keep")

## Making dataframe as input for cytoscape and correaltion network plot
correlation_cytoscape <- correlation_filtered %>% 
  select(-c(statistic, method, p)) %>% 
  filter(p_adjust < 0.05) %>% 
  mutate(genus = str_extract(var2, "^[^_]*")) %>% 
  mutate(direction = ifelse(cor > 0, "positive", "negative")) %>% 
  left_join(pair_counts, by = c("var1", "var2")) %>% 
  mutate(group = ifelse(count == 3, "all", NA),
         group = ifelse(count == 2, "curated_crcbiome", group),
         group = ifelse(count == 2 & (dataset == "norccap_maaslin_MG" | dataset == "CRCbiome_meta_maaslin_MG" ) & var1 == "Phascolarctobacterium_succinatutens" & var2 == "Olsenella_scatoligenes",
                        "norccap_crcbiome", group),
         group = ifelse(count == 2 & (dataset == "norccap_maaslin_MG" | dataset == "curated_maaslin_MG" ) & var1 == "Phascolarctobacterium_sp_CAG_266" & var2 == "Alistipes_inops",
                        "norccap_curated", group),
         group = ifelse(count == 1 & dataset == "curated_maaslin_MG", "curated", group),
         group = ifelse(count == 1 & dataset == "CRCbiome_meta_maaslin_MG", "crcbiome", group),
         group = ifelse(count == 1 & dataset == "norccap_maaslin_MG", "norccap", group)) %>% 
  mutate(phasco = var1) %>% 
  mutate(corr = ifelse(count > 1, var2, "")) %>% 
  group_by(var1, var2) %>%
  mutate(var3 = ifelse(direction == lead(direction) | direction == lag(direction), "same", "opposite")) %>%
  ungroup() %>% 
  filter(var3 == "same" & dataset == "curated_maaslin_MG") %>% 
  left_join(taxonomy, by=c("var2"="Species"))

############################################

###################### beta diversity groups
############################################

## filtering on those samples that contain reads to phasco
CRCbiome_phasco <- c(names_faec, names_succi, names_sp) %>% as.data.frame() %>% 
  dplyr::rename("sample_id" = ".") %>% 
  mutate(sample_id = str_replace(sample_id, "S-", "S_"),
         sample_id = str_extract(sample_id, "S_\\d+")) %>% 
  left_join(CRCbiome_meta_maaslin_MG) %>% 
  column_to_rownames("sample_id")

CRCbiome_phasco_meta <- c(names_faec, names_succi, names_sp) %>% as.data.frame() %>% 
  dplyr::rename("sample_id" = ".") %>% 
  mutate(species = ifelse(sample_id %in% names_succi, "succi", NA)) %>% 
  mutate(species = ifelse(sample_id %in% names_sp, "sp", species)) %>% 
  mutate(species = ifelse(sample_id %in% names_faec, "faec", species)) %>% 
  mutate(sample_id = str_replace(sample_id, "S-", "S_"),
         sample_id = str_extract(sample_id, "S_\\d+")) %>% 
  left_join(CRCbiome_metadata) 

## beta diversity
source("utils.R")
### only including species ones

CRCbiome_metadata_adonis <- CRCbiome_metadata_adonis %>% filter(phasco_group != "none" & phasco_group != "2grp")
CRCbiome_meta_adonis <- CRCbiome_meta_adonis %>% rownames_to_column("sample_id") %>% inner_join(CRCbiome_metadata_adonis %>% rownames_to_column("sample_id") %>% select(sample_id)) %>% 
  column_to_rownames("sample_id")

set.seed(1234)
adonis_in <-  CRCbiome_meta_adonis %>% 
  as.matrix() %>% 
  vegan::vegdist(method= "bray") 

set.seed(1234)
adonis_out <- vegan::adonis2(adonis_in ~ phasco_group + sex + age + geography, by = "margin" , data=CRCbiome_metadata_adonis)

set.seed(1234)
bray_dist <- CRCbiome_meta_adonis %>%
  as.matrix() %>% 
  vegan::vegdist(method= "bray") 

p <- 
  bray_dist %>% 
  dist_to_PCoA(group_var = CRCbiome_meta_adonis %>% 
                 rownames_to_column("sample_id") %>% 
                 dplyr::select(sample_id) %>% 
                 left_join(CRCbiome_metadata_adonis %>% rownames_to_column("sample_id") %>% dplyr::select(sample_id, phasco_group)) %>% 
                 pull(phasco_group)) %>% 
  plot_pcoa(dim_1 = "PCoA1", dim_2 = "PCoA2") +
  ## 
  annotate("text",
           label = paste("Permanova p =", round(adonis_out$`Pr(>F)`[1], 3), "R2", round(adonis_out$R2[1], 3)),
           x = 0,
           y=0.4,
           color = "grey30",
           size = 4) 
p

############################################

################################# qPCR plots
############################################

(p <- qpcr_all %>% 
    select(sample_id, pcr_succinatutens, Phascolarctobacterium_succinatutens, dataset) %>% 
    mutate(pcr_succinatutens = as.numeric(pcr_succinatutens), Phascolarctobacterium_succinatutens = as.numeric(Phascolarctobacterium_succinatutens)) %>%
    filter(pcr_succinatutens > 0 | Phascolarctobacterium_succinatutens > 0) %>% 
    mutate(log_succinatutens = LOG(Phascolarctobacterium_succinatutens), log_qpcr = LOG(pcr_succinatutens)) %>% 
    mutate(pcr_succinatutens = ifelse(is.na(pcr_succinatutens), 0, pcr_succinatutens)) %>% 
    ggplot(aes(x=pcr_succinatutens, y=Phascolarctobacterium_succinatutens)) +
    stat_poly_line()+
    geom_point() + 
    labs(title = expression(italic("P. succinatutens")),
         x = "Relative abundance qPCR", y = "Relative abundance NGS") + 
    facet_wrap("dataset", scales = "free", ncol=1, strip.position = "left") +
    theme_bw() + 
    theme(strip.text = element_text(size = 12),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text = element_text(size = 11))) 


(p1 <- qpcr_all %>% 
    select(sample_id, pcr_sp, Phascolarctobacterium_sp_CAG_266, dataset) %>% 
    mutate(pcr_sp = as.numeric(pcr_sp), Phascolarctobacterium_sp_CAG_266 = as.numeric(Phascolarctobacterium_sp_CAG_266)) %>%
    filter(pcr_sp > 0 | Phascolarctobacterium_sp_CAG_266 > 0) %>% 
    mutate(log_sp = LOG(Phascolarctobacterium_sp_CAG_266), log_qpcr = LOG(pcr_sp)) %>% 
    mutate(pcr_sp = ifelse(is.na(pcr_sp), 0, pcr_sp)) %>% 
    ggplot(aes(x=pcr_sp, y=Phascolarctobacterium_sp_CAG_266)) +
    stat_poly_line()+
    geom_point() + 
    labs(title = expression(italic("P. sp")),
         x = "Relative abundance qPCR", y = "Relative abundance NGS") + 
    facet_wrap("dataset", scales = "free", ncol=1) +
    theme_bw()+
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text = element_text(size = 11)))


(p2 <- qpcr_all %>% 
    select(sample_id, pcr_faecium, Phascolarctobacterium_faecium, dataset) %>% 
    mutate(pcr_faecium = as.numeric(pcr_faecium), Phascolarctobacterium_faecium = as.numeric(Phascolarctobacterium_faecium)) %>%
    filter(pcr_faecium > 0 | Phascolarctobacterium_faecium > 0) %>% 
    mutate(log_faecium = LOG(Phascolarctobacterium_faecium), log_qpcr = LOG(pcr_faecium)) %>% 
    mutate(pcr_faecium = ifelse(is.na(pcr_faecium), 0, pcr_faecium)) %>% 
    ggplot(aes(x=pcr_faecium, y=Phascolarctobacterium_faecium)) +
    stat_poly_line()+
    geom_point() +
    labs(title = expression(italic("P. faecium")),
         x = "Relative abundance qPCR", y = "Relative abundance NGS") + 
    facet_wrap("dataset", scales = "free", ncol=1) +
    theme_bw() + 
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text = element_text(size = 11)))   

ggarrange(p, p1, p2, ncol=3)

correlation_succ <- lapply(c("norccap_qpcr", "norccap_qpcr_MG", "crcbiome_qpcr", "crcahus_qpcr"), function(dataset){
  if (dataset == "norccap_qpcr") tmp_dataset <- norccap_qpcr
  if (dataset == "norccap_qpcr_MG") tmp_dataset <- norccap_qpcr_MG
  if (dataset == "crcahus_qpcr") tmp_dataset <- crcahus_qpcr
  if (dataset == "crcbiome_qpcr") tmp_dataset <- crcbiome_qpcr
  
  
  tmp <- tmp_dataset %>% 
    cor.test(~ pcr_succinatutens + Phascolarctobacterium_succinatutens , data=., method="spearman") %>% 
    tidy() %>% 
    tibble() %>% 
    mutate(dataset=dataset, taxa="succ")
}) %>% 
  bind_rows()

correlation_faecium <- lapply(c("norccap_qpcr", "norccap_qpcr_MG", "crcbiome_qpcr", "crcahus_qpcr"), function(dataset){
  if (dataset == "norccap_qpcr") tmp_dataset <- norccap_qpcr
  if (dataset == "norccap_qpcr_MG") tmp_dataset <- norccap_qpcr_MG
  if (dataset == "crcahus_qpcr") tmp_dataset <- crcahus_qpcr
  if (dataset == "crcbiome_qpcr") tmp_dataset <- crcbiome_qpcr
  
  
  tmp <- tmp_dataset %>% 
    cor.test(~ pcr_faecium + Phascolarctobacterium_faecium, data=., method="spearman") %>% 
    tidy() %>% 
    tibble() %>% 
    mutate(dataset=dataset, taxa="faec")
}) %>% 
  bind_rows()

correlation_sp <- lapply(c("norccap_qpcr", "norccap_qpcr_MG", "crcbiome_qpcr", "crcahus_qpcr"), function(dataset){
  if (dataset == "norccap_qpcr") tmp_dataset <- norccap_qpcr
  if (dataset == "norccap_qpcr_MG") tmp_dataset <- norccap_qpcr_MG
  if (dataset == "crcahus_qpcr") tmp_dataset <- crcahus_qpcr
  if (dataset == "crcbiome_qpcr") tmp_dataset <- crcbiome_qpcr
  
  
  tmp <- tmp_dataset %>% 
    cor.test(~ pcr_sp + Phascolarctobacterium_sp_CAG_266, data=., method="spearman") %>% 
    tidy() %>% 
    tibble() %>% 
    mutate(dataset=dataset, taxa="sp")
}) %>% 
  bind_rows()


# cohens kappa
kappa_succi <- lapply(c("norccap_qpcr", "norccap_qpcr_MG", "crcbiome_qpcr", "crcahus_qpcr"), function(dataset){
  if (dataset == "norccap_qpcr") tmp_dataset <- norccap_qpcr
  if (dataset == "norccap_qpcr_MG") tmp_dataset <- norccap_qpcr_MG
  if (dataset == "crcahus_qpcr") tmp_dataset <- crcahus_qpcr
  if (dataset == "crcbiome_qpcr") tmp_dataset <- crcbiome_qpcr
  
  x <- tmp_dataset %>% mutate(pres_succ_ngs = ifelse(Phascolarctobacterium_succinatutens > 0, 1, 0)) %>% 
    mutate(pres_succ_qpcr = ifelse(pcr_succinatutens > 0, 1, 0)) %>% 
    select(pres_succ_ngs, pres_succ_qpcr) %>% 
    as.data.frame() %>% 
    cohen.kappa() %>% 
    tidy() %>% 
    mutate(dataset = dataset, species = "P.succi")
}) %>% 
  bind_rows()

kappa_faec <- lapply(c("norccap_qpcr", "norccap_qpcr_MG", "crcbiome_qpcr", "crcahus_qpcr"), function(dataset){
  if (dataset == "norccap_qpcr") tmp_dataset <- norccap_qpcr
  if (dataset == "norccap_qpcr_MG") tmp_dataset <- norccap_qpcr_MG
  if (dataset == "crcahus_qpcr") tmp_dataset <- crcahus_qpcr
  if (dataset == "crcbiome_qpcr") tmp_dataset <- crcbiome_qpcr
  
  x <- tmp_dataset %>% mutate(pres_faec_ngs = ifelse(Phascolarctobacterium_faecium > 0, 1, 0)) %>% 
    mutate(pres_faec_qpcr = ifelse(pcr_faecium > 0, 1, 0)) %>% 
    select(pres_faec_ngs, pres_faec_qpcr) %>% 
    as.data.frame() %>% 
    cohen.kappa() %>% 
    tidy() %>% 
    mutate(dataset = dataset, species = "P.faec")
}) %>% 
  bind_rows()

kappa_sp <- lapply(c("norccap_qpcr", "norccap_qpcr_MG", "crcbiome_qpcr", "crcahus_qpcr"), function(dataset){
  if (dataset == "norccap_qpcr") tmp_dataset <- norccap_qpcr
  if (dataset == "norccap_qpcr_MG") tmp_dataset <- norccap_qpcr_MG
  if (dataset == "crcahus_qpcr") tmp_dataset <- crcahus_qpcr
  if (dataset == "crcbiome_qpcr") tmp_dataset <- crcbiome_qpcr
  
  x <- tmp_dataset %>% mutate(pres_sp_ngs = ifelse(Phascolarctobacterium_sp_CAG_266 > 0, 1, 0)) %>% 
    mutate(pres_sp_qpcr = ifelse(pcr_sp > 0, 1, 0)) %>% 
    select(pres_sp_ngs, pres_sp_qpcr) %>% 
    as.data.frame() %>% 
    cohen.kappa() %>% 
    tidy() %>% 
    mutate(dataset = dataset, species = "P.sp")
}) %>% 
  bind_rows()

############################################

########## looking into lifestyle/demography
############################################
# lifestyle and demography variables 

variables_cont <- 
  c("Energi_kcal", 
    "Prot_energi", 
    "Karboh_energi", 
    "Sukker_energi", 
    "Fiber", 
    "Fett_energi", 
    "Mettet_energi", 
    "C_enum_energi",
    "C_flerum_energi", 
    "Trans_u_energi", 
    "Alko")

variables_cat <- 
  c("kjonn", 
    "geography",
    "BMI",
    "PhysAct_Score",
    "Smoking", 
    "Snus",
    "Utdanning", 
    "Sivilstatus_cat2",
    "Arbeid_lump",
    "Nasj_cat2",
    "wcrf_index_main", 
    "Antibiotics",
    "Antacids") 

lm_adjusted_diet_log <- lapply(c("log_succinatutens", "log_sp", "log_faecium"), function(taxa) {
  lapply(variables_log, function(diet) {
    
    tmp <- CRCbiome_metadata_all %>% 
      dplyr::rename(tmp_taxa = all_of(taxa)) %>% 
      dplyr::rename(tmp_diet = all_of(diet)) %>% 
      lm(tmp_taxa ~ LOG(tmp_diet) + sex + age, data = ., na.action = na.omit) %>% 
      tidy() %>% 
      mutate(diet = diet, taxa = taxa) %>% 
      filter(term != "(Intercept)", term != "sexMale", term != "age")
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows()
lm_adjusted_diet_log <- lm_adjusted_diet_log %>% mutate(padj = p.adjust(p.value))

lm_adjusted_diet_nonlog <- lapply(c("log_succinatutens", "log_sp", "log_faecium"), function(taxa) {
  lapply(variables_nonlog, function(diet) {
    
    tmp <- CRCbiome_metadata_all %>% 
      dplyr::rename(tmp_taxa = all_of(taxa)) %>% 
      dplyr::rename(tmp_diet = all_of(diet)) %>% 
      lm(tmp_taxa ~ tmp_diet + sex + age, data = ., na.action = na.omit) %>% 
      tidy() %>% 
      mutate(diet = diet, taxa = taxa) %>% 
      filter(term != "(Intercept)", term != "sexMale", term != "age")
  }) %>% 
    bind_rows()
}) %>% 
  bind_rows()

lm_adjusted_diet_nonlog <- lm_adjusted_diet_nonlog %>% mutate(padj = p.adjust(p.value))

## four groups of FIT values
CRCbiome_metadata_all <- CRCbiome_metadata_all %>% mutate(group_FIT = FIT_value/5) %>%  
  mutate(group_FIT = case_when(
    group_FIT <= 20 ~ "group1",
    group_FIT >= 20 & group_FIT < 35 ~ "group2",
    group_FIT >= 35 & group_FIT < 70 ~ "group3",
    group_FIT >= 70 ~ "group4")) 

require(MASS)
model <- polr(as.factor(group_FIT) ~ LOG(Phascolarctobacterium_sp_CAG_266) + sex + age + geography, data = CRCbiome_metadata_all, Hess = TRUE)
ctable <- coef(summary(model))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
cbind(ctable, "p value" = p)
confint(model)

############################################

##################################### roary
###########################################
roary_succi <- read.csv("Phascolarctobacterium_A_succinatutens/gene_presence_absence.csv")
roary_mat_succi <- read.table("Phascolarctobacterium_A_succinatutens/gene_presence_absence.Rtab", header=T)
roary_succi_long <- read.delim("Phascolarctobacterium_A_succinatutens/gene_presence_long.tsv")
colnames(roary_mat_succi) <- str_replace(colnames(roary_mat_succi), "^S\\.", "S-")
names_succi <- colnames(roary_mat_succi[2:53]) 

roary_faec <- read.csv("Phascolarctobacterium_faecium/gene_presence_absence.csv")
roary_mat_faec <- read.table("Phascolarctobacterium_faecium/gene_presence_absence.Rtab", header=T)
roary_faec_long <- read.delim("Phascolarctobacterium_faecium/gene_presence_long.tsv")
colnames(roary_mat_faec) <- str_replace(colnames(roary_mat_faec), "^S\\.", "S-")
names_faec <- colnames(roary_mat_faec[2:132])

roary_sp <- read.csv("Phascolarctobacterium_sp/gene_presence_absence.csv")
roary_mat_sp <- read.table("Phascolarctobacterium_sp/gene_presence_absence.Rtab", header=T)
roary_sp_long <- read.delim("Phascolarctobacterium_sp/gene_presence_long.tsv")
colnames(roary_mat_sp) <- str_replace(colnames(roary_mat_sp), "^S\\.", "S-")
names_sp <- colnames(roary_mat_sp[2:33])

roary_all <- read.csv("Phascolarcto/gene_presence_absence.csv")
roary_mat_all <- read.table("Phascolarcto/gene_presence_absence.Rtab", header=T)
roary_all_long <- read.delim("Phascolarcto/gene_presence_long.tsv")
colnames(roary_mat_all) <- str_replace(colnames(roary_mat_all), "^S\\.", "S-")
names_all <- colnames(roary_mat_all[2:222])

## counting average number of gene clusters per genome
sample_sums <- colSums(roary_mat_all[, -1]) 
mean(sample_sums)

## also for each species

sample_sums_succi <- roary_mat_all %>% select(all_of(names_succi)) 
sample_sums_succi <- colSums(sample_sums_succi[, -1]) 
mean(sample_sums_succi)

sample_sums_sp <- roary_mat_all %>% select(all_of(names_sp)) 
sample_sums_sp <- colSums(sample_sums_sp[, -1]) 
mean(sample_sums_sp)

sample_sums_faec <- roary_mat_all %>% select(all_of(names_faec))
sample_sums_faec <- colSums(sample_sums_faec[, -1]) 
mean(sample_sums_faec)

## counting fraction of core genes per genome
tmp <- roary_mat_all %>% rowwise() %>%  mutate(percentage = sum(c_across(-1)) / 5 * 100, 
                                               all_group = case_when(
                                                 percentage >= 95 ~ "genus_core",
                                                 percentage >= 15 ~ "shell",
                                                 TRUE ~ "cloud")) %>% 
  filter(all_group == "genus_core") %>% 
  select(-c(percentage, all_group))

sample_sums_core <- colSums(tmp[, -1]) 
mean(sample_sums_core)

## all
roary_mat_all <- roary_mat_all %>% 
  rowwise() %>%
  mutate(perc_succi = sum(c_across(all_of(names_succi)) == 1) / length(names_succi) * 100) %>% 
  rowwise() %>%
  mutate(perc_sp = sum(c_across(all_of(names_sp)) == 1) / length(names_sp) * 100) %>% 
  rowwise() %>%
  mutate(perc_faec = sum(c_across(all_of(names_faec)) == 1) / length(names_faec) * 100) 

## counting number of genes for each species
roary_all_long %>% filter(fasta %in% names_succi) %>% dplyr::count(Identifier) %>% dim()
roary_all_long %>% filter(fasta %in% names_sp) %>% dplyr::count(Identifier) %>% dim()
roary_all_long %>% filter(fasta %in% names_faec) %>% dplyr::count(Identifier) %>% dim()

########## kegg id 

# first printing all of kegg ids that are distinct 
all_kegg <- roary_all_long %>%
  select(ko_id, fasta) %>% 
  mutate(presence = 1) %>% 
  filter(ko_id != "") %>% 
  select(ko_id, presence) %>% 
  distinct(ko_id, .keep_all = TRUE) %>% 
  write_tsv("all_kegg.tsv")

### checking all on ccstatus
roary_cc_status <- roary_all_long %>%
  select(ko_id, fasta) %>% 
  mutate(presence = 1) %>% 
  filter(ko_id != "") %>% 
  distinct(ko_id, fasta, .keep_all = TRUE) %>% 
  pivot_wider(id_cols = ko_id, names_from = fasta, values_from = presence, values_fill = 0) %>%  
  rowwise() %>%
  mutate(perc_cont = sum(c_across(all_of(names_cont)) == 1) / length(names_cont) * 100) %>% 
  rowwise() %>%
  mutate(perc_hra = sum(c_across(all_of(names_hra)) == 1) / length(names_hra) * 100) %>% 
  rowwise() %>%
  mutate(perc_crc = sum(c_across(all_of(names_crc)) == 1) / length(names_crc) * 100) 

################### kegg id 
kegg_succi <- roary_all_long %>%
  select(ko_id, fasta) %>% 
  mutate(presence = 1) %>% 
  filter(ko_id != "") %>% 
  distinct(ko_id, fasta, .keep_all = TRUE) %>% 
  pivot_wider(id_cols = ko_id, names_from = fasta, values_from = presence, values_fill = 0) %>% 
  column_to_rownames("ko_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  mutate(species = ifelse(sample_id %in% names_succi, "succi", NA)) %>% 
  mutate(species = ifelse(sample_id %in% names_sp, "other", species)) %>% 
  mutate(species = ifelse(sample_id %in% names_faec, "other", species)) %>% 
  mutate(sample_id = str_replace(sample_id, "S-", "S_"),
         sample_id = str_extract(sample_id, "S_\\d+")) %>% 
  column_to_rownames("sample_id") %>% 
  tbl_summary(by="species") %>% 
  add_p(
    pvalue_fun = function(x) style_pvalue(x, digits=3)) %>% 
  add_stat_label() %>% 
  modify_caption("Table 1. Chi2 test") %>% 
  as_gt() %>% 
  .$`_data` %>% 
  mutate(padj = p.adjust(p.value)) %>% 
  dplyr::rename(kegg_id = variable) %>% 
  select(-c(var_type, var_label, row_type, label, stat_label, test_result, parameter, statistic, alternative))  %>% 
  filter(padj < 0.05) %>% 
  dplyr::rename(succi = stat_2) %>% 
  write_tsv("kegg_succi_chi.tsv")

kegg_faec <- roary_all_long %>%
  select(ko_id, fasta) %>% 
  mutate(presence = 1) %>% 
  filter(ko_id != "") %>% 
  distinct(ko_id, fasta, .keep_all = TRUE) %>% 
  pivot_wider(id_cols = ko_id, names_from = fasta, values_from = presence, values_fill = 0) %>% 
  column_to_rownames("ko_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  mutate(species = ifelse(sample_id %in% names_succi, "other", NA)) %>% 
  mutate(species = ifelse(sample_id %in% names_sp, "other", species)) %>% 
  mutate(species = ifelse(sample_id %in% names_faec, "faec", species)) %>% 
  mutate(sample_id = str_replace(sample_id, "S-", "S_"),
         sample_id = str_extract(sample_id, "S_\\d+")) %>% 
  column_to_rownames("sample_id") %>% 
  tbl_summary(by="species") %>% 
  add_p(
    pvalue_fun = function(x) style_pvalue(x, digits=3)) %>% 
  add_stat_label() %>% 
  modify_caption("Table 1. Chi2 test") %>% 
  as_gt() %>% 
  .$`_data` %>% 
  mutate(padj = p.adjust(p.value)) %>% 
  dplyr::rename(kegg_id = variable) %>% 
  select(-c(var_type, var_label, row_type, label, stat_label, test_result, parameter, statistic, alternative))  %>% 
  filter(padj < 0.05) %>% 
  dplyr::rename(faec = stat_1) %>% 
  write_tsv("kegg_faec_chi.tsv")

kegg_sp <- roary_all_long %>%
  select(ko_id, fasta) %>% 
  mutate(presence = 1) %>% 
  filter(ko_id != "") %>% 
  distinct(ko_id, fasta, .keep_all = TRUE) %>% 
  pivot_wider(id_cols = ko_id, names_from = fasta, values_from = presence, values_fill = 0) %>% 
  column_to_rownames("ko_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  mutate(species = ifelse(sample_id %in% names_succi, "other", NA)) %>% 
  mutate(species = ifelse(sample_id %in% names_sp, "sp", species)) %>% 
  mutate(species = ifelse(sample_id %in% names_faec, "other", species)) %>% 
  mutate(sample_id = str_replace(sample_id, "S-", "S_"),
         sample_id = str_extract(sample_id, "S_\\d+")) %>% 
  column_to_rownames("sample_id") %>% 
  tbl_summary(by="species") %>% 
  add_p(
    pvalue_fun = function(x) style_pvalue(x, digits=3)) %>% 
  add_stat_label() %>% 
  modify_caption("Table 1. Chi2 test") %>% 
  as_gt() %>% 
  .$`_data` %>% 
  mutate(padj = p.adjust(p.value)) %>% 
  dplyr::rename(kegg_id = variable) %>% 
  select(-c(var_type, var_label, row_type, label, stat_label, test_result, parameter, statistic, alternative))  %>% 
  filter(padj < 0.05) %>% 
  dplyr::rename(sp = stat_2) %>% 
  write_tsv("kegg_sp_chi.tsv")


## testing if there is a difference in which kegg id the different ccstatus they have
kegg_ccstatus <- roary_all_long %>%
  select(ko_id, fasta) %>% 
  mutate(presence = 1) %>% 
  filter(ko_id != "") %>% 
  distinct(ko_id, fasta, .keep_all = TRUE) %>% 
  pivot_wider(id_cols = ko_id, names_from = fasta, values_from = presence, values_fill = 0) %>% 
  column_to_rownames("ko_id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  mutate(sample_id = str_replace(sample_id, "S-", "S_"),
         sample_id = str_extract(sample_id, "S_\\d+")) %>% 
  left_join(CRCbiome_metadata %>% select(c(sample_id, cc_status))) %>% 
  column_to_rownames("sample_id") %>% 
  tbl_summary(by="cc_status") %>% 
  add_p(
    pvalue_fun = function(x) style_pvalue(x, digits=3)) %>% 
  add_stat_label() %>% 
  modify_caption("Table 1. Chi2 test") %>% 
  as_gt() %>% 
  .$`_data` %>% 
  mutate(padj = p.adjust(p.value)) %>% 
  filter(padj < 0.05)


################### cazy id 
cazy_sp <- roary_all_long %>%
  dplyr::select(cazy_best_hit, fasta) %>% 
  mutate(presence = 1) %>% 
  filter(cazy_best_hit != "") %>% 
  distinct(cazy_best_hit, fasta, .keep_all = TRUE) %>% 
  pivot_wider(id_cols = cazy_best_hit, names_from = fasta, values_from = presence, values_fill = 0) %>% 
  column_to_rownames("cazy_best_hit") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  mutate(species = ifelse(sample_id %in% names_succi, "other", NA)) %>% 
  mutate(species = ifelse(sample_id %in% names_sp, "sp", species)) %>% 
  mutate(species = ifelse(sample_id %in% names_faec, "other", species)) %>% 
  mutate(sample_id = str_replace(sample_id, "S-", "S_"),
         sample_id = str_extract(sample_id, "S_\\d+")) %>% 
  column_to_rownames("sample_id") %>% 
  tbl_summary(by="species") %>% 
  add_p(
    pvalue_fun = function(x) style_pvalue(x, digits=3)) %>% 
  add_stat_label() %>% 
  modify_caption("Table 1. Chi2 test") %>% 
  as_gt() %>% 
  .$`_data` %>% 
  mutate(padj = p.adjust(p.value)) %>% 
  dplyr::rename(cazy_best_hit = variable) %>% 
  dplyr::select(-c(var_type, var_label, row_type, label, stat_label, test_result, parameter, statistic, alternative))  %>% 
  filter(padj < 0.05) %>% 
  dplyr::rename(sp = stat_2) %>%  ### add and stop here if wanting to write out a file write_tsv("chisq_sp_cazy.tsv")
  dplyr::select(cazy_best_hit, stat_1, sp) %>% 
  mutate(spec = str_extract(sp, "(?<=\\()\\d+\\.?\\d*(?=\\%)"),
         other = str_extract(stat_1, "(?<=\\()\\d+\\.?\\d*(?=\\%)")) %>% 
  dplyr::select(-c(sp, stat_1)) %>% 
  mutate(spec = as.numeric(spec), other = as.numeric(other)) %>% 
  mutate(log2FC = log2((spec + 1) / (other + 1))) %>% 
  mutate(cazy_best_hit = str_remove(cazy_best_hit, ".hmm")) %>% 
  mutate(species = "P. sp")

cazy_succi <- roary_all_long %>%
  dplyr::select(cazy_best_hit, fasta) %>% 
  mutate(presence = 1) %>% 
  filter(cazy_best_hit != "") %>% 
  distinct(cazy_best_hit, fasta, .keep_all = TRUE) %>% 
  pivot_wider(id_cols = cazy_best_hit, names_from = fasta, values_from = presence, values_fill = 0) %>% 
  column_to_rownames("cazy_best_hit") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  mutate(species = ifelse(sample_id %in% names_succi, "succi", NA)) %>% 
  mutate(species = ifelse(sample_id %in% names_sp, "other", species)) %>% 
  mutate(species = ifelse(sample_id %in% names_faec, "other", species)) %>% 
  mutate(sample_id = str_replace(sample_id, "S-", "S_"),
         sample_id = str_extract(sample_id, "S_\\d+")) %>% 
  column_to_rownames("sample_id") %>% 
  tbl_summary(by="species") %>% 
  add_p(
    pvalue_fun = function(x) style_pvalue(x, digits=3)) %>% 
  add_stat_label() %>% 
  modify_caption("Table 1. Chi2 test") %>% 
  as_gt() %>% 
  .$`_data` %>% 
  mutate(padj = p.adjust(p.value)) %>% 
  dplyr::rename(cazy_best_hit = variable) %>% 
  dplyr::select(-c(var_type, var_label, row_type, label, stat_label, test_result, parameter, statistic, alternative))  %>% 
  filter(padj < 0.05) %>% 
  dplyr::rename(succi = stat_2) %>% #write_tsv("cazy_succi_chi_pairwise.tsv") %>% 
  dplyr::select(cazy_best_hit, stat_1, succi) %>% 
  mutate(spec = str_extract(succi, "(?<=\\()\\d+\\.?\\d*(?=\\%)"),
         other = str_extract(stat_1, "(?<=\\()\\d+\\.?\\d*(?=\\%)")) %>% 
  select(-c(succi, stat_1)) %>% 
  mutate(spec = as.numeric(spec), other = as.numeric(other)) %>% 
  mutate(log2FC = log2((spec + 1) / (other + 1))) %>% 
  mutate(cazy_best_hit = str_remove(cazy_best_hit, ".hmm")) %>% 
  mutate(species = "P. succinatutens")

cazy_faec <- roary_all_long %>%
  dplyr::select(cazy_best_hit, fasta) %>% 
  mutate(presence = 1) %>% 
  filter(cazy_best_hit != "") %>% 
  distinct(cazy_best_hit, fasta, .keep_all = TRUE) %>% 
  pivot_wider(id_cols = cazy_best_hit, names_from = fasta, values_from = presence, values_fill = 0) %>% 
  column_to_rownames("cazy_best_hit") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  mutate(species = ifelse(sample_id %in% names_succi, "other", NA)) %>% 
  mutate(species = ifelse(sample_id %in% names_sp, "other", species)) %>% 
  mutate(species = ifelse(sample_id %in% names_faec, "faec", species)) %>% 
  mutate(sample_id = str_replace(sample_id, "S-", "S_"),
         sample_id = str_extract(sample_id, "S_\\d+")) %>% 
  column_to_rownames("sample_id") %>% 
  tbl_summary(by="species") %>% 
  add_p(
    pvalue_fun = function(x) style_pvalue(x, digits=3)) %>% 
  add_stat_label() %>% 
  modify_caption("Table 1. Chi2 test") %>% 
  as_gt() %>% 
  .$`_data` %>% 
  mutate(padj = p.adjust(p.value)) %>% 
  dplyr::rename(cazy_best_hit = variable) %>% 
  dplyr::select(-c(var_type, var_label, row_type, label, stat_label, test_result, parameter, statistic, alternative))  %>% 
  filter(padj < 0.05) %>%  ##  write_tsv("cazy_faec_chi_pairwise.tsv")
  dplyr::rename(faec = stat_1) %>% 
  dplyr::select(cazy_best_hit, stat_2, faec) %>% 
  mutate(spec = str_extract(faec, "(?<=\\()\\d+\\.?\\d*(?=\\%)"),
         other = str_extract(stat_2, "(?<=\\()\\d+\\.?\\d*(?=\\%)")) %>% 
  dplyr::select(-c(faec, stat_2)) %>% 
  mutate(spec = as.numeric(spec), other = as.numeric(other)) %>% 
  mutate(log2FC = log2((spec + 1) / (other + 1))) %>% 
  mutate(cazy_best_hit = str_remove(cazy_best_hit, ".hmm")) %>% 
  mutate(species = "P. faecium")


(p_cazy_log2fc <- cazy_faec %>% 
    rbind(cazy_sp) %>% 
    rbind(cazy_succi) %>% 
    dplyr::select(cazy_best_hit, log2FC, species, spec, other) %>% 
    mutate(species = factor(species, levels=c("P. succinatutens", "P. sp", "P. faecium"))) %>%
    mutate(cazy_best_hit = fct_reorder(cazy_best_hit, log2FC)) %>% 
    mutate(cazy_best_hit = factor(cazy_best_hit, levels = c("CBM13", "GT32", "GT4", "GT8", "GT113", "GT112", "GH171", "GH3", "GH33","GH126"))) %>% 
    ggplot(., aes(x = log2FC, y = cazy_best_hit, fill = species)) +
    geom_errorbarh(aes(xmin = 0, xmax = log2FC, color = species), height = 0, position = position_dodge(width = .5)) +
    geom_point(aes(color = species), position = position_dodge(widt = 0.5), size = 4) + 
    geom_vline(lty = 2, xintercept = 0) +
    xlab("Log2FC") +
    ylab("") +
    scale_color_manual(values = c("P. succinatutens" = "#00A087FF", "P. sp" = "#F39B7FFF", "P. faecium" = "#3C5488FF")) +
    theme_light() +
    theme(axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size = 14),
          axis.title.y = element_text(size = 12),
          axis.text = element_text(size = 12)))

#########################################

##################### uniref ----


## what is unique for faecium
faec_dram <- roary_mat_all %>% 
  filter(perc_faec > 95 & perc_sp == 0 & perc_succi == 0) %>% 
  select(-c(perc_succi, perc_faec, perc_sp)) %>% 
  pivot_longer(cols = -Gene, names_to = "sample", values_to = "presence") %>% 
  filter(presence == 1) %>% 
  inner_join(roary_all_long, by=c("Gene" = "Identifier", "sample" = "fasta")) %>% 
  distinct(Gene, uniref_id, .keep_all = TRUE) %>% 
  select(c(Gene, sample, uniref_id, uniref_hit, uniref_taxonomy, ko_id, kegg_hit, cazy_ids, cazy_hits)) %>% 
  group_by(Gene) %>%
  summarise(
    sample = paste(unique(sample), collapse = ", "),
    uniref_id = paste(unique(uniref_id), collapse = ", "),
    uniref_hit = paste(unique(uniref_hit), collapse = ", "),
    uniref_taxonomy = paste(unique(uniref_taxonomy), collapse = ", "),
    ko_id = paste(unique(ko_id), collapse = ", "),
    kegg_hit = paste(unique(kegg_hit), collapse = ", "),
    cazy_ids = paste(unique(cazy_ids), collapse = ", "),
    cazy_hits = paste(unique(cazy_hits), collapse = ", ")) %>%
  ungroup() %>% 
  select(-c(sample)) %>% 
  write_tsv("faec_dram.tsv")

## what is unique for succi
succi_dram <- roary_mat_all %>% 
  filter(perc_faec == 0 & perc_sp == 0 & perc_succi > 95) %>% 
  select(-c(perc_succi, perc_faec, perc_sp)) %>% 
  pivot_longer(cols = -Gene, names_to = "sample", values_to = "presence") %>% 
  filter(presence == 1) %>% 
  inner_join(roary_all_long, by=c("Gene" = "Identifier", "sample" = "fasta")) %>% 
  distinct(Gene, uniref_id, .keep_all = TRUE) %>% 
  select(c(Gene, sample, uniref_id, uniref_hit, uniref_taxonomy, ko_id, kegg_hit, cazy_ids, cazy_hits)) %>% 
  group_by(Gene) %>%
  summarise(
    sample = paste(unique(sample), collapse = ", "),
    uniref_id = paste(unique(uniref_id), collapse = ", "),
    uniref_hit = paste(unique(uniref_hit), collapse = ", "),
    uniref_taxonomy = paste(unique(uniref_taxonomy), collapse = ", "),
    ko_id = paste(unique(ko_id), collapse = ", "),
    kegg_hit = paste(unique(kegg_hit), collapse = ", "),
    cazy_ids = paste(unique(cazy_ids), collapse = ", "),
    cazy_hits = paste(unique(cazy_hits), collapse = ", ")) %>%
  ungroup() %>% 
  select(-c(sample)) %>% 
  write_tsv("succi_dram.tsv")

## what is unique for sp
sp_dram <- roary_mat_all %>% 
  filter(perc_faec == 0 & perc_sp > 95 & perc_succi == 0) %>% 
  select(-c(perc_succi, perc_faec, perc_sp)) %>% 
  pivot_longer(cols = -Gene, names_to = "sample", values_to = "presence") %>% 
  filter(presence == 1) %>% 
  inner_join(roary_all_long, by=c("Gene" = "Identifier", "sample" = "fasta")) %>% 
  distinct(Gene, uniref_id, .keep_all = TRUE) %>% 
  select(c(Gene, sample, uniref_id, uniref_hit, uniref_taxonomy, ko_id, kegg_hit, cazy_ids, cazy_hits)) %>% 
  group_by(Gene) %>%
  summarise(
    sample = paste(unique(sample), collapse = ", "),
    uniref_id = paste(unique(uniref_id), collapse = ", "),
    uniref_hit = paste(unique(uniref_hit), collapse = ", "),
    uniref_taxonomy = paste(unique(uniref_taxonomy), collapse = ", "),
    ko_id = paste(unique(ko_id), collapse = ", "),
    kegg_hit = paste(unique(kegg_hit), collapse = ", "),
    cazy_ids = paste(unique(cazy_ids), collapse = ", "),
    cazy_hits = paste(unique(cazy_hits), collapse = ", ")) %>%
  ungroup() %>% 
  select(-c(sample)) %>% 
  write_tsv("sp_dram.tsv")

#########################################
##### performing enrichment analyses on KEGG

faecium_higher <- list_of_data_frames[[1]] %>% 
  separate(faec, into = c("faec_nr", "faec_perc"), sep = " \\(", remove = FALSE) %>%
  mutate(faec_perc = gsub("\\%)", "", faec_perc)) %>% 
  separate(stat_2, into = c("other_nr", "other_perc"), sep = " \\(", remove = FALSE) %>%
  mutate(other_perc = gsub("\\%)", "", other_perc)) %>% 
  mutate(other_perc = as.numeric(other_perc), faec_perc = as.numeric(faec_perc)) %>% 
  mutate(result = if_else(faec_perc > other_perc, "higher", "lower")) %>% 
  filter(result == "higher")
kegg_faecium_higher <- as_vector(faecium_higher$kegg_id)
enr_faecium_higher <- enrichKO(kegg_faecium_higher, 
                               universe = all,
                               pvalueCutoff = 10, 
                               pAdjustMethod = "fdr", 
                               minGSSize = 10, 
                               maxGSSize = 500, 
                               qvalueCutoff = 0.05)
summary_faecium_higher <- enr_faecium_higher@result

succi_higher <- list_of_data_frames[[3]] %>% 
  separate(succi, into = c("succi_nr", "succi_perc"), sep = " \\(", remove = FALSE) %>%
  mutate(succi_perc = gsub("\\%)", "", succi_perc)) %>% 
  separate(stat_1, into = c("other_nr", "other_perc"), sep = " \\(", remove = FALSE) %>%
  mutate(other_perc = gsub("\\%)", "", other_perc)) %>% 
  mutate(other_perc = as.numeric(other_perc), succi_perc = as.numeric(succi_perc)) %>% 
  mutate(result = if_else(succi_perc > other_perc, "higher", "lower")) %>% 
  filter(result == "higher")
kegg_succi_higher <- as_vector(succi_higher$kegg_id)

enr_succi_higher <- enrichKO(kegg_succi_higher, 
                             universe = all,
                             pvalueCutoff = 10, 
                             pAdjustMethod = "fdr", 
                             minGSSize = 10, 
                             maxGSSize = 500, 
                             qvalueCutoff = 0.05)
summary_succi_higher <- enr_succi_higher@result


sp_higher <- list_of_data_frames[[2]] %>% 
  separate(sp, into = c("sp_nr", "sp_perc"), sep = " \\(", remove = FALSE) %>%
  mutate(sp_perc = gsub("\\%)", "", sp_perc)) %>% 
  separate(stat_1, into = c("other_nr", "other_perc"), sep = " \\(", remove = FALSE) %>%
  mutate(other_perc = gsub("\\%)", "", other_perc)) %>% 
  mutate(other_perc = as.numeric(other_perc), sp_perc = as.numeric(sp_perc)) %>% 
  mutate(result = if_else(sp_perc > other_perc, "higher", "lower")) %>% 
  filter(result == "higher")
kegg_sp_higher <- as_vector(sp_higher$kegg_id)
enr_sp_higher <- enrichKO(kegg_sp_higher, 
                          universe = all,
                          pvalueCutoff = 10, 
                          pAdjustMethod = "BH", 
                          minGSSize = 10, 
                          maxGSSize = 500, 
                          qvalueCutoff = 0.05)
summary_sp_higher <- enr_sp_higher@result


## Making plot for all combined
#### combining plots
all_summary <- summary_sign_faecium %>% mutate(dataset = "faecium") %>% 
  rbind(summary_sign_sp %>% mutate(dataset = "sp")) %>% 
  rbind(summary_sign_succi %>% mutate(dataset = "succi")) %>% 
  separate(GeneRatio, into = c("gene_set", "gene_all"), sep = "/", remove = FALSE) %>% 
  mutate(gene_set = as.numeric(gene_set), gene_all = as.numeric(gene_all)) %>% 
  mutate(generatio = gene_set/gene_all) %>% 
  filter(p.adjust < 0.05) %>% 
  mutate(Count = as.numeric(Count)) %>% 
  mutate(Description = factor(Description, levels = c("Phosphonate and phosphinate metabolism", "Histidine metabolism", "Two-component system", "ABC transporters", 
                                                      "Biosynthesis of amino acids", "Glycine, serine and threonine metabolism", 
                                                      "Glyoxylate and dicarboxylate metabolism","Porphyrin and chlorophyll metabolism")))
pdf("ORA.pdf", width = "7", height = "5")
p1 <- all_summary %>% ggplot(aes(generatio, Description, colour=dataset, size=Count)) +
  geom_point() +
  scale_size_continuous(range = c(3,7)) +
  scale_color_manual(values = c("sp" = "#F39B7FFF", "succi"="#00A087FF",  "faecium" = "#3C5488FF")) +
  xlim(0.047, 0.195) +
  theme_bw() +
  theme(axis.text = element_text(size=12))
dev.off()
