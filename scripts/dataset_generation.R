############### importing data
############################################

## NORCCAP metagenome
norccap_MG <- read_delim("norccap_MG") %>% 
  select(-c(clade_taxid)) %>%
  mutate(clade_name = gsub("[|]", "_", clade_name)) %>% 
  filter(grepl("s__", clade_name)) %>% 
  mutate(clade_name=sapply(strsplit(clade_name, split="s__"), function(x) x[2])) %>%
  column_to_rownames("clade_name") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  mutate(sample_id = gsub("_mp_profile", "", sample_id)) %>% 
  filter(sample_id!="unmatched" & sample_id!="Blank_NORCCAP" & sample_id!="ZYMO_std") %>% 
  as.data.frame() %>% 
  select(c(sample_id, Phascolarctobacterium_succinatutens, Phascolarctobacterium_sp_CAG_266, Phascolarctobacterium_faecium)) %>% 
  mutate(pres_succinatutens = ifelse(Phascolarctobacterium_succinatutens > 0, "yes", "no")) %>% 
  mutate(pres_faecium = ifelse(Phascolarctobacterium_faecium > 0, "yes", "no")) %>% 
  mutate(pres_sp = ifelse(Phascolarctobacterium_sp_CAG_266 > 0, "yes", "no")) %>% 
  mutate(log_succinatutens = LOG(Phascolarctobacterium_succinatutens)) %>% 
  mutate(log_faecium = LOG(Phascolarctobacterium_faecium)) %>% 
  mutate(log_sp = LOG(Phascolarctobacterium_sp_CAG_266)) %>% 
  left_join(metadata_MG %>% select(c(sample_id, cc_status, sex, age, senterid))) %>% 
  dplyr::rename(geography = senterid) %>% 
  mutate(cc_status = case_when(
    cc_status == "1" ~ "Control",
    cc_status == "2" ~ "Adenoma",
    cc_status == "3" ~ "CRC",
    TRUE ~ as.character(cc_status))) %>%
  mutate(sex = case_when(
    sex == "0" ~ "Female",
    sex == "1" ~ "Male",
    TRUE ~ as.character(sex))) %>% 
  mutate(cc_status=factor(cc_status, levels=c("Control", "Adenoma", "CRC"))) %>% 
  mutate(age = as.numeric(age), geography = as.factor(geography), sex = as.factor(sex)) %>% 
  as.data.frame() %>% 
  mutate(phasco_group = case_when(
    pres_succinatutens == "yes" & pres_faecium == "no" &  pres_sp == "no" ~ "succ",
    pres_faecium == "yes" & pres_sp == "no" &  pres_succinatutens == "no" ~ "faec",
    pres_sp == "yes" & pres_succinatutens == "no" &  pres_faecium == "no" ~ "sp",
    pres_succinatutens == "yes" & pres_faecium == "yes" | pres_succinatutens == "yes" & pres_sp == "yes" | pres_sp == "yes" & pres_faecium == "yes" ~ "2grp",
    TRUE ~"none"
  ))

norccap_maaslin_MG <- read_delim("norccap_maaslin_MG") %>% 
  select(-c(clade_taxid)) %>%
  mutate(clade_name = gsub("[|]", "_", clade_name)) %>% 
  filter(grepl("s__", clade_name)) %>% 
  mutate(clade_name=sapply(strsplit(clade_name, split="s__"), function(x) x[2])) %>%
  column_to_rownames("clade_name") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  mutate(sample_id = gsub("_mp_profile", "", sample_id)) %>% 
  filter(sample_id!="unmatched" & sample_id!="Blank_NORCCAP" & sample_id!="ZYMO_std") %>% 
  as.data.frame() 

norccap_maaslin_MG_metadata <- metadata_MG %>% 
  select(c(sample_id, cc_status, sex, age, senterid)) %>% 
  dplyr::rename(geography = senterid) %>% 
  mutate(cc_status = case_when(
    cc_status == "1" ~ "Control",
    cc_status == "2" ~ "Adenoma",
    cc_status == "3" ~ "CRC",
    TRUE ~ as.character(cc_status))) %>%
  mutate(sex = case_when(
    sex == "0" ~ "Female",
    sex == "1" ~ "Male",
    TRUE ~ as.character(sex))) %>% 
  mutate(cc_status=factor(cc_status, levels=c("Control", "Adenoma", "CRC"))) %>% 
  mutate(age = as.numeric(age), geography = as.factor(geography), sex = as.factor(sex)) %>% 
  left_join(norccap_MG %>% select(c(pres_sp, pres_succinatutens, pres_faecium, sample_id))) %>% 
  mutate(phasco_group = case_when(
    pres_succinatutens == "yes" & pres_faecium == "no" &  pres_sp == "no" ~ "succ",
    pres_faecium == "yes" & pres_sp == "no" &  pres_succinatutens == "no" ~ "faec",
    pres_sp == "yes" & pres_succinatutens == "no" &  pres_faecium == "no" ~ "sp",
    pres_succinatutens == "yes" & pres_faecium == "yes" | pres_succinatutens == "yes" & pres_sp == "yes" | pres_sp == "yes" & pres_faecium == "yes" ~ "2grp",
    TRUE ~"none"
  )) %>% 
  mutate(phasco_group = factor(phasco_group, levels=c("none", "succ", "sp", "faec", "2grp")))


## NORCCAP ASV
norccap_ASV <- physeq_ASV %>%
  transform_sample_counts(., function(x) x / sum(x)*100) %>%
  subset_taxa(Genus=="Phascolarctobacterium") %>% 
  otu_table() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  rowwise() %>% 
  mutate(Phascolarctobacterium_succinatutens = sum(c_across(c(2:25, 27:39)))) %>%
  mutate(Phascolarctobacterium_faecium = sum(fb7a1c0b3625f00ee42dbfcbaa001f12, a3e7c20c1249c83ca701a4a18d5d58f6)) %>%
  dplyr::rename(Phascolarctobacterium_sp_CAG_266 = "524e9ca6a89d33deb85f0c2ba1ade0ff") %>% 
  as.data.frame() %>% 
  left_join(physeq_ASV %>% sample_data() %>% as.matrix() %>% as.data.frame() %>% rownames_to_column(var="sample_id")) %>% 
  select(c(sample_id, Phascolarctobacterium_succinatutens, Phascolarctobacterium_faecium, Phascolarctobacterium_sp_CAG_266, cc_status, sex, age, senterid)) %>% 
  mutate(pres_succinatutens = ifelse(Phascolarctobacterium_succinatutens > 0, "yes", "no")) %>% 
  mutate(pres_faecium = ifelse(Phascolarctobacterium_faecium > 0, "yes", "no")) %>% 
  mutate(pres_sp = ifelse(Phascolarctobacterium_sp_CAG_266 > 0, "yes", "no")) %>% 
  mutate(log_succinatutens = LOG(Phascolarctobacterium_succinatutens)) %>% 
  mutate(log_faecium = LOG(Phascolarctobacterium_faecium)) %>% 
  mutate(log_sp = LOG(Phascolarctobacterium_sp_CAG_266)) %>% 
  dplyr::rename(geography = senterid) %>% 
  mutate(cc_status = case_when(
    cc_status == "1" ~ "Control",
    cc_status == "2" ~ "Adenoma",
    cc_status == "3" ~ "CRC",
    TRUE ~ as.character(cc_status))) %>%
  mutate(sex = case_when(
    sex == "0" ~ "Female",
    sex == "1" ~ "Male",
    TRUE ~ as.character(sex))) %>% 
  mutate(cc_status=factor(cc_status, levels=c("Control", "Adenoma", "CRC"))) %>% 
  mutate(age = as.numeric(age), geography = as.factor(geography), sex = as.factor(sex)) %>% 
  as.data.frame() %>% 
  mutate(phasco_group = case_when(
    pres_succinatutens == "yes" & pres_faecium == "no" &  pres_sp == "no" ~ "succ",
    pres_faecium == "yes" & pres_sp == "no" &  pres_succinatutens == "no" ~ "faec",
    pres_sp == "yes" & pres_succinatutens == "no" &  pres_faecium == "no" ~ "sp",
    pres_succinatutens == "yes" & pres_faecium == "yes" | pres_succinatutens == "yes" & pres_sp == "yes" | pres_sp == "yes" & pres_faecium == "yes" ~ "2grp",
    TRUE ~"none"
  ))

## collapse all ASVs belonging to p. succinatutens based on tree proximity
succinatutens_ASVs <- physeq_ASV %>%
  subset_taxa(Genus=="Phascolarctobacterium") %>% 
  otu_table() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  select(all_of(c(2:25, 27:39))) %>% 
  colnames() %>% 
  unlist()

norccap_maaslin_ASV <- physeq_ASV %>%
  transform_sample_counts(., function(x) x / sum(x)*100) %>%
  otu_table() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>%  
  mutate(Phascolarctobacterium_faecium = sum(fb7a1c0b3625f00ee42dbfcbaa001f12, a3e7c20c1249c83ca701a4a18d5d58f6)) %>% 
  mutate(Phascolarctobacterium_succinatutens = rowSums(select(., all_of(succinatutens_ASVs)))) %>%  
  select(-c(all_of(succinatutens_ASVs))) %>% 
  select(-c(fb7a1c0b3625f00ee42dbfcbaa001f12, a3e7c20c1249c83ca701a4a18d5d58f6)) %>% 
  dplyr::rename(Phascolarctobacterium_sp_CAG_266 = "524e9ca6a89d33deb85f0c2ba1ade0ff") %>% 
  as.data.frame()

norccap_maaslin_ASV_metadata <- physeq_ASV %>% 
  sample_data() %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column(var="sample_id") %>% 
  select(c(sample_id, cc_status, sex, age, senterid)) %>% 
  dplyr::rename(geography = senterid) %>% 
  mutate(cc_status = case_when(
    cc_status == "1" ~ "Control",
    cc_status == "2" ~ "Adenoma",
    cc_status == "3" ~ "CRC",
    TRUE ~ as.character(cc_status))) %>%
  mutate(sex = case_when(
    sex == "0" ~ "Female",
    sex == "1" ~ "Male",
    TRUE ~ as.character(sex))) %>% 
  mutate(cc_status=factor(cc_status, levels=c("Control", "Adenoma", "CRC"))) %>% 
  mutate(age = as.numeric(age), geography = as.factor(geography), sex = as.factor(sex)) %>% 
  as.data.frame() %>% 
  left_join(norccap_ASV %>% select(c(pres_sp, pres_succinatutens, pres_faecium, sample_id))) %>% 
  mutate(phasco_group = case_when(
    pres_succinatutens == "yes" & pres_faecium == "no" &  pres_sp == "no" ~ "succ",
    pres_faecium == "yes" & pres_sp == "no" &  pres_succinatutens == "no" ~ "faec",
    pres_sp == "yes" & pres_succinatutens == "no" &  pres_faecium == "no" ~ "sp",
    pres_succinatutens == "yes" & pres_faecium == "yes" | pres_succinatutens == "yes" & pres_sp == "yes" | pres_sp == "yes" & pres_faecium == "yes" ~ "2grp",
    TRUE ~"none"
  )) %>% 
  mutate(phasco_group = factor(phasco_group, levels=c("none", "succ", "sp", "faec", "2grp")))

## CRCbiome metadata

CRCbiome_metadata <- sample_meta %>% 
  left_join(screening_data) %>% 
  mutate(final_result_cat7 = as.character(final_result_cat7)) %>% 
  mutate(final_result_cat4 = as.character(final_result_cat4)) %>% 
  left_join(lifestyle) %>% 
  left_join(diet, by=c("ffq_ref_nr" = "id")) %>% 
  left_join(wcrf, by=c("ffq_ref_nr" = "id")) %>% 
  mutate(Cur_smoking_status_3cat = as.character(Cur_smoking_status_3cat)) %>% 
  mutate(smoking2 = case_when(Cur_smoking_status_3cat %in% c("Smoked earlier", "Current smoker") ~ "Ever smoked", TRUE ~ Cur_smoking_status_3cat)) %>% 
  mutate(final_result_cat5 = case_when(final_result_cat7 %in% c("Negative", "Non-advanced serrated/other lesion", "Non-advanced adenoma (<3)") ~ "Negative", TRUE ~ final_result_cat7)) %>% 
  filter(., Prøvetype == "Baseline" & (reads_proc*2*151) >= 1e9) %>% ## subset on qc_atlas over 1 milliard og baseline  
  mutate(wcrf_index_main_a_clincut4 = as.factor(as.character(wcrf_index_main_a_clincut4))) %>% 
  mutate(localization = as.character(localization)) %>% 
  mutate(localization = ifelse(is.na(localization), "Negative", localization)) %>% 
  mutate(location2 = case_when(localization %in% c("Proximal", "Both") ~ "any Proximal", TRUE ~ localization)) %>%
  mutate(Bowel_disorder_merged = as.character(Bowel_disorder_merged)) %>% 
  mutate(bowel_disorder2 = case_when(Bowel_disorder_merged %in% c("Other", "IBD", "IBS", "Celiac disease") ~ "Any disorder", TRUE ~ Bowel_disorder_merged)) %>% 
  mutate(across(c("Tarmkreft_Familie", "smoking2", "Antibiotics", "bowel_disorder2", "Intolerance"), ~ifelse(. %in% c("Unknown", "Missing"), NA_character_, as.character(.)))) %>% 
  mutate(wcrf_index_main_a = if_else(wcrf_index_main_a %in% c("Unknown", "Missing"), NA_real_, wcrf_index_main_a)) %>% 
  mutate(shannon_tertile = cut_number(shannon, n=3, labels=c("Low", "Medium", "High"))) %>% 
  mutate(invsimp_tertile = cut_number(invsimpson, n=3, labels=c("Low", "Medium", "High"))) %>% 
  inner_join(metaphlan_species %>% select(sample_id)) %>% 
  select(c(sample_id, kjonn, age_invitation, senter, final_result_cat3_serr, final_result_cat7, detect_worthy_lesions)) %>% 
  mutate(final_result_cat5 = case_when(final_result_cat7 %in% c("Negative", "Non-advanced serrated/other lesion", "Non-advanced adenoma (<3)", 
                                                             "Non-advanced adenoma (>=3)") ~ "Negative", TRUE ~ final_result_cat7)) %>% 
  mutate(cc_status = case_when(
    final_result_cat5 %in% c("Advanced adenoma", "Advanced serrated") ~ "Advanced adenoma", TRUE ~ final_result_cat5)) %>% 
  mutate(cc_status = case_when(
    cc_status == "Negative" ~ "Control",
    cc_status == "Advanced adenoma" ~ "Adenoma",
    cc_status == "CRC" ~ "CRC",
    TRUE ~ as.factor(cc_status))) %>%
  mutate(sex = case_when(
    kjonn == "Female" ~ "Female",
    kjonn == "Male" ~ "Male",
    TRUE ~ as.factor(kjonn))) %>% 
  mutate(cc_status=factor(cc_status, levels=c("Control", "Adenoma", "CRC"))) %>%
  mutate(final_result_cat7=factor(final_result_cat7, levels=c("Negative","Non-advanced serrated/other lesion", "Non-advanced adenoma (>=3)", "Advanced serrated", "Advanced adenoma", "CRC"))) %>%
  select(-c(kjonn, final_result_cat5, final_result_cat3_serr)) %>% 
  dplyr::rename(geography = senter) %>% 
  dplyr::rename(age = age_invitation) %>% 
  mutate(age = as.numeric(age), geography = as.factor(geography), sex = as.factor(sex)) %>% 
  as.data.frame() %>% 
  left_join(CRCbiome_mags_MG %>% select(c(pres_sp, pres_succinatutens, pres_faecium, sample_id))) %>% 
  mutate(phasco_group = case_when(
    pres_succinatutens == "yes" & pres_faecium == "no" &  pres_sp == "no" ~ "succ",
    pres_faecium == "yes" & pres_sp == "no" &  pres_succinatutens == "no" ~ "faec",
    pres_sp == "yes" & pres_succinatutens == "no" &  pres_faecium == "no" ~ "sp",
    pres_succinatutens == "yes" & pres_faecium == "yes" | pres_succinatutens == "yes" & pres_sp == "yes" | pres_sp == "yes" & pres_faecium == "yes" ~ "2grp",
    TRUE ~"none"
  )) %>% 
  mutate(phasco_group = factor(phasco_group, levels=c("none", "succ", "sp", "faec", "2grp")))


## CRCbiome mags
# first extracting those MAG_ID that belong to phascolarctobacterium
mag_taxonomy %>% filter(grepl("Phascolarctobacterium", genus)) %>% select(MAG_id)

CRCbiome_mags_MG <- readRDS("MAG_abundance.Rds") %>% 
  select(c("sample_id", "MAG1658", "MAG1826", "MAG1177", "MAG0822", "MAG0258")) %>% 
  left_join(CRCbiome_metadata) %>% 
  dplyr::rename(Phascolarctobacterium_faecium = MAG1658) %>% 
  dplyr::rename(Phascolarctobacterium_sp_CAG_266 = MAG1826) %>%
  mutate(pres_succinatutens = ifelse(MAG0258 > 0 | MAG0822 > 0 | MAG1177 > 0, "yes", "no")) %>% 
  mutate(pres_faecium = ifelse(Phascolarctobacterium_faecium > 0, "yes", "no")) %>% 
  mutate(pres_sp = ifelse(Phascolarctobacterium_sp_CAG_266 > 0, "yes", "no")) %>% 
  mutate(log_mag1177 = LOG(MAG1177)) %>%
  mutate(log_mag0258 = LOG(MAG0258)) %>%
  mutate(log_mag0822 = LOG(MAG0822)) %>%
  rowwise() %>% 
  mutate(Phascolarctobacterium_succinatutens = sum(MAG0822, MAG1177, MAG0258)) %>% 
  as.data.frame() %>% 
  mutate(log_succinatutens = LOG(Phascolarctobacterium_succinatutens)) %>% 
  mutate(log_faecium = LOG(Phascolarctobacterium_faecium)) %>%
  mutate(log_sp = LOG(Phascolarctobacterium_sp_CAG_266)) %>% 
  as.data.frame()

CRCbiome_mags_maaslin_MG <- readRDS("MAG_abundance.Rds") %>% 
  dplyr::rename(Phascolarctobacterium_faecium = MAG1658) %>% 
  dplyr::rename(Phascolarctobacterium_sp_CAG_266 = MAG1826) %>%
  rowwise() %>% 
  mutate(Phascolarctobacterium_succinatutens = sum(MAG0822, MAG1177, MAG0258)) %>% 
  as.data.frame() %>% 
  select(-c(MAG0822, MAG1177, MAG0258))%>% 
  as.data.frame()


## CRCbiome meta
CRCbiome_meta_MG <- readRDS("metaphlan_species.Rds") %>% 
  rename_with(~ str_remove(., "s__"), starts_with("s__")) %>% 
  select(c(sample_id, Phascolarctobacterium_succinatutens, Phascolarctobacterium_faecium, Phascolarctobacterium_sp_CAG_266)) %>% 
  left_join(CRCbiome_metadata) %>% 
  mutate(pres_succinatutens = ifelse(Phascolarctobacterium_succinatutens > 0, "yes", "no")) %>% 
  mutate(pres_faecium = ifelse(Phascolarctobacterium_faecium > 0, "yes", "no")) %>% 
  mutate(pres_sp = ifelse(Phascolarctobacterium_sp_CAG_266 > 0, "yes", "no")) %>% 
  mutate(log_succinatutens = LOG(Phascolarctobacterium_succinatutens)) %>%
  mutate(log_faecium = LOG(Phascolarctobacterium_faecium)) %>%
  mutate(log_sp = LOG(Phascolarctobacterium_sp_CAG_266)) %>% 
  as.data.frame()

CRCbiome_meta_maaslin_MG <- metaphlan_species %>%
  rename_with(~ str_remove(., "s__"), starts_with("s__"))


## curatedmetagenome
curatedMG_data_rel <- read.delim("20220322_mg_crc_count_data_rel.tsv")
curatedMG_sample_data <- read.delim("20220322_mg_crc_sample_data.tsv")
curatedMG_samples_caco <- curatedMG_sample_data %>% 
  filter(study_condition == "CRC" | study_condition == "control") %>% 
  column_to_rownames("sample_id")

curated_maaslin_MG_metadata <- curatedMG_sample_data %>% 
  select(c(sample_id, gender, country, study_condition, age, study_name)) %>% 
  filter(!is.na(gender), !is.na(country)) %>% 
  dplyr::rename(sex = gender) %>% 
  dplyr::rename(cc_status = study_condition) %>% 
  dplyr::rename(geography = study_name) %>% 
  filter(!is.na(sex), !is.na(country)) %>% 
  inner_join(curated_MG %>% select("sample_id")) %>% 
  mutate(cc_status = case_when(
    cc_status == "control" ~ "Control",
    cc_status == "adenoma" ~ "Adenoma",
    cc_status == "CRC" ~ "CRC",
    TRUE ~ as.factor(cc_status))) %>%
  mutate(sex = case_when(
    sex == "female" ~ "Female",
    sex == "male" ~ "Male",
    TRUE ~ as.character(sex))) %>% 
  mutate(cc_status=factor(cc_status, levels=c("Control", "Adenoma", "CRC"))) %>% 
  mutate(sex = as.factor(sex), age = as.numeric(age), geography = as.factor(geography)) %>% 
  as.data.frame() %>% 
  left_join(curated_MG %>% select(c(pres_sp, pres_succinatutens, pres_faecium, sample_id))) %>% 
  mutate(phasco_group = case_when(
    pres_succinatutens == "yes" & pres_faecium == "no" &  pres_sp == "no" ~ "succ",
    pres_faecium == "yes" & pres_sp == "no" &  pres_succinatutens == "no" ~ "faec",
    pres_sp == "yes" & pres_succinatutens == "no" &  pres_faecium == "no" ~ "sp",
    pres_succinatutens == "yes" & pres_faecium == "yes" | pres_succinatutens == "yes" & pres_sp == "yes" | pres_sp == "yes" & pres_faecium == "yes" ~ "2grp",
    TRUE ~"none"
  )) %>% 
  mutate(phasco_group = factor(phasco_group, levels=c("none", "succ", "sp", "faec", "2grp")))

curated_MG <- curatedMG_data_count %>%
  tidyr::separate(taxon, into= c("name1", "species"), sep="\\|s__") %>%
  filter(startsWith(name1, "k__Bacteria")) %>% 
  dplyr::select(-c("name1")) %>% 
  column_to_rownames("species") %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  select(sample_id, Phascolarctobacterium_faecium, Phascolarctobacterium_succinatutens, Phascolarctobacterium_sp_CAG_266) %>% 
  inner_join(curated_maaslin_MG_metadata) %>% 
  mutate(pres_succinatutens = ifelse(Phascolarctobacterium_succinatutens > 0, "yes", "no")) %>% 
  mutate(pres_faecium = ifelse(Phascolarctobacterium_faecium > 0, "yes", "no")) %>% 
  mutate(pres_sp = ifelse(Phascolarctobacterium_sp_CAG_266 > 0, "yes", "no")) %>% 
  mutate(log_succinatutens = LOG(Phascolarctobacterium_succinatutens)) %>% 
  mutate(log_faecium = LOG(Phascolarctobacterium_faecium)) %>% 
  mutate(log_sp = LOG(Phascolarctobacterium_sp_CAG_266)) %>% 
  filter(!is.na(sex), !is.na(geography)) %>% 
  mutate(cc_status=factor(cc_status, levels=c("Control", "Adenoma", "CRC"))) %>% 
  mutate(sex = as.factor(sex), age = as.numeric(age), geography = as.factor(geography)) %>% 
  as.data.frame()

curated_maaslin_MG <- curatedMG_data_count %>%
  tidyr::separate(taxon, into= c("name1", "species"), sep="\\|s__") %>%
  filter(startsWith(name1, "k__Bacteria")) %>% 
  dplyr::select(-c("name1")) %>% 
  column_to_rownames("species") %>% 
  t() %>%
  as.data.frame() %>% 
  rownames_to_column("sample_id") %>% 
  inner_join(curated_maaslin_MG_metadata %>% filter(!is.na(sex), !is.na(geography)) %>% select(sample_id)) %>% 
  as.data.frame()

############################################