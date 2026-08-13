library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

# Load in orthogroups matrix
orthos <- readr::read_tsv("../data/gene_CNV_analysis/orthoFinder/Orthogroups.tsv")

# Load N2 gene names and seq IDs
n2_names <- readr::read_tsv("../data/gene_CNV_analysis/gene_names/N2_geneName_seqName.tsv", col_names = c("gene_name", "seq_name"))

# Load in N2 genes of interest
n2_goi <- readr::read_tsv("../data/gene_CNV_analysis/gene_names/GENE_LIST.tsv", col_names = "gene_name") %>% dplyr::pull()

test <- n2_names %>% dplyr::filter(gene_name %in% n2_goi)

# Decompress N2 orthogroups, substitute N2 seq IDs for names, and filter for genes of interest
orthos_filt <- orthos %>%
  dplyr::rename_with(~ sub("\\..*", "", .x)) %>% dplyr::select(1,6,everything()) %>%
  tidyr::separate_rows(N2, sep = ", ") %>%
  dplyr::mutate(N2 = gsub("transcript_","", N2)) %>%
  dplyr::mutate(N2 = sub("\\.[^.]*$", "", N2)) %>% # removing trailing isoform numbers
  dplyr::mutate(N2 = sub("[A-Za-z]+$", "", N2)) %>% # removing trailing letters
  dplyr::left_join(n2_names, by = c("N2" = "seq_name")) %>%
  dplyr::filter(gene_name %in% n2_goi) %>%
  dplyr::group_by(Orthogroup) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ paste(unique(.x), collapse = ", ")), .groups = "drop") %>%
  dplyr::select(-N2) %>%
  dplyr::select(Orthogroup, N2 = gene_name, everything())

# Add CGC2 gene names
cgc2_names <- readr::read_tsv("../data/gene_CNV_analysis/gene_names/CGC2_seqName_geneName.tsv", col_names = c("seq_name", "gene_name"))

# Final orthogroups matrix with gene names
orthos_filt_CGC2_gene_names <- orthos_filt %>%
  dplyr::mutate(CGC2 = gsub("transcript_", "", CGC2)) %>%
  dplyr::left_join(cgc2_names, by = c("CGC2" = "seq_name")) %>%
  dplyr::select(-CGC2) %>%
  dplyr::select(Orthogroup, N2, CGC2 = gene_name, everything())

# Save the results as a table:
write.table(orthos_filt_CGC2_gene_names, "../tables/CNV_genesOfInterest.tsv", quote = F, row.names = F, col.names = T, sep = '\t')