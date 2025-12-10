# 16S
sample_data(exp1_ps_16s_r)$P_Day <- with(sample_data(exp1_ps_16s_r), paste(P.level, Day, sep = "_"))

Merge_by_TU.phylo <- merge_samples(exp1_ps_16s_r, "P_Day")

phylum_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Phylum"), add_meta_data = F)
class_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Class"), add_meta_data = F)
order_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Order"), add_meta_data = F)
family_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Family"), add_meta_data = F)
genus_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Genus"), add_meta_data = F)
species_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Species"), add_meta_data = F)


unique_counts <- data.frame(
  Sample = "Total",
  Phylum = length(unique(phylum_otu_counts$Phylum)),
  Class = length(unique(class_otu_counts$Class)),
  Order = length(unique(order_otu_counts$Order)),
  Family = length(unique(family_otu_counts$Family)),
  Genus = length(unique(genus_otu_counts$Genus)),
  Species = length(unique(species_otu_counts$Species))
)


taxa_levels_summary<-phylum_otu_counts %>%
  dplyr::count(levels(Phylum),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Phylum")

taxa_levels_summary<-class_otu_counts %>%
  dplyr::count(levels(Class),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Class")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-order_otu_counts %>%
  dplyr::count(levels(Order),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Order")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-family_otu_counts %>%
  dplyr::count(levels(Family),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Family")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-genus_otu_counts %>%
  dplyr::count(levels(Genus),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Genus")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-species_otu_counts %>%
  dplyr::count(levels(Species),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Species")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary <- bind_rows(taxa_levels_summary, unique_counts)


taxa_levels_summary_16S <- taxa_levels_summary %>%
  mutate(Community = "Prokaryotes") %>%
  relocate(Community, .before = Sample) %>%
  separate(Sample, into = c("SP (kg/ha)", "Time (d)"), sep = "_")



# ITS
sample_data(exp1_ps_ITS_r)$P_Day <- with(sample_data(exp1_ps_ITS_r), paste(P.level, Day, sep = "_"))
Merge_by_TU.phylo <- merge_samples(exp1_ps_ITS_r, "P_Day")

phylum_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Phylum"), add_meta_data = F)
class_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Class"), add_meta_data = F)
order_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Order"), add_meta_data = F)
family_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Family"), add_meta_data = F)
genus_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Genus"), add_meta_data = F)
species_otu_counts <- phyloseq_ntaxa_by_tax(Merge_by_TU.phylo, TaxRank = c("Species"), add_meta_data = F)


unique_counts <- data.frame(
  Sample = "Total",
  Phylum = length(unique(phylum_otu_counts$Phylum)),
  Class = length(unique(class_otu_counts$Class)),
  Order = length(unique(order_otu_counts$Order)),
  Family = length(unique(family_otu_counts$Family)),
  Genus = length(unique(genus_otu_counts$Genus)),
  Species = length(unique(species_otu_counts$Species))
)


taxa_levels_summary<-phylum_otu_counts %>%
  dplyr::count(levels(Phylum),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Phylum")

taxa_levels_summary<-class_otu_counts %>%
  dplyr::count(levels(Class),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Class")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-order_otu_counts %>%
  dplyr::count(levels(Order),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Order")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-family_otu_counts %>%
  dplyr::count(levels(Family),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Family")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-genus_otu_counts %>%
  dplyr::count(levels(Genus),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Genus")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary<-species_otu_counts %>%
  dplyr::count(levels(Species),Sample,sort = FALSE)%>%
  rename_with(.cols = 2, ~"Species")%>%
  merge(x=taxa_levels_summary,by="Sample")

taxa_levels_summary <- bind_rows(taxa_levels_summary, unique_counts)

# make a nice summary table ITS
taxa_levels_summary_ITS <- taxa_levels_summary %>%
  mutate(Community = "Fungi") %>%
  relocate(Community, .before = Sample) %>%
  separate(Sample, into = c("SP (kg/ha)", "Time (d)"), sep = "_")


all_taxa <- bind_rows(taxa_levels_summary_16S, taxa_levels_summary_ITS)