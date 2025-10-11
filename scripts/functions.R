
# function to read .cover
# all files within a directory will be loaded
cover_parser <- function(path = "../coverage", patt = ".cover"){
  
  cols <- c("Sample", "contig", "startpos", "endpos", 
            "numreads","covbases",
            "coverage", "meandepth", "meanbaseq", "meanmapq")
  
  data.frame(sample.id = paste(dir(path, 
                                   pattern = patt))) %>%
    # read each file from the directory (current path) and then unnest
    mutate(file_contents = map(sample.id, ~ read_tsv(
      file.path(path, .), 
      col_names = F))) %>%
    unnest() %>%
    mutate(sample.id = gsub(patt, "", sample.id)) %>%
    filter(X2 != "startpos") -> df
  colnames(df) <- cols
  
  return(df)
}


# function to select donor-specific features
DonFeature <- function(feature_c = "../data/clean/binCover.txt", 
                       donor = "DonA",
                       cutoff_pres = 0,
                       method = "Assembly") {
  
  # feature_c= a long format dataframe of feature coverage
  #            including Sample, feature, and coverage columns
  # donor = donor_name/ label: options: "DonA" or "DonB" 
  #         by default DonA is selected
  # cutoff_pres = feature must covered in donor sample >= of cutoff_pres
  # method = Reads or Asssemnly: in Assembly method, filtering to feature 
  #                             directly assembled from the donor.
  #                              in Reads method just needs abudance cutoff
  
  # Read the data
  df <- read.delim(feature_c, header = TRUE, stringsAsFactors = FALSE)
  
  # Determine donor sample
  donor_sample <- switch(donor,
                         "DonA" = "X608",
                         "DonB" = "Wild116",
                         stop("Invalid donor specified. Choose 'DonA' or 'DonB'."))
  
  # Filter based on method
  if (method == "Assembly") {
    df_filtered <- df %>%
      filter(str_detect(feature, donor)) %>%
      filter(Sample == donor_sample & coverage >= cutoff_pres)
  } else if (method == "Reads") {
    df_filtered <- df %>%
      filter(Sample == donor_sample & coverage >= cutoff_pres)
  } else {
    stop("Invalid method specified. Choose 'Assembly' or 'Reads'.")
  }
  
  # get the list of features that passed the filter
  f_feature <- df_filtered %>% pull(feature) %>% unique()
  
  # pull those features across all samples
  donor_feature <- df %>% filter(feature %in% f_feature)
  
  return(donor_feature)
}


# a function to count the number of colonization in each sample
countColoniz <- function(donor_feature,
                          cutoff_pres = 0,
                          mapfile = "../data/clean/mapfile.txt"){
  
  # Read mapfile
  map <- read.delim(mapfile, header = TRUE, 
                        stringsAsFactors = FALSE)
  
  coloniz_count = donor_feature %>%
    filter(coverage >= cutoff_pres) %>%
    group_by(Sample) %>%
    summarise(n_coloniz = n(), .groups = "drop") %>%
    right_join(map %>% select(Sample, Group, Category), by = "Sample") %>%
    mutate(n_coloniz = replace_na(n_coloniz, 0))
  
  return(coloniz_count)
}

# a function to compare the colonized feature in pups
identColoniz <- function(donor_feature,
                         cutoff_pres = 0,
                         mapfile = "../data/clean/mapfile.txt"){
  
  # Read mapfile
  map <- read.delim(mapfile, header = TRUE, 
                    stringsAsFactors = FALSE)
  
  coloniz_ident = donor_feature %>%
    filter(coverage >= cutoff_pres) %>%
    left_join(map %>% select(Sample, Group, Category)) %>%
    #filter(Category == "Pup") %>% 
    select(-Category)
  
  return(coloniz_ident)
}

# compare colonized feature in pups- wrapper for two above functions
CompFig <- function(feature_c = "../data/clean/binCover.txt", 
                        cutoff_pres = 0,
                        method = "Assembly",
                        mapfile = "../data/clean/mapfile.txt") {
  # DonorA
  donorA_feature <- DonFeature(feature_c = feature_c, donor = "DonA",
                               cutoff_pres = cutoff_pres, method = method)
  finalA_identity <- identColoniz(donor_feature = donorA_feature, 
                               cutoff_pres = cutoff_pres,
                               mapfile = mapfile) %>%
                    mutate(Donor = paste("DonorA"))
  # DonorB
  donorB_feature <- DonFeature(feature_c = feature_c, donor = "DonB",
                               cutoff_pres = cutoff_pres, method = method)
  finalB_identity <- identColoniz(donor_feature = donorB_feature, 
                                  cutoff_pres = cutoff_pres,
                                  mapfile = mapfile) %>%
                    mutate(Donor = paste("DonorB"))
  
  table <- finalA_identity %>%
    rbind(finalB_identity) 
  
  # compute feature importance
  feature_order <- table %>%
    group_by(Donor, Group, feature) %>%
    summarise(
      total_coverage = sum(coverage, na.rm = TRUE),
      detection_count = sum(coverage > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Donor, Group, desc(detection_count), desc(total_coverage)) %>%
    group_by(Donor, Group) %>%
    mutate(feature_rank = row_number())
  
  # join ranks back to original data
  table_ranked <- table %>%
    left_join(feature_order, by = c("Donor", "Group", "feature")) %>%
    mutate(feature = fct_reorder2(feature, Donor, -feature_rank))
  
  # Plot
  compfig <- ggplot(table_ranked, aes(Sample, feature, fill = coverage)) +
    geom_tile() +
    facet_grid(Donor ~ Group, space = "free", scales = "free") +
    scale_fill_gradient(low = "#eff3ff", high = "#08519c") +
    theme_classic() +
    theme(axis.text.x = element_blank())
  
  return(compfig)
}

#####################################################################

CompFigWithTree <- function(feature_c = "../data/clean/binCover.txt",
                            cutoff_pres = 0,
                            method = "Assembly",
                            mapfile = "../data/clean/mapfile.txt",
                            treeA_file = "../data/tree/DonorA_tree.nwk") {
  
  # Load tree
  treeA <- read.tree(treeA_file)
  
  # DonorA data
  donorA_feature <- DonFeature(feature_c = feature_c, donor = "DonA",
                               cutoff_pres = cutoff_pres, method = method)
  
  finalA_identity <- identColoniz(donor_feature = donorA_feature, 
                                  cutoff_pres = cutoff_pres,
                                  mapfile = mapfile) %>%
    mutate(Donor = "DonorA")
  
  # Filter to match tree tip labels
  tableA <- finalA_identity %>%
    filter(feature %in% treeA$tip.label)
  
  # Reshape for gheatmap
  heatmap_data <- tableA %>%
    select(feature, Sample, coverage) %>%
    spread(Sample, coverage, fill = 0) %>%
    column_to_rownames("feature")
  
  # Create annotation data for samples
  sample_group <- tableA %>%
    distinct(Sample, Group) %>%
    arrange(Sample) %>%
    column_to_rownames("Sample")
  
  
  # Plot tree + heatmap
  pA <- gheatmap(ggtree(treeA), heatmap_data,
                 offset = 0.4, width = 0.8,
                 font.size = 2,
                 #colnames_offset_y = 0.5,
                 colnames_angle = 90, colnames_position = "bottom") +
    geom_tiplab(size=2, align=TRUE, linesize=.5) + 
    #theme_tree2() +
    scale_fill_gradient(low = "white", high = "#08519c")
  
  return(pA)
}

####################################################################

CompFig_pheatmap <- function(table, flabel) {
  
  # table: comparison of features, use CompFeature() to generate this.
  # flabel: feature labels, how you want to be shown in the final fig
  # must contains feature, and label columns
 
  table %>%
    left_join(flabel) -> table
  
  # Reshape to wide format
  heatmap_data <- table %>%
    select(label, Sample, coverage) %>%
    pivot_wider(names_from = Sample, values_from = coverage, values_fill = 0) %>%
    column_to_rownames("label")
  
  # Create annotation data for samples
  sample_group <- table %>%
    distinct(Sample, Group) %>%
    column_to_rownames("Sample")
  label_group <- table %>%
    distinct(label, Donor) %>%
    column_to_rownames("label")
  
  annotation_colors_manual <- list(
    Group = c(DonorA = "#40004b", ParentA = "#9970ab", AB = "#e7d4e8",
              DonorB = "#00441b", ParentB = "#5aae61", BA = "#d9f0d3",
              ParentAB = "#e08214", ABAB = "#fdb863"),
    Donor = c(DonorA = "#7b3294", DonorB = "#008837")
  )
  
  # Plot heatmap with sample names and group annotation
  fig <- pheatmap(mat = heatmap_data,
                  #scale = "row",
                  cutree_cols=3,
                  annotation_col = sample_group,
                  annotation_row = label_group,
                  annotation_colors = annotation_colors_manual,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  show_colnames = FALSE,       # keeps sample names
                  show_rownames = TRUE,
                  color = colorRampPalette(c("white", "#0868ac"))(100))
                  #main = "Colonization Pattern")
  return(fig)
}

########################################################################
# count features. this is to generate a table of all counts. 
# this is a slight modification to CountFig function. only removing the 
# ggplot part of it.

CountFeature <- function(feature_c = "../data/clean/binCover.txt", 
                     cutoff_pres = 0,
                     method = "Assembly",
                     mapfile = "../data/clean/mapfile.txt") {
  
  # DonorA
  donorA_feature <- DonFeature(feature_c = feature_c, donor = "DonA",
                               cutoff_pres = cutoff_pres, method = method)
  finalA_count <- countColoniz(donor_feature = donorA_feature, 
                               cutoff_pres = cutoff_pres,
                               mapfile = mapfile) %>%
    mutate(Donor = paste("DonorA"))
  # DonorB
  donorB_feature <- DonFeature(feature_c = feature_c, donor = "DonB",
                               cutoff_pres = cutoff_pres, method = method)
  finalB_count <- countColoniz(donor_feature = donorB_feature, 
                               cutoff_pres = cutoff_pres,
                               mapfile = mapfile) %>%
    mutate(Donor = paste("DonorB"))
  
  
  # Combine and summarize data
  summary_df <- bind_rows(finalA_count, finalB_count)
 
  return(summary_df)
}


#########################################################################
# compare features. this is to generate a table of all pup features. 
# this is a slight modification to CompFig function. only removing the 
# ggplot part of it.

# compare colonized feature in pups- wrapper for two above functions
        CompFeature <- function(feature_c = "../data/clean/binCover.txt", 
                    cutoff_pres = 0,
                    method = "Assembly",
                    mapfile = "../data/clean/mapfile.txt") {
  # DonorA
  donorA_feature <- DonFeature(feature_c = feature_c, donor = "DonA",
                               cutoff_pres = cutoff_pres, method = method)
  finalA_identity <- identColoniz(donor_feature = donorA_feature, 
                                  cutoff_pres = cutoff_pres,
                                  mapfile = mapfile) %>%
    mutate(Donor = paste("DonorA"))
  # DonorB
  donorB_feature <- DonFeature(feature_c = feature_c, donor = "DonB",
                               cutoff_pres = cutoff_pres, method = method)
  finalB_identity <- identColoniz(donor_feature = donorB_feature, 
                                  cutoff_pres = cutoff_pres,
                                  mapfile = mapfile) %>%
    mutate(Donor = paste("DonorB"))
  
  table <- finalA_identity %>%
    rbind(finalB_identity) 
  
  return(table)
}

# a function to visualize colonization counts,
# based on a parsed table of counts
# it know where the input parsed file is located, 
# just select the feature of interest to visualize
CountFig <- function(ftype = "Species") {
  # select your feature type
  table <- read.csv("../data/clean/coloniz_visual.txt", sep = "\t")
  table <- table %>% filter(feature == ftype)
  
  sort_orders <- c("donor", "non-donor",
                   "parent", "non-parent", "mixed-parent",
                   "first", "second", "mixed")
  sort_groups <- c("DonorA", "DonorB",
                   "ParentA", "ParentB", "ParentAB",
                   "AB", "BA", "ABAB")
  table$order = factor(table$order, levels = sort_orders)
  table$Group = factor(table$Group, levels = sort_groups)
  # Combine and summarize data
  summary_df <- table %>%
    group_by(Donor, order, Category) %>%
    summarise(
      mean_coloniz = mean(n_coloniz, na.rm = TRUE),
      sd_coloniz = sd(n_coloniz, na.rm = TRUE),
      .groups = "drop"
    )
  
  p1 <- ggplot(summary_df, aes(x = order, y = mean_coloniz, 
                               group = Donor)) +
    geom_line(size = 1, aes(color = Donor)) +
    geom_point(size = 1, aes(color = Donor)) +
    geom_ribbon(aes(ymin = mean_coloniz - sd_coloniz, 
                    ymax = mean_coloniz + sd_coloniz,
                    fill = Donor), alpha = 0.2) +
    scale_color_manual(values = c("DonorA" = "#7b3294",
                                  "DonorB" = "#008837")) +
    scale_fill_manual(values = c("DonorA" = "#7b3294",
                                 "DonorB" = "#008837")) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.title = element_blank(),
      legend.position = c(0.98, 0.98),  # Top right corner
      legend.justification = c(1, 1)  # Align legend box to top right
    )
  
  summary_df2 <- table %>%
    group_by(Donor, Group, Category) %>%
    summarise(
      mean_coloniz = mean(n_coloniz, na.rm = TRUE),
      sd_coloniz = sd(n_coloniz, na.rm = TRUE),
      .groups = "drop"
    )
  
  p2 <- ggplot(summary_df2, aes(x = Group, y = mean_coloniz, 
                                fill = Donor)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("DonorA" = "#7b3294",
                                 "DonorB" = "#008837")) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 90),
          legend.position = "none")
  
  fig <- cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(3, 1))
  
  return(fig)
}

# a function to visualize Family-level profile of metaphlan4
profiler.fam <- function(metaphlan = "../data/clean/speciesAbund.txt",
                         mapfile = "../data/clean/mapfile.txt") {
  # requires long format of metaphlan4 feature table
  read.csv(mapfile, sep = "\t") %>%
    select(Sample, Group, Category) -> mapfile
  read.csv(metaphlan, sep = "\t") -> table
  table %>%  
    distinct(feature) %>%
    mutate(temp=feature) %>%
    separate(temp, c("Kingdom", "Phylum", "Class", 
                     "Order", "Family", "Genus", 
                     "Species", "SGB"), sep = "[|]") -> lineage
  
  table %>%
    left_join(mapfile) %>%
    left_join(lineage %>% 
                select(feature, Family)) -> Ctab
  ################################
  #select top 10 families
  Ctab %>%
    count(Family, wt = coverage) %>%
    filter(!is.na(Family)) %>%
    arrange(desc(n)) %>%
    pull(Family) -> tfam
  tfam10 <- tfam[1:10]
  # edit a new family column
  Ctab %>%
    mutate(Family_lab= case_when(Family %in% tfam10 ~ paste(Family),
                                 TRUE ~ paste("other"))) -> Ctab
  # sum the abundance for each class, across all IDs, & sort the result
  Ctab %>%
    count(Family_lab, wt = coverage) %>%
    arrange(n) %>%
    pull(Family_lab) -> sort.family
  
  group_order = c("DonorA","ParentA", "AB",
                  "DonorB", "ParentB", "BA",
                  "ParentAB", "ABAB")
  
  Ctab$Family_lab <- factor(Ctab$Family_lab, levels = sort.family)
  Ctab$Group <- factor(Ctab$Group, levels = group_order)
  
  family_col= c(
    "f__Bacteroidales_unclassified" = "#6a3d9a",
    "f__Muribaculaceae" = "#33a02c",
    "f__Lachnospiraceae" = "#1f78b4",
    "f__Rikenellaceae" = "#ff7f00",
    "f__Lactobacillaceae" = "#b2df8a",
    "f__Odoribacteraceae" = "#fdbf6f",
    "f__Helicobacteraceae" = "#ffff99",
    "f__FGB8157" =  "#a6cee3",
    "f__FGB9644" = "#e31a1c",
    "f__FGB9318" = "#f781bf",
    "other" = "#808080"
  )
  
  fig <- Ctab %>%
    ggplot(aes(Sample,coverage, fill=Family_lab, color=Family_lab)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = family_col,
                      guide = guide_legend(reverse = TRUE)) +
    scale_color_manual(values = family_col,
                       guide = guide_legend(reverse = TRUE)) +
    facet_grid(~Group, scales = "free", space = "free") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) + 
    theme_classic() +
    theme(legend.position = "bottom", legend.title = element_blank(),
          #legend.text = element_text(size=7),
          axis.text.x = element_blank(),
          strip.text.x = element_text(angle = 90)) +
    xlab("Samples") + ylab("Relative Abundance")
  
  return(fig)
  
}

vis.shannon <- function(Shannon = "../data/clean/diversity_shannon.txt",
                        Mapfile = "../data/clean/mapfile.txt") {
  
  read.csv(Mapfile, sep = "\t") -> mapfile
  read.csv(Shannon, sep = "\t") %>%
    left_join(mapfile) -> table
  
  group_order = c("DonorA","ParentA", "AB",
                  "DonorB", "ParentB", "BA",
                  "ParentAB", "ABAB")
  
  table$Group <- factor(table$Group, levels = group_order)
  ggplot(table, aes(x = Group, y = Shannon)) +
    geom_violin(outlier.color = NA) +
    geom_point(size=3,shape=21, 
               position = position_dodge(0.2)) + 
    ggpubr::geom_signif(test = "t.test",
                        comparisons = list(c("ParentA", "AB"),
                                           c("ParentB", "BA"),
                                           c("ParentAB", "ABAB")),
                        map_signif_level = FALSE,
                        step_increase = 0.1) +
    theme_classic() +
    theme(legend.title = element_blank(),
          text = element_text(size = 10),
          legend.position = "none", axis.title.x = element_blank()) +
    ylab("Shannon")
}

vis.pairBray <- function(Paired_bray = "../data/clean/paired_bray.txt",
                         Mapfile = "../data/clean/mapfile.txt"){
  
  table <- read.csv(Paired_bray, sep = "\t")
  
  group_order = c("DonorA","ParentA", "AB",
                  "DonorB", "ParentB", "BA",
                  "ParentAB", "ABAB")
  
  table$Group_2 <- factor(table$Group_2, levels = group_order)
  
  table %>%
    ggplot(aes(Group_1, value, color = Group_1)) +
    geom_boxplot(outlier.color = NA) +
    geom_point(size=3,shape=21, 
               position = position_dodge(0.2)) + 
    facet_grid(~Group_2, space = "free", scales = "free") +
    scale_color_manual(values = c("DonorA" = "#7b3294",
                                  "DonorB" = "#008837")) +
    ggpubr::geom_signif(test = "t.test",
                        comparisons = list(c("DonorA", "DonorB")),
                        map_signif_level = FALSE,
                        step_increase = 0.1, color = "black") +
    theme_classic() +
    theme(legend.title = element_blank(),
          text = element_text(size = 10),
          legend.position = "none", axis.title.x = element_blank()) +
    ylab("Bray-Curtis")
  
  
}

vis.contigs <- function(
    cumulative_df = "../data/clean/assembly_quality.txt") {
  
  read.csv(cumulative_df, sep = "\t") %>%
  mutate(Method = factor(Method, levels = c("Single", "Co-assembly"))) %>%
    ggplot(aes(Contigid2, Cumulative, color=Sample)) +
    geom_line(size=2) +
    facet_grid(~Method+Sample, space = "free", scales = "free") +
    theme_bw() +
    scale_color_manual(values = c("DonorA" = "#7b3294",
                                  "DonorB" = "#008837")) +
    ggtitle("Cumulative length of assembled contigs over 1kb") +
    ylab("Cumulative length in Mbp") +
    xlab("Contigs are ordered from largest (contig #1) to smallest (x1000).") +
    theme(text = element_text(size= 10),
          axis.text.x = element_text(angle = 0),
          legend.position = "none", legend.title = element_blank())
}


