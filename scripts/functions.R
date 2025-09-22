
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

# get total coloniation count- wraper for two above functions
CountFig <- function(feature_c = "../data/clean/binCover.txt", 
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
  summary_df <- bind_rows(finalA_count, finalB_count) %>%
    group_by(Donor, Category, Group) %>%
    summarise(
      mean_coloniz = mean(n_coloniz, na.rm = TRUE),
      sd_coloniz = sd(n_coloniz, na.rm = TRUE),
      .groups = "drop"
    )
  
  summary_df$Group <- factor(summary_df$Group, 
                             levels = c("DonorA", "DonorB",
                               "ParentA", "ParentB", "ParentAB",
                               "AB", "BA", "ABAB"))

  countfig <- ggplot(summary_df, aes(x = Group, 
                                     y = mean_coloniz, 
                                     group = Donor)) +
    # Add shaded rectangle between DonorA and DonorB
    geom_rect(data = NULL, 
              aes(xmin = which(levels(summary_df$Group) == "DonorA"), 
                  xmax = which(levels(summary_df$Group) == "ParentAB"), 
                  ymin = -Inf, ymax = Inf),
              fill = "#252525", alpha = 0.009) +
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
      axis.text.x = element_text(angle = 90),
      legend.title = element_blank(),
      legend.position = c(0.98, 0.98),  # Top right corner
      legend.justification = c(1, 1)  # Align legend box to top right
    )
  
  
  return(countfig)
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
    filter(Category == "Pup") %>% select(-Category)
  
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

CompFigWithTree_pheatmap <- function(feature_c = "../data/clean/binCover.txt",
                                     cutoff_pres = 0,
                                     method = "Assembly",
                                     mapfile = "../data/clean/mapfile.txt") {
  
  # DonorA data
  donorA_feature <- DonFeature(feature_c = feature_c, donor = "DonA",
                               cutoff_pres = cutoff_pres, method = method)
  
  finalA_identity <- identColoniz(donor_feature = donorA_feature, 
                                  cutoff_pres = cutoff_pres,
                                  mapfile = mapfile) %>%
    mutate(Donor = "DonorA")
  
  # Reshape to wide format
  heatmap_data <- finalA_identity %>%
    select(feature, Sample, coverage) %>%
    pivot_wider(names_from = Sample, values_from = coverage, values_fill = 0) %>%
    column_to_rownames("feature")
  
  # Create annotation data for samples
  sample_group <- finalA_identity %>%
    distinct(Sample, Group) %>%
    column_to_rownames("Sample")
  
  # Plot heatmap with sample names and group annotation
  pheatmap(mat = heatmap_data,
           annotation_col = sample_group,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_colnames = TRUE,       # keeps sample names
           show_rownames = TRUE,
           color = colorRampPalette(c("#eff3ff", "#08519c"))(100),
           main = "DonorA Colonization Pattern")
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


