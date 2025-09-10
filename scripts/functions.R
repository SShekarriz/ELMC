
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
    #group_by(Sample) %>%
    #summarise(n_coloniz= n()) %>%
    #left_join(map %>% select(Sample, Group, Category))
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
  
  # Plot
  #countfig <- ggplot(summary_df, aes(x = Group, 
  #                                   y = mean_coloniz, 
  #                                   group = Donor)) +
  #  geom_line(size =1, aes(color = Donor)) +
  #  geom_point(size = 1, aes(color = Donor)) +
  #  geom_ribbon(aes(ymin = mean_coloniz - sd_coloniz, 
  #                    ymax = mean_coloniz + sd_coloniz,
  #                  fill = Donor), alpha = 0.2) +
  #  theme_bw() +
  #  theme(
  #        axis.title = element_blank(),
  #        axis.text.x = element_text(angle = 90))
  
  
  countfig <- ggplot(summary_df, aes(x = Group, 
                                     y = mean_coloniz, 
                                     group = Donor)) +
    # Add shaded rectangle between DonorA and DonorB
    geom_rect(data = NULL, 
              aes(xmin = which(levels(summary_df$Group) == "DonorA"), 
                  xmax = which(levels(summary_df$Group) == "ParentAB"), 
                  ymin = -Inf, ymax = Inf),
              fill = "grey", alpha = 0.05) +
    geom_line(size = 1, aes(color = Donor)) +
    geom_point(size = 1, aes(color = Donor)) +
    geom_ribbon(aes(ymin = mean_coloniz - sd_coloniz, 
                    ymax = mean_coloniz + sd_coloniz,
                    fill = Donor), alpha = 0.2) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 90),
      legend.position = c(1, 1),  # Top right corner
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
  compfig <- ggplot(table_ranked, aes(feature, Sample, fill = coverage)) +
    geom_tile() +
    facet_grid(Group ~ Donor, space = "free", scales = "free") +
    theme(axis.text.x = element_blank())
  
  return(compfig)
}


# compare colonized feature in pups- by adding phylogeny of features
CompPhyloFig <- function(feature_c = "../data/clean/binCover.txt", 
                    cutoff_pres = 0,
                    method = "Assembly",
                    mapfile = "../data/clean/mapfile.txt",
                    DonA_tree= "../data/clean/MAGs_DonA.tree") {
  
  # DonorA
  donorA_feature <- DonFeature(feature_c = feature_c, donor = "DonA",
                               cutoff_pres = cutoff_pres, method = method)
  finalA_identity <- identColoniz(donor_feature = donorA_feature, 
                                  cutoff_pres = cutoff_pres,
                                  mapfile = mapfile) %>%
    mutate(Donor = paste("DonorA"))

  table <- finalA_identity
  
  # read tree file 
  Atree <- read.tree(DonA_tree)
  Ap_tree <- ggtree(Atree) + geom_tiplab(align = TRUE)
  
  #order the feature by tip.label
  data_tile <- data.frame(feature = Atree$tip.label)
  table$feature <- factor(data_tile$feature, 
                   levels = Atree$tip.label[order(match(Atree$tip.label,
                   Atree$tip.label))])
  
  # Plot
  compfig <- ggplot(table, aes(feature, Sample, fill = coverage)) +
    geom_tile() +
    facet_grid(Group ~ Donor, space = "free", scales = "free") +
    scale_fill_viridis_c() +
    theme(axis.text.x = element_blank())
  
  return(compfig)
}








