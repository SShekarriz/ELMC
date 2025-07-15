
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
    summarise(n_coloniz= n()) %>%
    left_join(map %>% select(Sample, Group, Category))
  
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
  
  countfig <- finalA_count%>%
    rbind(finalB_count) %>%
    ggplot(
      aes(Group, n_coloniz)) +
    geom_boxplot() + geom_point() +
    facet_grid(Donor~Category, space = "free", scales = "free_x")
  
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
  
  compfig <- finalA_identity%>%
    rbind(finalB_identity) %>%
    ggplot(
      aes(feature, Sample, fill= coverage)) +
    geom_tile() +
    facet_grid(Group~Donor, space = "free", scales = "free")
  
  return(compfig)
}

