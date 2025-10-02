sum_bins = function(df){
  breaks = n_distinct(df$bins)
  labs = paste('bin', 1:breaks, sep = '')
  df = (df
        %>% group_by(bins)
        %>% summarize(mean = mean(raw_count),
                      sd = sd(raw_count),
                      n = length(raw_count))
        %>% ungroup()
        %>% mutate(bins = factor(bins),
                   bins = order_levs(bins, mean),
                   bins = factor(bins, 
                                 labels = paste('bin', 1:breaks, sep = ''))))
  return(df)
}

mk_brk_df = function(df, breaks){
  bin_df = (df
            %>% select(raw_count)
            %>% mutate(bins = cut(raw_count, breaks = breaks)))
  bin_df = sum_bins(bin_df)
  return(bin_df)
}

plt_mean_sd = function(df, loglog = 'n'){
  plt = ggplot(df, aes(mean, sd)) +
    geom_point(alpha = 0.4) +
    xlim(0, 5e4)
  if (loglog == 'y'){
    plt = plt + 
      scale_y_log10() +
      scale_x_log10(limits = c(1, 5e4))
  }
  return(plt)
}

plt_bs = function(df){
  plt = ggplot(df, aes(bins, n)) +  
    geom_point() +  rotate_ticks() +
    ylab('Samples per bin') +
    xlab('bins (sorted by mean)') +
    ylim(0, 200)
  
  return(plt)
}

mk_qq = function(df, grp){
  cols = paste(c('fit', 'res'), grp, sep = '_')
  plt = ggplot(df, aes(sample = .data[[cols[2]]])) +
    geom_qq() +
    geom_qq_line()
  return(plt)
}

mk_fr = function(df, grp, abs = FALSE, guide = 'legend'){
  cols = paste(c('fit', 'res'), grp, sep = '_')
  if (abs){
    plt = ggplot(df, aes(x = .data[[cols[1]]],
                         y = abs(.data[[cols[2]]])))
  } else {
    plt = ggplot(df, aes(x = .data[[cols[1]]],
                         y = .data[[cols[2]]]))
    
  }
  
  plt = plt + geom_point(aes(colour = feature)) +
    scale_colour_brewer(palette = 'Dark2', guide = guide) +
    geom_smooth(method = 'lm')
  return(plt)
}

mk_plt_row = function(df, grp){
  rw = plot_grid(mk_qq(fit_df, grp), mk_fr(fit_df, grp, guide = 'none'),
                 mk_fr(fit_df, grp,TRUE, guide = 'none'),
                 ncol = 3, 
                 labels = grp,
                 label_size = 10)
  return(rw)
}