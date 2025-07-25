
as_tibble_rnm <- function(dat, rowname_col = "rowname") {
  
  #Convert matrix or data.frame
  
    
    # Convert to tibble
    tibble_df <- as_tibble(dat) %>% 
      # Add the row names as a new column
      mutate(!!rowname_col := rownames(dat)) %>% 
      relocate(!!sym(rowname_col), .before = 1)
    
  
  return(tibble_df)
}

