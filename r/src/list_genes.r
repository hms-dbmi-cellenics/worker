#
#
# req has body with        {'selectFields'= {"gene_names", "dispersions"},
#                          'orderBy'= {"gene_names", "dispersions"},
#                          'geneNamesFilter' = string with r placeholders (^ before for endswith and $ after for startswith),
#                          'orderDirection'= {DESC, ASC} 
#                          'offset'= int: number of genes from the first one in the sorted and filtered dataframe to offset the results by
#                          'limit'= int: how many results per page.
#
# returns: df with variance.standarized, SYMBOL (gene_names) and full_count collumn with number of genes after filtering. 
# 
# We decided to return variance.standarized because it's a commonly used indicator of dispersion that accounts for the mean.[1] 
# 
# [1]Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W. M., ... & Satija, R. (2019). Comprehensive integration of single-cell data. Cell, 177(7), 1888-1902.
#
getList <- function(req){
  message(data)
  selectFields <- req$body$selectFields
  orderBy <- req$body$orderBy
  orderDirection <- req$body$orderDirection == "DESC"
  #Gene dispersion slot generated in data ingest script with the same info as the meta.features slot but with the annotated genes
  df <- data@misc$gene_dispersion
  #Saving the full count in case there's no filter set.
  df$full_count = nrow(df)
  if ("geneNamesFilter" %in% names(req$body)){
    #create logical vector with grepl. The placeholders for searching are set on the UI side.
    df <- df[grepl(req$body$geneNamesFilter, df$SYMBOL,ignore.case = TRUE),]
    #Update the number of rows for the UI to know how many pages to create.
    if (nrow(df)>0){
      df$full_count = nrow(df)
    }    
  }  
  #Sorting the dataset
  if (orderBy == "dispersions"){
    df <- df[order(df$variance.standardized, decreasing = orderDirection),]
  }else if (orderBy == "gene_names"){
    df <- df[order(df$SYMBOL, decreasing = orderDirection),]
  }
  #
  #  Offset 0 ->  offset + 1 = 1 ( r arrays start from 1)
  #  Limit 4 -> limit + offset -1 = 4
  #
  #  [ 0 ! 1 ! 2 ! 3 ! 4 ! 5 ! 6 ! 7 ! 8 ! 9 ]
  #    !           !       
  #    O           L
  #
  #  Offset 4 ->  offset + 1 = 5 
  #  Limit 4 -> limit + offset - 1 = 8
  #
  #  [ 0 ! 1 ! 2 ! 3 ! 4 ! 5 ! 6 ! 7 ! 8 ! 9 ]
  #                    !           !
  #                    O           L
  #
  offset <- req$body$offset + 1 
  limit <- req$body$limit - 1
  #Remove NA's and offset the df
  df <- na.omit(df[(offset):(offset+limit),])
  #clean out unwanted columns
  res <- data.frame(gene_names = df$SYMBOL, dispersions = df$variance.standardized, full_count = df$full_count)
  return(res)
}