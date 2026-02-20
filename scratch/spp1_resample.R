#load packages
library(limma)
library(Biobase)

## 1. Read and format data
#read expression data
load("C:/Users/matts/Desktop/projects/spp1/SWOG_GEX_ProteinCoding_without_duplica.Rdata", )

#convert to data frame
swog_expressions = as.data.frame(SWOG_matrix)

#remove original matrix
rm(SWOG_matrix)

#read in metadata
swog_meta = read.table("C:/Users/matts/Desktop/projects/spp1/S1314 metadata.txt", sep = "\t", header = TRUE)

## 2.
#subset expression matrix, on subtypes
#create helper function
expression_subs = function(subtypes = NULL,
                           expression_data = swog_expressions,
                           this_metadata = swog_meta,
                           calc_degs = FALSE,
                           exclude = NULL,
                           plot = FALSE,
                           out_path = NULL,
                           write_degs = FALSE,
                           this_variable = "LundTax2023_simple",
                           sample_id_col = "ALTPATID",
                           return_this = "degs"){

  if(write_degs && is.null(out_path)){
    stop("No output path provided...")
  }
  
  if(!is.null(subtypes)){
    #subset metadata to sample IDs based on the subtype
    message("Fiiltering samples based on selected subtype...")
    my_metadata = dplyr::filter(this_metadata, !!as.symbol(this_variable) %in% subtypes)
  }else{
    message("No subtype is defined, the function will use all sample IDs available in the metadata...")
    my_metadata = this_metadata
  }


  #return the sample IDs for the selected subtype
  my_samples = my_metadata %>%
    pull(sample_id_col)
  
  if(!is.null(exclude)){
    my_samples = my_samples[!my_samples %in% exclude]
    message(paste0("The following samples are excluded: ", exclude))
  }

  #return the expression matric for the selected subtypes
  my_expr = dplyr::select(expression_data, my_samples)

  #generate an ExpressionSet (Biobase) for the selected subtypes
  my_eset = ExpressionSet(assayData = as.matrix(my_expr))

  #get the design matrix
  my_pCR = factor(my_metadata$pT0)
  my_design = model.matrix(~my_pCR)

  message("Calculate DEGs (moderated t-test) for pCR overall")

  fit = lmFit(my_eset,
              my_design)

  fit = eBayes(fit,
               trend = TRUE,
               robust = TRUE)

  results = decideTests(fit)

  summary(results)

  my_degs = topTable(fit,
                     coef = "my_pCR1",
                     n = nrow(my_expr))
  if(write_degs){
    write.table(my_degs,
                file = paste0(out_path, "DEGs.txt"),
                sep = "\t",
                col.names = TRUE)
  }

  if(plot){
    plotMD(fit,
             coef = "my_pCR",
             status = results[,5],
             values = c(1,-1),
             hl.col = c("red","blue"))
  }

  #deal with returns
  if(return_this == "sample_ids"){
    message(paste0(length(my_samples), " sample IDs returned for ", subtypes))
    return(my_samples)
  }else if (return_this == "metadata"){
    message(paste0("Metadata returned for the following subtypes;  ", subtypes))
    return(my_metadata)
  }else if(return_this == "expressions"){
    message(paste0("Expression values returned for ", subtypes))
    return(my_expr)
  }else if(return_this == "eset"){
    message(paste0("ExpressionSet returned for ", subtypes))
    return(my_eset)
  }else if(return_this == "design_matrix"){
    message(paste0("Design matrix returned for ", subtypes))
    return(my_design)
  }else if(return_this == "degs"){
    message(paste0("DEGs returned for ", subtypes))
    return(my_degs)
  }else if(return_this == "all"){
    message(paste0("All data returned for ", subtypes))
    my_list = list(sample_ids = my_samples,
                   metadata = my_metadata,
                   expression_data = my_expr,
                   eset = my_eset,
                   design = my_design,
                   degs = my_degs)
  }else if(return_this == "nothing"){
    message("Nothing is returned")
    return()
  }else{
    stop("Possible values for return_this are; sample_ids, metadata, expressions, eset, design_matrix, degs, all, and nothing...")
  }
}

#return expression values for each subtype
test_df = expression_subs(expression_data = swog_expressions, 
                          this_metadata = swog_meta, 
                          subtypes = "Uro",
                          return_this = "degs")



