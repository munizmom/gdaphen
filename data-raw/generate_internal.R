i















## code to prepare `generate_internal` dataset goes here
phenoInput<- as.data.frame(read_excel(paste0(f_inputFile,"phenoParameters_Scn9a.xlsx"),sheet =1, col_names =T)); 

usethis::use_data(generate_internal, overwrite = TRUE)
