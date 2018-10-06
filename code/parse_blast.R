# load libraries
library(readr)
library(dplyr)

# read in file
data <- readr::read_tsv("Desktop/firstset.txt",
                        comment = "#", 
                        col_names=c("query","subject","evalue","alnlength",
                                    "identity","identical","mismatches",
                                    "qframe","sframe")) %>%
  dplyr::filter(!(query == subject)) %>% # query does not equal subject
  dplyr::filter(!(qframe == sframe)) %>% # frames are not equal
  dplyr::filter(identity >= 80.0) %>% # minimum percentage identity
  dplyr::filter(evalue <= 1e-10) # minimum E-value

# get maximum alignment length
aln_lookup_table <- data %>%
  dplyr::group_by(query, subject) %>% 
  dplyr::summarise(max_aln_length = max(alnlength))

# filter by maximum alignment length
data <- dplyr::inner_join(data, aln_lookup_table, 
                          by = c("query", "subject","alnlength" = "max_aln_length"))
  
# write standard file
data %>%
  dplyr::filter(alnlength >= 50) %>% # minimum sequence alignment length
  write_tsv(path="Desktop/standardset.tsv")

# write conservative file
data %>%
  dplyr::filter(alnlength >= 100) %>% # minimum sequence alignment length
  write_tsv(path="Desktop/conservativeset.tsv")
  
# print(sessionInfo())
# R version 3.4.2 (2017-09-28)
# Platform: x86_64-apple-darwin17.0.0 (64-bit)
# Running under: macOS High Sierra 10.13.1
# 
# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] bindrcpp_0.2 dplyr_0.7.4  readr_1.1.1 
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_0.12.13     lattice_0.20-35  assertthat_0.2.0 plyr_1.8.4       grid_3.4.2       R6_2.2.2         nlme_3.1-131     magrittr_1.5    
# [9] rlang_0.1.4      tools_3.4.2      glue_1.2.0       hms_0.3          yaml_2.1.14      compiler_3.4.2   pkgconfig_2.0.1  bindr_0.1       
# [17] tibble_1.3.4   