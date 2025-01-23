#
#
#
#
#
#
#
#
#
#
#
#
#

#-
library(tidyverse)
library(tidyr)
# load library to read excel files
library(readxl)
library(ggvenn)
library(knitr)
library(glue)
setwd("/Users/lmz306/Library/CloudStorage/OneDrive-UniversityofCopenhagen/delta_ccr5/2024_07_15_latest_results_8del/")
biorxiv  <- read_excel("ccr5_biorxiv_suppl.xlsx", sheet = "Table S4 BP_Ancient samples ") %>%
arrange(Sample)



new  <- read_excel("rs333_kirstine_18_07_24.xlsx")
setwd("/Users/lmz306/Library/CloudStorage/OneDrive-UniversityofCopenhagen/delta_ccr5/2024_07_15_latest_results_8del/2024_07_18_Newfiles_Kirstine/PyScriptDelres_no_duplicated_RISE")
# View(biorxiv)
# View(new)

deletions  <- list.files(pattern = "_E"); deletions
deletions  <- str_remove(deletions, "_E.xlsx")

deletion  <- deletions[7]; deletion
df  <- read_excel(glue("Py{deletion}/results_{deletion}.xlsx"))

View(df %>% 
filter(Permissive_filter != No_filter) %>% 
arrange(No_filter) %>% 
select(Sample:Permissive_filter, N_reads_ref:N_reads_del, sum_x, referencecount, purereferencecount, haplocount, purehaplocount, starts_with("SNP")))


View(df %>% 
# filter(Strict_filter != Permissive_filter) %>%
# filter(Strict_filter != No_filter) %>%
arrange(No_filter) %>% 
select(Sample:Strict_filter, N_reads_ref:N_reads_del, sum_x, referencecount, purereferencecount, haplocount, purehaplocount, starts_with("SNP")))



# deletion 7 RISE98 weird change from RR to RD. Do not agree but I'll keep it like that 

#
#
#
