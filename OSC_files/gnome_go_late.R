################################################################################
#
## Load required libraries
#

library(rbioapi)
library(tidyverse)

################################################################################
#
## Import data
#

late <- read.csv("Output Files/gnome_female_edgeR_null_DEG_late.csv")

postfilter <- read.csv("Output Files/gnome_female_edgeR_null_DEG_postfilter.csv")

################################################################################
#
## Set up results lists
#

go_bp_late <- setNames(vector('list', 1000), 1:1000)

go_mf_late <- setNames(vector('list', 1000), 1:1000)

go_cc_late <- setNames(vector('list', 1000), 1:1000)

################################################################################
#
## Late Loop
#

for (i in 782:1000) {
  genes <- late[[i]]
  background <- postfilter[[i]]
  bp_test <-
    rba_panther_enrich(genes = genes, organism = 9606, annot_dataset = "GO:0008150",
                       test_type = "FISHER", correction = "FDR", ref_genes = background,
                       ref_organism = 9606, cutoff = 0.05)
  go_bp_late[[i]] <- bp_test$result$term.label

  mf_test <-
    rba_panther_enrich(genes = genes, organism = 9606, annot_dataset = "GO:0003674",
                       test_type = "FISHER", correction = "FDR", ref_genes = background,
                       ref_organism = 9606, cutoff = 0.05)
  go_mf_late[[i]] <- mf_test$result$term.label

  cc_test <-
    rba_panther_enrich(genes = genes, organism = 9606, annot_dataset = "GO:0005575",
                       test_type = "FISHER", correction = "FDR", ref_genes = background,
                       ref_organism = 9606, cutoff = 0.05)
  go_cc_late[[i]] <- cc_test$result$term.label
}

################################################################################
#
## Save late files
#

max_length <- max(unlist(lapply(go_bp_late, length)))
go_bp_late_filled <-
  lapply(go_bp_late, function(x) {ans <- rep(NA,length=max_length);
  ans[0:length(x)]<- x;
  return(ans)})
go_bp_late_final <- do.call(cbind, go_bp_late_filled)
write.csv(go_bp_late_final, file = "gnome_female_edgeR_null_GO_bp_late.csv", row.names = FALSE)

max_length <- max(unlist(lapply(go_mf_late, length)))
go_mf_late_filled <-
  lapply(go_mf_late, function(x) {ans <- rep(NA,length=max_length);
  ans[0:length(x)]<- x;
  return(ans)})
go_mf_late_final <- do.call(cbind, go_mf_late_filled)
write.csv(go_mf_late_final, file = "gnome_female_edgeR_null_GO_mf_late.csv", row.names = FALSE)

max_length <- max(unlist(lapply(go_cc_late, length)))
go_cc_late_filled <-
  lapply(go_cc_late, function(x) {ans <- rep(NA,length=max_length);
  ans[0:length(x)]<- x;
  return(ans)})
go_cc_late_final <- do.call(cbind, go_cc_late_filled)
write.csv(go_cc_late_final, file = "gnome_female_edgeR_null_GO_cc_late.csv", row.names = FALSE)
