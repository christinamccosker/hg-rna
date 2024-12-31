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

acute <- read.csv("gnome_female_edgeR_DEG_null_acute.csv")
peak <- read.csv("gnome_female_edgeR_DEG_null_peak.csv")
late <- read.csv("gnome_female_edgeR_DEG_null_late.csv")

postfilter <- read.csv("gnome_female_edgeR_DEG_null_postfilter.csv")

################################################################################
#
## Set up results lists
#

go_bp_acute <- setNames(vector('list', 1000), 1:1000)
go_bp_peak <- setNames(vector('list', 1000), 1:1000)
go_bp_late <- setNames(vector('list', 1000), 1:1000)

go_mf_acute <- setNames(vector('list', 1000), 1:1000)
go_mf_peak <- setNames(vector('list', 1000), 1:1000)
go_mf_late <- setNames(vector('list', 1000), 1:1000)

go_cc_acute <- setNames(vector('list', 1000), 1:1000)
go_cc_peak <- setNames(vector('list', 1000), 1:1000)
go_cc_late <- setNames(vector('list', 1000), 1:1000)

################################################################################
#
## Acute Loop
#

for (i in 1:1000) {
  genes <- acute[[i]]
  background <- postfilter[[i]]
  bp_test <-
    rba_panther_enrich(genes = genes, organism = 9606, annot_dataset = "GO:0008150",
                       test_type = "FISHER", correction = "FDR", ref_genes = background,
                       ref_organism = 9606, cutoff = 0.05)
  go_bp_acute[[i]] <- bp_test$result$term.label

  mf_test <-
    rba_panther_enrich(genes = genes, organism = 9606, annot_dataset = "GO:0003674",
                       test_type = "FISHER", correction = "FDR", ref_genes = background,
                       ref_organism = 9606, cutoff = 0.05)
  go_mf_acute[[i]] <- mf_test$result$term.label

  cc_test <-
    rba_panther_enrich(genes = genes, organism = 9606, annot_dataset = "GO:0005575",
                       test_type = "FISHER", correction = "FDR", ref_genes = background,
                       ref_organism = 9606, cutoff = 0.05)
  go_cc_acute[[i]] <- cc_test$result$term.label
}

################################################################################
#
## Save acute files
#

max_length <- max(unlist(lapply(go_bp_acute, length)))
go_bp_acute_filled <-
  lapply(go_bp_acute, function(x) {ans <- rep(NA,length=max_length);
  ans[0:length(x)]<- x;
  return(ans)})
go_bp_acute_final <- do.call(cbind, go_bp_acute_filled)
write.csv(go_bp_acute_final, file = "gnome_female_edgeR_GO_bp_acute.csv", row.names = FALSE)

max_length <- max(unlist(lapply(go_mf_acute, length)))
go_mf_acute_filled <-
  lapply(go_mf_acute, function(x) {ans <- rep(NA,length=max_length);
  ans[0:length(x)]<- x;
  return(ans)})
go_mf_acute_final <- do.call(cbind, go_mf_acute_filled)
write.csv(go_mf_acute_final, file = "gnome_female_edgeR_GO_mf_acute.csv", row.names = FALSE)

max_length <- max(unlist(lapply(go_cc_acute, length)))
go_cc_acute_filled <-
  lapply(go_cc_acute, function(x) {ans <- rep(NA,length=max_length);
  ans[0:length(x)]<- x;
  return(ans)})
go_cc_acute_final <- do.call(cbind, go_cc_acute_filled)
write.csv(go_cc_acute_final, file = "gnome_female_edgeR_GO_cc_acute.csv", row.names = FALSE)

################################################################################
#
## Peak Loop
#

for (i in 1:1000) {
  genes <- peak[[i]]
  background <- postfilter[[i]]
  bp_test <-
    rba_panther_enrich(genes = genes, organism = 9606, annot_dataset = "GO:0008150",
                       test_type = "FISHER", correction = "FDR", ref_genes = background,
                       ref_organism = 9606, cutoff = 0.05)
  go_bp_peak[[i]] <- bp_test$result$term.label

  mf_test <-
    rba_panther_enrich(genes = genes, organism = 9606, annot_dataset = "GO:0003674",
                       test_type = "FISHER", correction = "FDR", ref_genes = background,
                       ref_organism = 9606, cutoff = 0.05)
  go_mf_peak[[i]] <- mf_test$result$term.label

  cc_test <-
    rba_panther_enrich(genes = genes, organism = 9606, annot_dataset = "GO:0005575",
                       test_type = "FISHER", correction = "FDR", ref_genes = background,
                       ref_organism = 9606, cutoff = 0.05)
  go_cc_peak[[i]] <- cc_test$result$term.label
}

################################################################################
#
## Save peak files
#

max_length <- max(unlist(lapply(go_bp_peak, length)))
go_bp_peak_filled <-
  lapply(go_bp_peak, function(x) {ans <- rep(NA,length=max_length);
  ans[0:length(x)]<- x;
  return(ans)})
go_bp_peak_final <- do.call(cbind, go_bp_peak_filled)
write.csv(go_bp_peak_final, file = "gnome_female_edgeR_GO_bp_peak.csv", row.names = FALSE)

max_length <- max(unlist(lapply(go_mf_peak, length)))
go_mf_peak_filled <-
  lapply(go_mf_peak, function(x) {ans <- rep(NA,length=max_length);
  ans[0:length(x)]<- x;
  return(ans)})
go_mf_peak_final <- do.call(cbind, go_mf_peak_filled)
write.csv(go_mf_peak_final, file = "gnome_female_edgeR_GO_mf_peak.csv", row.names = FALSE)

max_length <- max(unlist(lapply(go_cc_peak, length)))
go_cc_peak_filled <-
  lapply(go_cc_peak, function(x) {ans <- rep(NA,length=max_length);
  ans[0:length(x)]<- x;
  return(ans)})
go_cc_peak_final <- do.call(cbind, go_cc_peak_filled)
write.csv(go_cc_peak_final, file = "gnome_female_edgeR_GO_cc_peak.csv", row.names = FALSE)

################################################################################
#
## Late Loop
#

for (i in 1:1000) {
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
write.csv(go_bp_late_final, file = "gnome_female_edgeR_GO_bp_late.csv", row.names = FALSE)

max_length <- max(unlist(lapply(go_mf_late, length)))
go_mf_late_filled <-
  lapply(go_mf_late, function(x) {ans <- rep(NA,length=max_length);
  ans[0:length(x)]<- x;
  return(ans)})
go_mf_late_final <- do.call(cbind, go_mf_late_filled)
write.csv(go_mf_late_final, file = "gnome_female_edgeR_GO_mf_late.csv", row.names = FALSE)

max_length <- max(unlist(lapply(go_cc_late, length)))
go_cc_late_filled <-
  lapply(go_cc_late, function(x) {ans <- rep(NA,length=max_length);
  ans[0:length(x)]<- x;
  return(ans)})
go_cc_late_final <- do.call(cbind, go_cc_late_filled)
write.csv(go_cc_late_final, file = "gnome_female_edgeR_GO_cc_late.csv", row.names = FALSE)
