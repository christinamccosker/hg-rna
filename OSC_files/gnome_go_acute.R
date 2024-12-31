################################################################################
#
## Load required libraries
#

library(rbioapi)

################################################################################
#
## Import data
#

peak <- read.csv("Output Files/gnome_female_edgeR_null_DEG_peak.csv")

peak1 <- peak[,1:250]
peak2 <- peak[,251:500]
peak3 <- peak[,501:750]
peak4 <- peak[,751:1000] 

postfilter <- read.csv("Output Files/gnome_female_edgeR_null_DEG_postfilter.csv")

################################################################################
#
## Set up results lists
#


go_cc_peak1 <- setNames(vector('list', 250), 1:250)
go_cc_peak2 <- setNames(vector('list', 250), 1:250)
go_cc_peak3 <- setNames(vector('list', 250), 1:250)
go_cc_peak4 <- setNames(vector('list', 250), 1:250)

################################################################################
#
## peak Loops (loop over peak1-4)
#




### Start here ###


# Cellular Component
for (i in 1:250) {
  genes <- peak4[[i]]
  background <- postfilter[[i]]
  cc_test <-
    rba_panther_enrich(genes = genes, organism = 9606, annot_dataset = "GO:0005575",
                       test_type = "FISHER", correction = "FDR", ref_genes = background,
                       ref_organism = 9606, cutoff = 0.05)
  go_cc_peak4[[i]] <- cc_test$result$term.label
}






# Save lists!
cc_peak_all <- 
  c(go_cc_peak1, go_cc_peak2, go_cc_peak3,go_cc_peak4)




max_length <- max(unlist(lapply(cc_peak_all, length)))
cc_peak_all <-
  lapply(cc_peak_all, function(x) {ans <- rep(NA,length=max_length);
  ans[0:length(x)]<- x;
  return(ans)})
cc_peak_all <- do.call(cbind, cc_peak_all)
write.csv(cc_peak_all, file = "Output Files/gnome_female_edgeR_null_GO_cc_peak.csv", row.names = FALSE)
