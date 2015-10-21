
# Run the run_kallisto.R script in the /data/kallisto/ directory

library(sleuth)
library(edgeR)

kdirs <- dir("../data/kallisto","^[NT]",full=TRUE)
names(kdirs) <- basename(kdirs)

df <- data.frame(sample=kdirs, group=substr(names(kdirs),1,1),
                 stringsAsFactors=FALSE)

so <- sleuth_prep(kdirs, df, ~group)

# function to extract from sleuth/kallisto object
extract_sleuth <- function(slo, col="est_counts") {
  x <- sapply(slo$kal, function(u) u$abundance[[col]])
  rownames(x) <- slo$kal[[1]]$abundance$target_id
  x
}


countsT <- extract_sleuth(so,"est_counts")
tpms <- extract_sleuth(so,"tpm")

library(rtracklayer)
anno <- import("../data/anno/transcripts.gtf")

tid <- rownames(countsT)
m <- match(tid, anno$transcript_id)
gid <- anno$gene_id[m]

inds <- split(1:nrow(tpms),gid)

# calculate transcript proportions by gene
propsT <- lapply(inds, function(u) { 
  z <- tpms[u,,drop=FALSE]; 
  t(t(z)/colSums(z)) 
})
propsT <- do.call("rbind",propsT)

# find the max proportion over all samples
library(matrixStats)
rm <- rowMaxs(propsT,na.rm=TRUE)

# filter the annotation, export it
annos <- anno[ anno$transcript_id %in% rownames(propsT)[rm>.05] ]
export(annos, "../data/anno/transcripts_kallisto_filtered.gtf")


## -------------------------
## goto dexseq_count.R
## -------------------------

