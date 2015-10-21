
# To get kallisto "expected counts", run the run_kallisto.R 
# script in the /data/kallisto/ directory

library(sleuth)
library(edgeR)

kdirs <- dir("../data/kallisto","^[NT]",full=TRUE)
names(kdirs) <- basename(kdirs)

df <- data.frame(sample=kdirs, group=substr(names(kdirs),1,1),
                 stringsAsFactors=FALSE)

so <- sleuth_prep(kdirs, df, ~group)

# function to extract transcript counts from sleuth/kallisto object
extract_sleuth <- function(slo, col="est_counts") {
  x <- sapply(slo$kal, function(u) u$abundance[[col]])
  colnames(x) <- basename(sapply(slo$kal,function(u) u$abundance$sample[1]))
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

# sum TPMs at gene level
tpmsG <- t(sapply(inds, function(u) colSums(tpms[u,,drop=FALSE])))

# create gene-level counts, scaled to the library sizes
countsG <- t(t(tpmsG)*colSums(countsT)/1e6)

# standard edgeR pipeline (other options exist, of course)
d <- DGEList(countsG, group=df$group)
d <- calcNormFactors(d)

mm <- model.matrix(~group,data=df)

# diagnostic 1
plotMDS(d)

# filter lowly expressed transcripts
k <- rowSums(cpm(d)>1) >= 3
d <- d[k,]

d <- estimateGLMCommonDisp(d, mm)
d <- estimateGLMTrendedDisp(d, mm)
d <- estimateGLMTagwiseDisp(d, mm)

# diagnostic 2
plotBCV(d)

cps <- cpm(d)
d$genes <- data.frame(gene_id=rownames(d), round(cps,1))

f <- glmFit(d, mm)
lrt <- glmLRT(f, coef=2)

# diagnostic 3
hist(lrt$table$PValue, 50)


options(width=160)
topTags(lrt)

# diagnostic 4
barplot( cps["ENSG00000227857",] )


