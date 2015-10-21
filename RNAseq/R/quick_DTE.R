
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

counts <- extract_sleuth(so,"est_counts")

# standard edgeR pipeline (other options exist, of course)
d <- DGEList(counts, group=df$group)
d <- calcNormFactors(d)

mm <- model.matrix(~group,data=df)

d <- DGEList(counts, group=df$group)
d <- calcNormFactors(d)

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
d$genes <- data.frame(transcript_id=rownames(d), round(cps,1))

f <- glmFit(d, mm)
lrt <- glmLRT(f, coef=2)

# diagnostic 3
hist(lrt$table$PValue, 50)


options(width=160)
topTags(lrt)

# diagnostic 4
barplot( cps["ENST00000361886",] )


