
# 1. Run
# 2. Run 

library(DEXSeq)

ff <- "../data/anno/transcripts_kallisto_filtered.gff"

dxf <- dir("../data/dexseq",".counts$",full=TRUE)
names(dxf) <- basename(dxf)

df <- data.frame(countfile=dxf, group=substr(names(dxf),1,1),
                 stringsAsFactors=FALSE)

sample_table <- data.frame(condition = conditions)

dxd <- DEXSeqDataSetFromHTSeq(countfiles=dxf, sampleData=df, 
                              design=~ sample + exon + group:exon, flattenedfile=ff)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM=MulticoreParam(workers=4))
dxd <- testForDEU(dxd, BPPARAM=MulticoreParam(workers=4))
res <- DEXSeqResults(dxd)


pdf("dexseq_top.pdf", w=10,h=10)
plotDEXSeq(res, "ENSG00000016490", fitExpToVar="group", 
           displayTranscripts=TRUE); 
dev.off()

