## prepare annotation (at command prompt)
#python dexseq_prepare_annotation.py --aggregate='no' ../data/anno/transcripts_kallisto_filtered.gtf ../data/anno/transcripts_kallisto_filtered.gff

bf <- dir("../data/BAM","_s.bam$", full=TRUE)
names(bf) <- gsub("_STAR_s.bam","",basename(bf))

gff <- "../data/anno/transcripts_kallisto_filtered.gff"


for(i in 1:length(bf)) {
  cmd <- sprintf("python dexseq_count.py --format='bam' --paired='yes' --order='pos' %s %s ../data/dexseq/%s.counts",
                 gff, bf[i], names(bf)[i])
  cat(cmd,"\n")
  system(cmd,wait=FALSE)
}

