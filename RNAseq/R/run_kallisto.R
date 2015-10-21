
cDNA.fasta <- "../anno/transcripts.fa"
index <- "../indices/kallisto/kallisto_index"

fq <- dir("../FASTQ",".fq.gz$", full=TRUE)
names(fq) <- gsub("_[12].fq.gz","",basename(fq))

fqs <- split(fq, names(fq))

for(i in 1:length(fqs)) {
  cmd <- sprintf("kallisto quant -i %s -o %s %s %s",
                 index, names(fqs)[i], fqs[[i]][1], fqs[[i]][2])
  cat(cmd,"\n"); 
  st <- system.time(system(cmd))
  print(st)
}


