library(GenomicRanges)
library(rtracklayer)

bw.file="http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/E003-H3K4me3.fc.signal.bigwig"
k4me3=import(bw.file, 
       which=seqinfo(BigWigFile(bw.file))["chr20"])

bw.file="http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/E003-DNase.fc.signal.bigwig"
dnase=import(bw.file, 
             which=seqinfo(BigWigFile(bw.file))["chr20"])

bw.file="http://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/E003-H3K4me1.fc.signal.bigwig"
k4me1=import(bw.file, 
             which=seqinfo(BigWigFile(bw.file))["chr20"])
k4me3

export.bw(k4me3,"H1.ESC.H3K4me3.chr20.bw")
export.bw(k4me1,"H1.ESC.H3K4me1.chr20.bw")
export.bw(dnase,"H1.ESC.dnase.chr20.bw")
