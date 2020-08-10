#!/usr/bin/Rscript
library(dada2)
library(data.table)

#Read the path to fastq files from command line arugments
args=(commandArgs(TRUE))
path = args[1]

#List of forward reads
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
#List of reverse reads
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

# Extract sample names from filenames
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Filtered output filenames for forward reads
filtFs <- file.path(path, "trimmed", paste0(sample.names, "_1_filt.fastq.gz"))

#Filtered output filenames for reverse reads
filtRs <- file.path(path, "trimmed", paste0(sample.names, "_2_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

#trim primers from the fastq files (Forward trim (index+primer =27) and Reverse trim (index+primer =27)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(27,28), compress=TRUE, multithread=TRUE)

#Create Error models
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Infer unique sequences
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Merged Sequences
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)
write.csv(seqtab, file.path(path,"DADA2SeqTable.csv"))


#Generate merged files for each sample
Data <- fread(file.path(path,"DADA2SeqTable.csv"), header=FALSE, sep=",")
names(Data) <- c("SampleName", paste0("seq", 1:(length(Data)-1)))
SampleNames <- Data$SampleName[-1]

system(paste0("mkdir ",file.path(path,"Merged")))

#split the sequenance by sample
for(sample in SampleNames){
    SampleData <-rbind(Data[1,],Data[ which(Data$SampleName== sample)])
    SampleData <- t(SampleData[,-1])
    SampleData <- SampleData[which(SampleData[,2] != 0),]
    fileConn<-file(file.path(path, "Merged", paste0(sample,".fna")))
    lines <- list()
    for (seq in rownames(SampleData)){
      line <- paste(rep(paste0(">",seq, "\n" ,SampleData[seq,1], "\n"), SampleData[seq,2]),  collapse = '')
      lines = append(lines, list(line))
    }
    lines <- paste(lines, collapse = '')
    writeLines(lines, fileConn)
    close(fileConn)
}

#Save a summary of sequence abundance in each file
Abundance <- Data[-1,]
write.csv(Abundance, file.path(path,"SequenceAbundanceSummary.csv"), row.names = FALSE)


#Save sequence id file
Seqid<- data.frame(t(Data[1,-1]))
names(Seqid) <-"Seq"
write.csv(Seqid, file.path(path,"SequenceID.csv"), row.names = TRUE)

 
