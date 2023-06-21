library(ExomeDepth)
library(methods)
args <- (commandArgs(trailingOnly = TRUE))
#El path al bedfile modificado se sustituye en este caso con el argumento número $10 de bash que se
#señala como args[10] es decir, el número 10 en la lista de argumentos:
arg10 <- args[10]
print(args)

#==============================INPUTS=========================#
result_files <- as.character(args[11])
target.file <- read.table(args[10], header = TRUE, sep = "\t") # set path to BED file
reference.file <- args[5]
bamFile_control <- list.files(args[2], patter=".bam$", full.names=TRUE)
#==============================================================#

my.refcounts <- getBamCounts(bed.frame = target.file,
                             bam.files = bamFile_control,
                             include.chr = FALSE,
                             referenceFasta = reference.file)
output = paste(result_files, "ref_counts.txt", sep = "/")
print(head(my.refcounts))
write.table(my.refcounts, file = output)

bamFile_muestra <- list.files(args[1], patter=".bam$", full.names=TRUE)
print(head(bamFile_muestra))

output = paste(result_files, "orden_proceso.txt", sep = "/")
write.table(bamFile_muestra, file= output)

my.counts <- getBamCounts(bed.frame = target.file,
                          bam.files = bamFile_muestra,
                          include.chr = FALSE,
                          referenceFasta = reference.file)

output = paste(result_files, "my_counts.txt", sep = "/")
write.table(my.counts,file = output)

my.counts.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')

#my.counts.dafr$chromosome <- gsub(as.character(my.counts.dafr$space),
#                                 pattern = 'chr',
#                                replacement = '') ##remove the annoying chr letters

print(head(my.counts.dafr))
samplecounts.mat<-as.matrix(my.counts.dafr[,grep(names(my.counts.dafr),pattern='*.bam')])
nsamples<-ncol(samplecounts.mat)
my.refcounts.dafr <- as(my.refcounts[, colnames(my.refcounts)], 'data.frame')
#my.refcounts.dafr$chromosome <- gsub(as.character(my.refcounts.dafr$space),
#                                    pattern = 'chr',
#                                   replacement = '') ##remove the annoying chr letters
my.ref.samples<-colnames(my.refcounts.dafr)[7:(ncol(my.refcounts.dafr)-1)]
print(head(my.refcounts.dafr))
print('=====================================================')
#loop over samples in my.counts
for (i in 1:nsamples) {{
  my.current.samplename <-colnames(my.counts.dafr[4+i])
  print(my.current.samplename)
  #  my.current.sample<-samplecounts.mat$my.current.sample
  my.reference.set <- as.matrix(my.refcounts.dafr[,my.ref.samples])
  my.choice<-select.reference.set(test.counts=samplecounts.mat[,i],
                                  reference.counts=(my.reference.set),
                                  bin.length=(my.counts.dafr$end - my.counts.dafr$start)/1000,
                                  n.bins.reduced = 10000)
  
  print(my.choice[[1]])
  my.matrix <- as.matrix( my.refcounts.dafr[, my.choice$reference.choice, drop = FALSE])
  my.reference.selected <- apply(X = my.matrix,
                                 MAR = 1,
                                 FUN = sum)
  
  #CNV calling
  all.exons <- new('ExomeDepth',
                   test = samplecounts.mat[,i],
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')
  
  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        start = my.counts.dafr$start,
                        end = my.counts.dafr$end,
                        chromosome = my.counts.dafr$chromosome,
                        name = my.counts.dafr$exon)
  
  
  output.file <- paste(my.current.samplename,'.exome_calls.csv',sep = ",")
  save(all.exons,file=paste(my.current.samplename,".all.exons.txt",sep = "\t"))
  write.csv(file = output.file,
            x = all.exons@CNV.calls,
            row.names = FALSE)
}
  ###Imprimimos en un artchivo de salida las CNVs encontradas######
  if(nrow(all.exons@CNV.calls)>0){
    output.file <- paste0(i ,".cnv.txt")
    write.table(file = file.path(result_files,output.file),x = cbind(i,all.exons@CNV.calls),row.names = FALSE,quote=F,sep="\t") > result_files
  }
}
