##################################################################
##### Instalar y cargar la librería de panelcn.mops
library(panelcn.mops)
args <- (commandArgs(trailingOnly = TRUE))
#################### CARGAR ARCHIVO BED ##########################
##### Read bedfile y countWindows
bed <- args[12]
countWindows <- getWindows(bed, chr = TRUE)
##### Read Bamfiles - Test

message("Before line 10....")
message(args[1])
BAMFiles <- list.files(args[1], patter=".bam$", full.names=TRUE)
message("After line 10....")
message(BAMFiles)
test <- countBamListInGRanges(countWindows = countWindows,
                              bam.files = BAMFiles, read.width = 150)
message("After countBamListInGRanges ...")
##### Read Bamfiles - Control
BAMFiles_control <- list.files(args[2] ,patter=".bam$", full.names=TRUE)
control <- countBamListInGRanges(countWindows = countWindows,
                                 bam.files = BAMFiles_control, read.width = 150)
#####
panelcnmopstest <- test
elementMetadata(panelcnmopstest) <- cbind(elementMetadata(panelcnmopstest),
                                          elementMetadata(control))
##### Read sampleNames
sampleNames <- colnames(elementMetadata(test))
print(sampleNames)
##### Select Genes
selectedGenes <- unlist(strsplit(args[13], split=","))
print(selectedGenes)
##### RESULTS, con este loop te genera un archivo con terminación "CNVPJA.txt" los CNV
##### Puedes cambiarlo al nombre de tus archivos.
##### resultlistPJA y resulttablePJA, lo puedes cambiar de acuerdo a la teminación
##### que tu quieras, por ejemplo resultlistGuatemala o lo puedes dejar así
##### Porque trabaja en el loop y al fina te da los archivos txt que se ocupan
##### en el siguiente comando
for (i in c(1:length(BAMFiles)))
{resultlist <- runPanelcnMops(panelcnmopstest, testiv=c(i),
                                 countWindows = countWindows, selectedGenes = selectedGenes,
                                 I = c(0.025, 0.57, 1, 1.46,2), normType = "quant", sizeFactor = "quant",
                                 qu = 0.25, quSizeFactor = 0.75, norm = 1, priorImpact = 1, minMedianRC = 30,
                                 maxControls = 25, sex = "female")
message("before resulttable...")
resulttable <- createResultTable(resultlist, panelcnmopstest,
                                    countWindows = countWindows, selectedGenes = selectedGenes, sampleNames)
message("after resulttable...")
message("before write table")
write.table(resulttable, file = paste(BAMFiles[i], args[14], sep =""))
message("after write table...")
}
##### Join the individual results of each sample in a single .txt file
##### Establezca la ruta donde se guardaron los archivos txt generado del
##### loop y con los siguientes comandos unes todos los archivos
##### en una sola hoja de txt

#setwd(args[15])
files <- list.files(pattern=".txt")
allresults <- lapply(files, function(I)read.table(I))
out <- do.call(rbind, allresults)

output_save <- paste(args[15],'/panelcnmops_CNVresults.txt', sep='')

write.table(out, output_save, sep = "\t")
