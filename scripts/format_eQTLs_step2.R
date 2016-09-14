library(data.table)
source("/hpc/users/giambc02/scripts/COLOC/pipeline/functions_coloc_pipeline.R")

############################################################
# input
indir = "/sc/orga/projects/psychgen/resources/COLOC2/files/eQTLs/"
#eqtl_files= list.files(indir,pattern="*.tab",full.names=T)
eqtl.datasets = c("ACC", "DLPFC", "CMC_eGenes_1Mb", "CMC_eRNA_1Mb", "Liver_Schadt")
eqtl_files = c("/sc/orga/projects/epigenAD/coloc/data/CMC/eQTLs/matrixEQTL/ACC/ALLsummary2.tab", "/sc/orga/projects/epigenAD/coloc/data/CMC/eQTLs/matrixEQTL/DLPFC/ALLsummary2.tab", "/sc/orga/projects/epigenAD/coloc/data/CMC/eQTLs/matrixEQTL/eRNA_NoGenes_1e6/ALLsummary2.tab", "/sc/orga/projects/epigenAD/coloc/data/CMC/eQTLs/matrixEQTL/eRNA_OnlyGenes/ALLsummary2.tab", "/sc/orga/projects/epigenAD/coloc/data/liver_Schadt/eQTLs/matrixEQTL/Liver/ALLsummary2.tab")

gtex_files = list.files(indir, pattern="*.coloc.txt", full.names=T)
gtex_annotations = read.table("/sc/orga/projects/psychgen/resources/COLOC2/data/gtex_annot.txt", header=T)
gtex_eqtl.datasets = basename(gtex_files)

# output
odir = "/sc/orga/projects/psychgen/resources/COLOC2/files/eQTLs2/"
###########################################################
other_eqtl=FALSE
if (other_eqtl) {
for (i in 5:length(eqtl_files)) {
     eqtl.fname = eqtl_files[i]
     eqtl.dataset = eqtl.datasets[i]
     outFile= paste(odir, eqtl.dataset, "_formatted", sep="")
     eqtl_sample_size = NULL # in the data

     cat("\n-----------------------------------------------------------------\nSTUDY:",eqtl.dataset,"file:",eqtl.fname,"output:", outFile,"\n-----------------------------------------------------------------\n")

     eqtl.df = tryCatch(formatColoc(fname = eqtl.fname, type="quant", N=eqtl_sample_size, Ncases=NULL, info_filter=0, maf_filter=0, fread=T, eqtl=TRUE), error=function(e) NULL )

     if (is.null(eqtl.df)) stop("******************************************!!!!! Study could not be processed")
 
     # take out "chr" from SNPID
     eqtl.df$SNPID= gsub("chr", "", eqtl.df$SNPID)
     # take out "chr" from CHR
     eqtl.df$CHR= gsub("chr", "", eqtl.df$CHR)
     n_cutoff = quantile(eqtl.df$N, c(.50))
     info = data.frame(data=eqtl.dataset, NsnpsCompleteData_before_filter = nrow(eqtl.df[complete.cases(eqtl.df),]), PVAL_range = paste(range(eqtl.df$PVAL), collapse=" , "), NsamplesMax= max(eqtl.df$N), NsamplesMin=min(eqtl.df$N), Ncut=as.numeric(n_cutoff), filtered_snps=nrow(eqtl.df)-nrow(eqtl.df[eqtl.df$N >=n_cutoff,]))

     write.table(info, file=paste(odir, "_info.txt", sep=""), append = TRUE, sep="\t")
     write.table(eqtl.df, file=outFile, row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
     write.table(eqtl.df[eqtl.df$N >=n_cutoff,], file=paste(odir, eqtl.dataset, "_formatted_Nfiltered", sep=""), row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
     message("File for ", eqtl.dataset, " saved in ", outFile)
   }
}
###########################################################
gtex=TRUE
if (gtex) {
for (eqtl.fname in gtex_files[13:44]) {

     eqtl.dataset = basename(eqtl.fname)
     outFile= paste(odir, eqtl.dataset, "_formatted", sep="")

#if(eqtl.fname %in% gtex_files){
    eqtl_sample_size = gtex_annotations$RNASeq_AND_Genotyped_samples[which(gtex_annotations$file.name == basename(eqtl.fname))]
#}else{
#    eqtl_sample_size = NULL
#}

     cat("\n-----------------------------------------------------------------\nSTUDY:",eqtl.dataset,"file:",eqtl.fname,"output:", outFile,"\n-----------------------------------------------------------------\n")

     eqtl.df = tryCatch(formatColoc(fname = eqtl.fname, type="quant", N=eqtl_sample_size, Ncases=NULL, info_filter=0, maf_filter=0, fread=T, eqtl=TRUE), error=function(e) NULL )

     if (is.null(eqtl.df)) stop("******************************************!!!!! Study could not be processed")

     # take out "chr" from SNPID
     eqtl.df$SNPID= gsub("chr", "", eqtl.df$SNPID)
     # take out "chr" from CHR
     eqtl.df$CHR= gsub("chr", "", eqtl.df$CHR)
     n_cutoff = quantile(eqtl.df$N, c(.50))
     info = data.frame(data=eqtl.dataset, NsnpsCompleteData_before_filter = nrow(eqtl.df[complete.cases(eqtl.df),]), PVAL_range = paste(range(eqtl.df$PVAL), collapse=" , "), NsamplesMax= max(eqtl.df$N), NsamplesMin=min(eqtl.df$N), Ncut=as.numeric(n_cutoff), filtered_snps=nrow(eqtl.df)-nrow(eqtl.df[eqtl.df$N >=n_cutoff,]))

     write.table(info, file=paste(odir, "_info.txt", sep=""), append = TRUE, sep="\t")
     write.table(eqtl.df, file=outFile, row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
     #write.table(eqtl.df[eqtl.df$N >=n_cutoff,], file=paste(odir, eqtl.dataset, "_formatted_Nfiltered", sep=""), row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
     message("File for ", eqtl.dataset, " saved in ", outFile)

   }
}
