library(data.table)

other_eqtl=TRUE
if (other_eqtl) {
#for (eqtl.dataset in c("ACC", "DLPFC", "CMC_eGenes_1Mb", "CMC_eRNA_1Mb", "Liver_Schadt")) {
for (eqtl.dataset in c("Liver_Schadt")) {

 if (eqtl.dataset == "ACC") {
     dir = "/sc/orga/projects/epigenAD/coloc/data/CMC/eQTLs/matrixEQTL/ACC/"
     eqtl.fnames = list.files(dir, pattern="_all_pval.tab", full.names=T)
     NCASE = NULL
     N=NULL
 }

 if (eqtl.dataset == "DLPFC") {
     dir = "/sc/orga/projects/epigenAD/coloc/data/CMC/eQTLs/matrixEQTL/DLPFC/"
     eqtl.fnames = list.files(dir, pattern="_all_pval.tab", full.names=T)
     NCASE = NULL
     N=NULL
 }

 if (eqtl.dataset == "CMC_eRNA_1Mb") {
     dir = "/sc/orga/projects/epigenAD/coloc/data/CMC/eQTLs/matrixEQTL/eRNA_NoGenes_1e6/"
     eqtl.fnames = list.files(dir, pattern="_all_pval.tab", full.names=T)
     NCASE = NULL
     N=NULL
 }

 if (eqtl.dataset == "CMC_eGenes_1Mb") {
     dir = "/sc/orga/projects/epigenAD/coloc/data/CMC/eQTLs/matrixEQTL/eRNA_OnlyGenes/"
     eqtl.fnames = list.files(dir, pattern="_all_pval.tab", full.names=T)
     NCASE = NULL
     N=NULL
 }

 if (eqtl.dataset == "Liver_Schadt") {
     dir = "/sc/orga/projects/epigenAD/coloc/data/liver_Schadt/eQTLs/matrixEQTL/Liver/"
     eqtl.fnames = list.files(dir, pattern="_all_pval.tab", full.names=T)
     NCASE = NULL
     N=NULL
 }


# make sure have all chroms
if (!all(paste("chr", 1:22, "_all_pval.tab", sep="") %in% basename(eqtl.fnames))) stop("Some chromosome files are missing from folder")
data = data.frame()

if (eqtl.dataset!="Liver_Schadt") {
for (i in 1:22) {
    fn = paste(dir, "chr", i, "_all_pval.tab", sep="")
    d = data.frame(fread(fn, header=T))
    d = d[,c("SNPID", "ProbeID", "beta", "t.stat", "PVAL", "FDR", "CHR","POS","F","REF","ALT","se.beta","N", "gene.chromosome", "gene.position.start", "gene.position.end")]
    message(i, " Ngenes: ", length(unique(d$ProbeID)))
    data=rbind.data.frame(data, d)
}
}

if (eqtl.dataset=="Liver_Schadt") {
for (i in 1:22) {
    fn = paste(dir, "chr", i, "_all_pval.tab", sep="")
    d = data.frame(fread(fn, header=T))
    d = d[,c("SNPID", "ProbeID", "beta", "t.stat", "PVAL", "CHR", "POS", "F", "allele.1", "allele.2","se.beta","N", "gene.chromosome", "ensemblID", "Gene.name","gene.position.start","gene.position.end")]
    message(i, " Ngenes: ", length(unique(d$ProbeID)))
    data=rbind.data.frame(data, d)
}
}

oFile = paste(dir, "ALLsummary2.tab", sep="")
write.table(data, file=oFile, row.names = FALSE, quote = FALSE, col.names = TRUE)


}

}

starnet=FALSE
if (starnet) {

for (eqtl.dataset in c("AOR", "Blood", "LIV", "MAM", "SF",  "SKLM",  "VAF")) {
     dir = paste("/sc/orga/projects/epigenAD/coloc/data/STARNET/eQTLs/matrixEQTL/", eqtl.dataset, "/", sep="")
     eqtl.fnames = list.files(dir, pattern="_all_pval.tab", full.names=T)
     if (!all(paste("chr", 1:22, "_all_pval.tab", sep="") %in% basename(eqtl.fnames))) stop("Some chromosome files are missing from folder")
     data = data.frame()
     for (i in 1:22) {
        fn = paste(dir, "chr", i, "_all_pval.tab", sep="")
        d = data.frame(fread(fn, header=T))
        d = d[,c("SNPID", "input_name", "ProbeID", "BETA", "T_STAT", "PVAL", "SE", "CHR", "POS", "REF", "ALT", "info.impute2", "F", "N", "gene.chromosome", "gene.position.start", "gene.position.end")]
       names(d)[1] = "chrpos"
       names(d)[2] ="SNPID"
    message(eqtl.dataset, ": chr", i, " Ngenes: ", length(unique(d$ProbeID)))
    data=rbind.data.frame(data, d)
}
oFile = paste(dir, "ALLsummary2.tab", sep="")
write.table(data, file=oFile, row.names = FALSE, quote = FALSE, col.names = TRUE)

}
}

