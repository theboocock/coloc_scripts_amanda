args <- commandArgs(trailingOnly=TRUE)
biom.fname=args[1]
eqtl.fname=args[2]
outfolder = args[3]
prefix = args[4]

source("/sc/orga/projects/epigenAD/coloc/coloc2_gitrepo/coloc_scripts/scripts/functions_coloc_pipeline.R")
source("/sc/orga/projects/epigenAD/coloc/coloc2_gitrepo/coloc_scripts/scripts/functions_coloc_likelihood_summary_integrated.R")

# filter by maf and imp quality when doing colocalization
# liver eQTL is messed up: data= fread(fname, select = colsAll[[1]], col.names=names(colsAll[[1]]), colClasses=c(SNPID="character", ProbeID="character", BETA="numeric", PVAL="numeric", CHR="character", POS="numeric", F="numeric", SE="numeric", N="numeric"))
# data= fread(fname, stringsAsFactors=FALSE)
#biom.df = formatColoc(fname = biom.fname, type=type, N=NA, Ncases=NA, info_filter=0, maf_filter=0, fread=T, eqtl=FALSE)
#eqtl.df = formatColoc(fname = eqtl.fname, type="quant", N=NA, Ncases=NA, info_filter=0, maf_filter=0, fread=T, eqtl=TRUE)

library(data.table)
biom.df = data.frame(fread(biom.fname, header=T))
eqtl.df = data.frame(fread(eqtl.fname, header=T))

p12=1e-6 # this is for coloc with set priors
plot = FALSE
useBETA=TRUE
save.coloc.output=FALSE

# estimate per locus likelihoods (for now also outputs the ppa using set priors
res1 = coloc.eqtl.biom(eqtl.df=eqtl.df, biom.df=biom.df, p12=p12, useBETA=useBETA, outfolder=outfolder, prefix=prefix, plot=plot, save.coloc.output=save.coloc.output, match_snpid=FALSE)


# optimization and posteriors using these priors estimates
outfname = paste(outfolder, prefix, '_summary.tab', sep='')
res.all = read.table(outfname, header=T, sep="\t")
# split the lklds
out <- strsplit(as.character(res.all$coloc.old.pval.lkl),',')
t= data.frame(res.all, do.call(rbind, out))
names(t)[(ncol(t)-4):ncol(t)] = c("lH0.abf", "lH1.abf", "lH2.abf", "lH3.abf", "lH4.abf") 

res2 = est_lkl(t)


  #res.all <- res.all[with(res.all, order(pp4, decreasing=T)),]
   #outfname = paste(outfolder, prefix, '_summary.tab', sep='')
   write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')

   # If Gene.name is missing, use ensemblID instead, then try to retrieve name from biomaRt. 
   if (length(res.all$ProbeID[grep("ENSG", res.all$ProbeID)]) >0  & !("Gene.name" %in% names(res.all))) addGeneName = TRUE
   addGeneName= FALSE
   if (addGeneName) {
   res.all$Gene.name = res.all$ProbeID
   # TODO ANnotation output filewith gene name
   biomart=FALSE # it doesn't work sometimes -- cannot connect etc
      if (biomart) {
      library(biomaRt)
      #if (length(res.all$Gene.name[grep("ENSG", res.all$Gene.name)]) >0 ) {
        mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
        res.gn <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = as.character(res.all$Gene.name[grep("ENSG", res.all$Gene.name)]), mart = mart)
        res.gn = res.gn[res.gn$hgnc_symbol!="",]
        res.all$Gene.name = res.gn[match(res.all$ProbeID, res.gn$ensembl_gene_id),"hgnc_symbol"]
        #res.all$Gene.name[which(res.all$Gene.name %in% res.gn$ensembl_gene_id)]= res.gn[match(res.all$Gene.name[which(res.all$Gene.name %in% res.gn$ensembl_gene_id)], res.gn$ensembl_gene_id), "hgnc_symbol"]
     } else {
        geneFileNames = "/sc/orga/projects/roussp01a/resources/Ensembl2HGNC/ENSEMBL_v70_TO_HGNC.tsv"
        genes = read.table(geneFileNames, header=F, stringsAsFactors=FALSE, col.names=c("ensembl_gene_id", "hgnc_symbol"))
        res.all$Gene.name = genes[match(res.all$Gene.name, genes$ensembl_gene_id), "hgnc_symbol"]
    }
   }
   write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')
   write.table(x = removed_snp_list, file = out_removed_snps, row.names = FALSE, quote = FALSE, sep = '\t')
   return(res.all)
}

