args <- commandArgs(trailingOnly=TRUE)
biom.fname=args[1]
eqtl.fname=args[2]
type=args[3]
outfolder = args[4]
prefix = args[5]
n_eqtl = args[6]

source("~/psychgen/resources/COLOC2/COLOC_scripts/scripts/functions_coloc_pipeline.R")
source("~/psychgen/resources/COLOC2/COLOC_scripts/scripts/functions_coloc_likelihood_summary_integrated.R")

# filter by maf and imp quality when doing colocalization
# liver eQTL is messed up: data= fread(fname, select = colsAll[[1]], col.names=names(colsAll[[1]]), colClasses=c(SNPID="character", ProbeID="character", BETA="numeric", PVAL="numeric", CHR="character", POS="numeric", F="numeric", SE="numeric", N="numeric"))
# data= fread(fname, stringsAsFactors=FALSE)
biom.df = formatColoc(fname = biom.fname, type=type, N=NA, Ncases=NA, info_filter=0, maf_filter=0, fread=T, eqtl=FALSE)
eqtl.df = formatColoc(fname = eqtl.fname, type="quant", N=n_eqtl, Ncases=NA, info_filter=0, maf_filter=0, fread=T, eqtl=TRUE)
biom.df$type = type
p12=1e-6 # this is for coloc with set priors
plot = TRUE
useBETA=TRUE
save.coloc.output=TRUE

res = coloc.eqtl.biom(eqtl.df=eqtl.df, biom.df=biom.df, p12=p12, useBETA=TRUE, outfolder=outfolder, prefix=prefix, plot=plot, save.coloc.output=TRUE, match_snpid=FALSE)

