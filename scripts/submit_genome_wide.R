################
source("/hpc/users/giambc02/scripts/COLOC/pipeline/functions_coloc_pipeline.R")
source("/hpc/users/giambc02/scripts/COLOC/functions_coloc_likelihood_summary_integrated.R")

biom.dataset = "SCZ"
eqtl.fname = '/sc/orga/projects/epigenAD/coloc/data/CMC/eQTLs/matrixEQTL/DLPFC/ALLsummary.tab'
biom.fname = '/sc/orga/projects/roussp01a/resources/PGC2/daner_PGC_SCZ52_0513a.hq2'
NCASE = 35476
NCONTROL = 46839
N = NCASE + NCONTROL
type="cc"

biom.df = formatColoc(fname = biom.fname, type=type, N=N, Ncases=NCASE, info_filter=0.3, maf_filter=0.001, fread=T, eqtl=FALSE)
eqtl.df = formatColoc(fname = eqtl.fname, type="quant", N=NA, Ncases=NA, info_filter=0.3, maf_filter=0.001, fread=T, eqtl=TRUE)
p12=1e-6 # this is for coloc with set priors
outfolder = "/sc/orga/projects/epigenAD/coloc/results/SCZ_CMC/DLPFC/coloc.output.loci.perSNPpriors"
plot = FALSE
useBETA=TRUE
prefix= "SCZ_CMC"
save.coloc.output=FALSE

res = coloc.eqtl.biom(eqtl.df=eqtl.df, biom.df=biom.df, p12=p12, useBETA=TRUE, outfolder=outfolder, prefix=prefix, plot=plot, save.coloc.output=FALSE)

