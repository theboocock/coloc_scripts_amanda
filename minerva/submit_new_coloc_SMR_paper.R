match_snpid = TRUE
submit = TRUE

GWAS_DIR = "/sc/orga/projects/psychgen/resources/COLOC2/files/GWAS/"
#EQTL_DIR = "/sc/orga/projects/psychgen/resources/COLOC2/files/eQTLs/"
EQTL_DIR = "/sc/orga/projects/roussp01a/Claudia_TMP/COLOC2/temp/"
OUT_DIR = "/sc/orga/projects/psychgen/resources/COLOC2/"


eqtl_data = c("ACC", "DLPFC", "CMC_eGenes_1Mb", "CMC_eRNA_1Mb", "STARNET_AOR", "STARNET_Blood", "STARNET_LIV", "STARNET_MAM", "STARNET_SF", "STARNET_SKLM", "STARNET_VAF", "liver_Schadt_MAX", "liver_Schadt")
eqtl_files = paste(EQTL_DIR, eqtl_data, "_ALLsummary_CURATED.tab", sep="")
# Geuvadis? /sc/orga/projects/roussp01a/resources/SMR/eQTLs/Geuvadis/files/EUR373.gene.cis.FDR5.all.rs137.txt.gz 

#biom_data = c("AD", "CD", "ASD", "BMI", "CAD_ADD", "CAD_REC", "EXTRAVERSION", "EduYears", "FEMORAL_BMD", "FOREARM_BMD", "GIANT_HIP", "GIANT_HIPadjBMI", "GIANT_WC", "GIANT_WCadjBMI", "GIANT_WHR", "GIANT_WHRadjBMI", "HEIGHT", "IBD", "LIPIDS_HDL", "LIPIDS_LDL", "LIPIDS_TC", "LIPIDS_TG", "MI_ADD", "NEUROTICISM", "PGC2_BD", "PGC2_MDD", "PGC2_SCZ", "SPINE_BMD", "UC", "WellBeing_DS", "WellBeing_Neuroticism", "WellBeing_SWB", "cIMT", "cIMT_PLQ")
biom_data = c("GIANT_WCadjBMI", "GIANT_WHR", "GIANT_WHRadjBMI", "HEIGHT", "IBD", "LIPIDS_HDL", "LIPIDS_LDL", "LIPIDS_TC", "LIPIDS_TG", "MI_ADD", "NEUROTICISM", "PGC2_BD", "PGC2_MDD", "PGC2_SCZ", "SPINE_BMD", "UC", "WellBeing_DS", "WellBeing_Neuroticism", "WellBeing_SWB", "cIMT", "cIMT_PLQ")
#biom_files = paste(GWAS_DIR, biom_data, ".SMR_INPUT.tsv", sep="")
biom_files = paste(GWAS_DIR, biom_data, ".INPUT.tsv", sep="")

addCHRPOS = FALSE
if (addCHRPOS) {
   source("/hpc/users/giambc02/scripts/COLOC/pipeline/functions_coloc_pipeline.R")
   library(data.table)
   map_file = "/sc/orga/projects/psychgen/resources/COLOC2/files/GWAS/SNP_MAP_POSITION_HG19.tsv"
   map = fread(map_file, header=F)
   map = data.frame(map)
   names(map) = c("CHR", "POS", "SNPID")

   for (biom.fname in biom_files) {
     if (length(grep("BMI|EXTRAVERSION|EduYears|FEMORAL_BMD|FOREARM_BMD|SPINE_BMD|GIANT_HIP|GIANT_HIPadjBMI|GIANT_WC|GIANT_WCadjBM|GIANT_WHR|GIANT_WHRadjBMI|HEIGHT|LIPIDS_HDL|LIPIDS_LDL|LIPIDS_TC|LIPIDS_TG|NEUROTICISM|WellBeing_DS|WellBeing_Neuroticism|WellBeing_SWB|cIMT",  biom.fname))>0) type = "quant"
     if (length(grep("AD|CD|ASD|CAD_ADD|UC|IBD|MI_ADD|PGC2_BD|PGC2_MDD|PGC2_SCZ|PLQ",  biom.fname))>0) type = "cc"

     biom.df = formatColoc(fname = biom.fname, type=type, N=NA, Ncases=NA, info_filter=0, maf_filter=0, fread=T, eqtl=FALSE)
     biom.df$SNPID = as.character(biom.df$SNPID)
     #biom.df = biom.df[complete.cases(biom.df),]
     message("Number of rows before matching: ", nrow(biom.df))
     biom.df$CHR = map[match(biom.df$SNPID, map$SNPID), "CHR"]
     biom.df$POS = map[match(biom.df$SNPID, map$SNPID), "POS"]
     message("Number of rows before matching: ", nrow(biom.df[!is.na(biom.df$CHR),]))
     message("Number of rows before matching: ", nrow(biom.df[!is.na(biom.df$POS),]))
     write.table(biom.df, file=biom.fname, row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
   }
}

table_pairs = expand.grid(biom_files, eqtl_files)
names(table_pairs)=c("biom_files", "eqtl_files")

for (i in 1:nrow(table_pairs)) {

biom.fname = as.character(table_pairs$biom_files[i])
eqtl.fname = as.character(table_pairs$eqtl_files[i])

if (length(grep("BMI|EXTRAVERSION|EduYears|FEMORAL_BMD|FOREARM_BMD|SPINE_BMD|GIANT_HIP|GIANT_HIPadjBMI|GIANT_WC|GIANT_WCadjBM|GIANT_WHR|GIANT_WHRadjBMI|HEIGHT|LIPIDS_HDL|LIPIDS_LDL|LIPIDS_TC|LIPIDS_TG|NEUROTICISM|WellBeing_DS|WellBeing_Neuroticism|WellBeing_SWB|cIMT",  biom.fname))>0) type = "quant"
if (length(grep("AD|CD|ASD|CAD_ADD|UC|IBD|MI_ADD|PGC2_BD|PGC2_MDD|PGC2_SCZ|PLQ",  biom.fname))>0) type = "cc"

main_script = '/hpc/users/giambc02/scripts/COLOC/new_coloc_SMR_paper.R'
biom.name = gsub(".INPUT.tsv", "", basename(biom.fname))
eqtl.name = gsub("_ALLsummary_CURATED.tab", "", basename(eqtl.fname))
#eqtl.name = as.character(lapply(strsplit(as.character(eqtl.fname), "/", fixed=TRUE), "[", 9))
#message("Using biom ", biom.name, " and eQTL ", eqtl.name)
prefix = paste(biom.name, eqtl.name, sep="_")
print(prefix)
outfolder = paste(OUT_DIR, "/results/", prefix, "/", sep="")
      if (!file.exists (outfolder)) dir.create(outfolder)

scriptname=paste(OUT_DIR, "/scripts/Submit_main_script_", prefix, ".sh", sep="")

     write(file=scriptname, paste("#!/bin/bash
          #BSUB -J coloc_", prefix,"
          #BSUB -q alloc
          #BSUB -P acc_epigenAD
          #BSUB -n 5
          #BSUB -R span[hosts=1]
          #BSUB -R 'rusage[mem=4000]'
          #BSUB -W 30:00 
          #BSUB -L /bin/bash
          #BSUB -oo ", OUT_DIR, "/log/", prefix, ".out
          #BSUB -eo ", OUT_DIR, "/log/", prefix, ".err", sep=""), append=F)
          write(file=scriptname, paste("Rscript", main_script, biom.fname, eqtl.fname, type, outfolder, prefix, sep=" "), append=T)

          message("Submit script: ", scriptname)
          if (submit) {
           system(paste("bsub < ", scriptname, sep=""))
          }


}
 
