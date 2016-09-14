match_snpid = F 
#submit = TRUE
submit=F
GWAS_DIR = "/sc/orga/projects/psychgen/resources/COLOC2/files/GWAS/"
#EQTL_DIR = "/sc/orga/projects/psychgen/resources/COLOC2/files/eQTLs/"
OUT_DIR = "/sc/orga/projects/psychgen/resources/COLOC2/temp_results"


eqtl_files= list.files("/sc/orga/projects/psychgen/resources/COLOC2/files/eQTLs/",pattern="*.tab",full.names=T)
eqtl_data = basename(eqtl_files) 

gtex_files = list.files("/sc/orga/projects/psychgen/resources/COLOC2/files/GTEX/", pattern="*.coloc.txt", full.names=T)
gtex_data = basename(gtex_files)

gtex_annotations = read.table("/sc/orga/projects/psychgen/resources/COLOC2/data/gtex_annot.txt", header=T)



# Geuvadis? /sc/orga/projects/roussp01a/resources/SMR/eQTLs/Geuvadis/files/EUR373.gene.cis.FDR5.all.rs137.txt.gz 

biom_data =c("CARDIoGRAMplusC4D", "CARDIoGRAMplusC4D_REC", "MI", "DAIGRAM", "RA", "ASD", "UC", "CD", "IBD", "SCZ", "BIP", "AD", "PLAQUE", "CIMT", "HDL", "LDL", "TG", "TC", "EXTRAVERSION", "NEUROTICISM", "EduYears", "FEMORAL", "FOREARM", "SPINE", "WellBeing_DS", "WellBeing_Neuroticism", "WellBeing_SWB", "BMI", "HEIGHT", "GIANT_HIP", "GIANT_HIPadjBMI", "GIANT_WC", "GIANT_WCadjBMI", "GIANT_WHR", "GIANT_WHRadjBMI", "MAGIC_Manning_et_al_FastingGlucose_MainEffect", "MAGIC_HbA1C", "MAGIC_Manning_et_al_lnFastingInsulin_MainEffect", "MAGIC_ln_HOMA-B", "MAGIC_ln_HOMA-IR", "URATE")
  
#biom_files = paste(GWAS_DIR, biom_data, ".SMR_INPUT.tsv", sep="")
biom_files = paste(GWAS_DIR, biom_data, "_formatted", sep="")
eqtl_files = c(eqtl_files, gtex_files)
table_pairs = expand.grid(biom_files, eqtl_files)
names(table_pairs)=c("biom_files", "eqtl_files")


for (i in 1:nrow(table_pairs)) {

biom.fname = as.character(table_pairs$biom_files[i])
eqtl.fname = as.character(table_pairs$eqtl_files[i])
if(eqtl.fname %in% gtex_files){
    eqtl_sample_size = gtex_annotations$RNASeq_AND_Genotyped_samples[which(gtex_annotations$file.name == basename(eqtl.fname))] 
}else{
    eqtl_sample_size = NA
}
if (length(grep("CIMT|HDL|LDL|TG|TC|EXTRAVERSION|NEUROTICISM|EduYears|FEMORAL|FOREARM|SPINE|WellBeing_DS|WellBeing_Neuroticism|WellBeing_SWB|BMI|HEIGHT|GIANT_HIP|GIANT_HIPadjBMI|GIANT_WC|GIANT_WCadjBMI|GIANT_WHR|GIANT_WHRadjBMI|MAGIC_Manning_et_al_FastingGlucose_MainEffect|MAGIC_HbA1C|MAGIC_Manning_et_al_lnFastingInsulin_MainEffect|MAGIC_ln_HOMA-B|MAGIC_ln_HOMA-IR|URATE",  biom.fname))>0) type = "quant"
if (length(grep("CARDIoGRAMplusC4D|CARDIoGRAMplusC4D_REC|MI|DAIGRAM|RA|ASD|UC|CD|IBD|SCZ|BIP|AD|PLAQUE", biom.fname))>0) type = "cc"

main_script = '/hpc/users/giambc02/scripts/COLOC/new_coloc_SMR_paper.R'
biom.name = gsub(".SMR_INPUT.tsv", "", basename(biom.fname))
eqtl.name = gsub(".tab", "", basename(eqtl.fname))
#eqtl.name = as.character(lapply(strsplit(as.character(eqtl.fname), "/", fixed=TRUE), "[", 9))
message("Using biom ", biom.name, " and eQTL ", eqtl.name)
prefix = paste(biom.name, eqtl.name, sep="_")
outfolder = paste(OUT_DIR, "/results/", prefix, "/", sep="")
      if (!file.exists (outfolder)) dir.create(outfolder)

scriptname=paste(OUT_DIR, "/scripts/Submit_main_script_", prefix, ".sh", sep="")

     write(file=scriptname, paste("#!/bin/bash
          #BSUB -J coloc_", prefix,"
          #BSUB -q premium 
          #BSUB -P acc_psychgen 
          #BSUB -n 20
          #BSUB -R span[hosts=1]
          #BSUB -R 'rusage[mem=50000]'
          #BSUB -W 30:00 
          #BSUB -L /bin/bash
          #BSUB -oo ", OUT_DIR, "/log/", prefix, ".out
          #BSUB -eo ", OUT_DIR, "/log/", prefix, ".err", sep=""), append=F)
          write(file=scriptname, paste("Rscript", main_script, biom.fname, eqtl.fname, type, outfolder, prefix,eqtl_sample_size, sep=" "), append=T)

          message("Submit script: ", scriptname)
          if (submit) {
           system(paste("bsub < ", scriptname, sep=""))
          }
}
 
