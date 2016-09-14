addCHRPOS_fn <- function(biom.df) {
     biom.df$SNPID = as.character(biom.df$SNPID)
     #biom.df = biom.df[complete.cases(biom.df),]
     message("Number of rows before matching: ", nrow(biom.df))
     biom.df$CHR = map[match(biom.df$SNPID, map$SNPID), "CHR"]
     biom.df$POS = map[match(biom.df$SNPID, map$SNPID), "POS"]
     message("Number of rows after matching: ", nrow(biom.df[!is.na(biom.df$CHR),]))
     message("Number of rows after matching: ", nrow(biom.df[!is.na(biom.df$POS),]))
     if (nrow(biom.df[is.na(biom.df$CHR),])>0) {
           message("There are still missing CHR")
           # try with other maps
          biom.df$CHR2 = map2[match(biom.df$SNPID, map2$SNPID), "CHR"]
          biom.df$CHR2 = gsub("chr", "", biom.df$CHR2)
          biom.df$POS2 = map2[match(biom.df$SNPID, map2$SNPID), "POS"]
          biom.df$CHR = ifelse(is.na(biom.df$CHR), biom.df$CHR2, biom.df$CHR)
          biom.df$POS = ifelse(is.na(biom.df$POS), biom.df$POS2, biom.df$POS)
          message("Number of rows after matching: ", nrow(biom.df[!is.na(biom.df$CHR),]))
          message("Number of rows after matching: ", nrow(biom.df[!is.na(biom.df$POS),]))
          biom.df = biom.df[,-which(names(biom.df)=="CHR2")]
          biom.df = biom.df[,-which(names(biom.df)=="POS2")]
      }
     # try splitting the chr:pos?
     if (length(grep("^[0-9]{1,2}[:][1-9][0-9]*$", biom.df$SNPID))>0) split_chr_pos=TRUE else split_chr_pos=FALSE
    if (split_chr_pos) {
       message("Try splitting chr:pos?")
       biom.df$chrpos=paste(biom.df$CHR, biom.df$POS, sep=":")
       biom.df$chrpos=ifelse(biom.df$chrpos=="NA:NA", biom.df$SNPID, biom.df$chrpos)
       biom.df$chrpos = gsub("chr", "", biom.df$chrpos)
       biom.df$CHR=as.numeric(lapply(strsplit(as.character(biom.df$chrpos), ":", fixed=TRUE), "[", 1))
       biom.df$POS=as.numeric(lapply(strsplit(as.character(biom.df$chrpos), ":", fixed=TRUE), "[", 2))
       message("Number of rows after matching: ", nrow(biom.df[!is.na(biom.df$CHR),]))
       message("There are ", nrow(biom.df[is.na(biom.df$CHR),]), " SNPs with no chr pos")
      }
    if (nrow(biom.df[is.na(biom.df$CHR),])>0) {
       biomart=FALSE
       if (biomart) {
       #missing = biom.df[is.na(biom.df$CHR),"SNPID"]
       source("/hpc/users/giambc02/scripts/GENERAL/functions_pipeline.R")
       try= from_rs_to_chrpos(data=biom.df)
       }
   }
   return(biom.df)
   }

addCHRPOS = TRUE
if (addCHRPOS) {
   library(data.table)
   #map_file = "/sc/orga/projects/psychgen/resources/COLOC2/files/GWAS/SNP_MAP_POSITION_HG19.tsv"
   map_file = "/sc/orga/projects/roussp01b/resources/1000G_Phase3_v5_20130502/SNP_MAP_1000G_NO_DUP_POSITION.tsv"
   map = fread(map_file, header=F)
   map = data.frame(map)
   #names(map) = c("CHR", "POS", "SNPID")
   names(map) = c("CHR_POS", "CHR", "POS", "SNPID")
}

if (addCHRPOS) {
   library(data.table)
   map_file2 = "/sc/orga/projects/roussp01b/resources/snp147Common.txt.gz"
   map2 = fread(paste("gunzip -c ", map_file2, sep=""), header=F)
   map2 = data.frame(map2)
   names(map2) = c("CHR", "POS", "SNPID")
}

odir = "/sc/orga/projects/psychgen/resources/COLOC2/files/GWAS5/"

#for (biom.dataset in c("CARDIoGRAMplusC4D", "CARDIoGRAMplusC4D_REC", "MI", "ASD", "UC", "CD", "IBD", "BIP", "AD", "PLAQUE", "CIMT", "GIANT_HIPadjBMI", "GIANT_WC", "GIANT_WCadjBMI", "GIANT_WHR", "GIANT_WHRadjBMI")) { 
for (biom.dataset in c("CARDIoGRAMplusC4D", "CARDIoGRAMplusC4D_REC", "MI", "DAIGRAM", "RA", "ASD", "UC", "CD", "IBD", "SCZ", "BIP", "AD", "PLAQUE", "CIMT", "HDL", "LDL", "TG", "TC", "EXTRAVERSION", "NEUROTICISM", "EduYears", "FEMORAL", "FOREARM", "SPINE", "WellBeing_DS", "WellBeing_Neuroticism", "WellBeing_SWB", "BMI", "HEIGHT", "GIANT_HIP", "GIANT_HIPadjBMI", "GIANT_WC", "GIANT_WCadjBMI", "GIANT_WHR", "GIANT_WHRadjBMI", "MAGIC_Manning_et_al_FastingGlucose_MainEffect", "MAGIC_HbA1C", "MAGIC_Manning_et_al_lnFastingInsulin_MainEffect", "MAGIC_ln_HOMA-B", "MAGIC_ln_HOMA-IR", "URATE")) {
# LIPID TESL

 if (biom.dataset == "ASD") { # autism
     biom.fname ="/sc/orga/projects/roussp01a/resources/GWAS/ASD/pgcasdeuro"
     NCASE = 3303
     NCONTROL = 3428
     N = NCASE + NCONTROL
     typeCC=TRUE
  }

 if (biom.dataset == "UC") { 
     biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/IBD/EUR.UC.gwas.assoc.gz"
     NCASE = 6968
     NCONTROL = 20464
     N = NCASE + NCONTROL
     typeCC=TRUE
  }

 if (biom.dataset == "CD") {
     biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/IBD/EUR.CD.gwas.assoc.gz"
     NCASE = 5956
     NCONTROL = 14927
     N = NCASE + NCONTROL
     typeCC=TRUE
  }


 if (biom.dataset == "IBD") {
     biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/IBD/EUR.IBD.gwas.assoc.gz"
     NCASE = 12882
     NCONTROL = 21770 
     N = NCASE + NCONTROL
     typeCC=TRUE
  }


 if (biom.dataset == "CARDIoGRAMplusC4D") {
     biom.fname="/sc/orga/projects/psychgen/resources/GWAS/CARDIoGRAMplusC4D/rawData/cad.add.160614.website.txt"
     NCASE = 60801
     NCONTROL = 123504
     N = NCASE + NCONTROL
     typeCC=TRUE
  }


 if (biom.dataset == "CARDIoGRAMplusC4D_REC") {
     biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/CAD/cad.rec.090715.web.txt"
     NCASE = 60801
     NCONTROL = 123504
     N = NCASE + NCONTROL
     typeCC=TRUE
  }

 if (biom.dataset == "MI") {
     biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/CAD/mi.add.030315.website.txt"
     NCASE =  43863 # overestimate from the percentage given in Table 1 in Supplementary of individuals who have been genotyped
     NCONTROL = 123504
     N = NCASE + NCONTROL
     typeCC=TRUE
  }

 if (biom.dataset == "SCZ") {
     biom.fname="/sc/orga/projects/roussp01a/resources/PGC2/daner_PGC_SCZ52_0513a.hq2"
     NCASE = 35476
     NCONTROL = 46839
     N = NCASE + NCONTROL
     typeCC=TRUE
  }

 if (biom.dataset == "BIP") {
     biom.fname="/sc/orga/projects/roussp01a/resources/PGC2_BD/daner_PGC_BIP32b_mds7a.gz" 
     NCASE = 20352
     NCONTROL = 31358
     N = NCASE + NCONTROL
     typeCC=TRUE
  }

 if (biom.dataset == "AD") {
     # Use the same one used for fGWAS becasue it is matched to the 1000G for Frequency:
     biom.fname = "/sc/orga/projects/roussp01a/resources/GWAS/AD/IGAP_stage_1_withMAF.txt"
     NCASE = 17008
     NCONTROL = 37154
     N = NCASE + NCONTROL
     typeCC=TRUE
  }

if (biom.dataset == "RA") {
     biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/RA/RA_GWASmeta_European_v2.txt"
     NCASE = 14361 # http://plaza.umin.ac.jp/~yokada/datasource/software.htm
     NCONTROL = 43923
     N = NCASE + NCONTROL
     typeCC=TRUE
  }
# Get SE:
#x$SE= (x$OR_95.CIup - x$OR_95.CIlow)/3.92 # used this
#?? or SE = intervention effect estimate / Z

 if (biom.dataset == "PLAQUE") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/PLAQUE/rawData/Plaque_meta_022016.csv" 
     NCASE = 16697
     NCONTROL = 24448 #total 41145
     N = NCASE + NCONTROL
     typeCC=TRUE
 }

##### QUANTITATIVE TRAITS typeCC=FALSE
## I don't have an example now of typeCC=FALSE, once there is must specify the N instead of the NCASES / NCONTROL
if (biom.dataset == "CIMT") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/CIMT/rawData/IMT.EA.META.MAF1.HetDF4_jun.csv"
     NCASE = NULL
     N = 57384 # but use the TotalSampleSize from the file if have it
     typeCC=FALSE
 }

 if (biom.dataset == "HDL") {
     biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/Lipids/jointGwasMc_HDL.txt"
     NCASE = NULL
     N = 93561
     typeCC=FALSE
 }

 if (biom.dataset == "LDL") {
     biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/Lipids/jointGwasMc_LDL.txt"
     NCASE = NULL
     N=89888
     typeCC=FALSE
 }

 if (biom.dataset == "TG") {
     biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/Lipids/jointGwasMc_TG.txt"
     NCASE = NULL
     N=90785
     typeCC=FALSE
 }

 if (biom.dataset == "TC") {
     biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/Lipids/jointGwasMc_TC.txt"
     NCASE = NULL
     N=94595
     typeCC=FALSE
 }

#### Lipids in Teslovich data # Reformat this or the RData file (b/c have N samples size in here) /sc/orga/projects/epigenAD/coloc/data/lipids_Teslovich/rawData/all.txt
# /sc/orga/projects/epigenAD/coloc/data/lipids_Teslovich/summaryStats
 if (biom.dataset == "HDL_tesl") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/lipids_Teslovich/rawData/hdl_tesl_with_Effect_with_chr_pos.txt"
     NCASE = NULL
     N = 96857
     typeCC=FALSE
 }

 if (biom.dataset == "LDL_tesl") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/lipids_Teslovich/rawData/ldl_tesl_with_Effect_with_chr_pos.txt"
     NCASE = NULL
     N=92454
     typeCC=FALSE
 }

 if (biom.dataset == "TG_tesl") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/lipids_Teslovich/rawData/logtg_tesl_with_Effect_with_chr_pos.txt"
     NCASE = NULL
     N=93512
     typeCC=FALSE
 }

 if (biom.dataset == "TC_tesl") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/lipids_Teslovich/rawData/tc_tesl_with_Effect_with_chr_pos.txt"
     NCASE = NULL
     N=98653
     typeCC=FALSE
 }

 if (biom.dataset == "EXTRAVERSION") {
     #biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/Personality/GPC-2.EXTRAVERSION.full.txt"
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/EXTRAVERSION/rawData/GPC-2.EXTRAVERSION.full.txt"
     NCASE = NULL
     N = 63030
     typeCC=FALSE
  }  

 if (biom.dataset == "NEUROTICISM") {
     #biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/Personality/GPC-2.NEUROTICISM.full.txt"
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/NEUROTICISM/rawData/GPC-2.NEUROTICISM.full.txt"
     NCASE = NULL
     N = 63661
     typeCC=FALSE
  }

 if (biom.dataset == "EduYears") {
     #biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/EducAt/EduYears_Main.txt.gz"
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/EduYears/rawData/EduYears_Main.txt.gz"
     NCASE = NULL
     N = 293723
     typeCC=FALSE
  }


 if (biom.dataset == "FEMORAL") {
     #biom.fname= "/sc/orga/projects/roussp01a/resources/GWAS/Osteoporosis/fn2stu.MAF0_.005.pos_.out_.gz"
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/FEMORAL/rawData/fn2stu.MAF0_.005.pos_.out_"
     NCASE = NULL
     N = 32735
     typeCC=FALSE
  }

 if (biom.dataset == "FOREARM") {
     #biom.fname= "/sc/orga/projects/roussp01a/resources/GWAS/Osteoporosis/fa2stu.MAF0_.005.pos_.out_.gz" 
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/FOREARM/rawData/fa2stu.MAF0_.005.pos_.out_"
     NCASE = NULL
     N = 8143
     typeCC=FALSE
  }

 if (biom.dataset == "SPINE") {
     #biom.fname= "/sc/orga/projects/roussp01a/resources/GWAS/Osteoporosis/ls2stu.MAF0_.005.pos_.out_.gz"
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/SPINE/rawData/ls2stu.MAF0_.005.pos_.out_"
     NCASE = NULL
     N = 28498
     typeCC=FALSE
  }


 if (biom.dataset == "WellBeing_DS") {
     biom.fname= "/sc/orga/projects/roussp01a/resources/GWAS/WellBeing/DS_Full.txt.gz"
     NCASE = NULL
     N = 161460
     typeCC=FALSE
  }

 if (biom.dataset == "WellBeing_Neuroticism") {
     biom.fname="/sc/orga/projects/roussp01a/resources/GWAS/WellBeing/Neuroticism_Full.txt.gz"
     NCASE = NULL
     N = 170911
     typeCC=FALSE
  }

 if (biom.dataset == "WellBeing_SWB") {
     biom.fname= "/sc/orga/projects/roussp01a/resources/GWAS/WellBeing/SWB_Full.txt.gz"
     NCASE = NULL
     N = 298420
     typeCC=FALSE
  }


##########
# DAIGRAM # T2DS
 if (biom.dataset == "DIAGRAM") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/DIAGRAM/rawData/DIAGRAMv3.2012DEC17.txt"
     NCASE = 11645
     N=44414
     typeCC=TRUE
 }

########## GIANT (BMI) hapmap
 if (biom.dataset == "BMI") {
     # European ancestry only
     # ref Locke 2015 Nature https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4382211/
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/BMI/rawData/SNP_gwas_mc_merge_nogc.tbl.uniq"
     NCASE = NULL
     #N=322079
     N = NULL
     typeCC=FALSE
 }

 if (biom.dataset == "HEIGHT") {
     # ref Wood 2014 Nature Genetics http://www.ncbi.nlm.nih.gov/pubmed/25282103
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/HEIGHT/rawData/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt"
     NCASE = NULL
     N=NULL
     typeCC=FALSE
 }


 if (biom.dataset == "GIANT_HIP") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/GIANT_HIP/rawData/GIANT_2015_HIP_COMBINED_EUR.txt"
     NCASE = NULL
     N=NULL
     typeCC=FALSE
 }

 if (biom.dataset == "GIANT_HIPadjBMI") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/GIANT_HIP/rawData/GIANT_2015_HIPadjBMI_COMBINED_EUR.txt"
     NCASE = NULL
     N=NULL
     typeCC=FALSE
 }

 if (biom.dataset == "GIANT_WC") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/GIANT_WC/rawData/GIANT_2015_WC_COMBINED_EUR.txt"
     NCASE = NULL
     N=NULL
     typeCC=FALSE
 }

 if (biom.dataset == "GIANT_WCadjBMI") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/GIANT_WC/rawData/GIANT_2015_WCadjBMI_COMBINED_EUR.txt"
     NCASE = NULL
     N=NULL
     typeCC=FALSE
 }

 if (biom.dataset == "GIANT_WHR") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/GIANT_WHR/rawData/GIANT_2015_WHR_COMBINED_EUR.txt"
     NCASE = NULL
     N=NULL
     typeCC=FALSE
 }

 if (biom.dataset == "GIANT_WHRadjBMI") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/GIANT_WHR/rawData/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt"
     NCASE = NULL
     N=NULL
     typeCC=FALSE
 }

##########
# /sc/orga/projects/epigenAD/coloc/data/MAGIC_README
#include: MAGIC_Manning_et_al_FastingGlucose_MainEffect, MAGIC_HbA1C, MAGIC_Manning_et_al_lnFastingInsulin_MainEffect, MAGIC_ln_HOMA-B, MAGIC_ln_HOMA-IR

 if (biom.dataset == "MAGIC_Manning_et_al_FastingGlucose_MainEffect") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/MAGIC_Manning_et_al_FastingGlucose_MainEffect/rawData/MAGIC_Manning_et_al_FastingGlucose_MainEffect.txt"
    # Not adjusted for BMI
     NCASE = NULL
     N=58074
     typeCC=FALSE
 }

 if (biom.dataset == "MAGIC_HbA1C") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/MAGIC_HbA1C/rawData/MAGIC_HbA1C.txt"
     NCASE = NULL
     N=46368
     typeCC=FALSE
 }

 if (biom.dataset == "MAGIC_Manning_et_al_lnFastingInsulin_MainEffect") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/MAGIC_Manning_et_al_lnFastingInsulin_MainEffect/rawData/MAGIC_Manning_et_al_lnFastingInsulin_MainEffect.txt"
    # Not adjusted for BMI
     NCASE = NULL
     N=51750
     typeCC=FALSE
 }

 if (biom.dataset == "MAGIC_ln_HOMA-B") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/MAGIC_ln_HOMA-B/rawData/MAGIC_ln_HOMA-B.txt"
     NCASE = NULL
     N=46186
     typeCC=FALSE
 }

 if (biom.dataset == "MAGIC_ln_HOMA-IR") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/MAGIC_ln_HOMA-IR/rawData/MAGIC_ln_HOMA-IR.txt"
     NCASE = NULL
     N=46186
     typeCC=FALSE
 }

## not included for now:
 if (biom.dataset == "MAGIC_2hrGlucose_AdjustedForBMI") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/MAGIC_2hrGlucose_AdjustedForBMI/rawData/MAGIC_2hrGlucose_AdjustedForBMI.txt"
     NCASE = NULL
     N=15234
     typeCC=FALSE
 }

 if (biom.dataset == "MAGIC_FastingGlucose") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/MAGIC_FastingGlucose/rawData/MAGIC_FastingGlucose.txt"
     NCASE = NULL
     N=46186
     typeCC=FALSE
 }

 if (biom.dataset == "MAGIC_ln_FastingInsulin") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/MAGIC_ln_FastingInsulin/rawData/MAGIC_ln_FastingInsulin.txt"
     NCASE = NULL
     N=46186
     typeCC=FALSE
 }


 if (biom.dataset == "MAGIC_ln_fastingProinsulin") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/MAGIC_ln_fastingProinsulin/rawData/MAGIC_ln_fastingProinsulin.txt.gz"
     NCASE = NULL
     N=10701
     typeCC=FALSE
 }

########## This includes other nine traits for glucose-stimulated insulin secretion (GSIS) indices during an oral glucose tolerance test
 if (biom.dataset == "MAGIC_InsulinSecretion_data_release_May2014") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/MAGIC_InsulinSecretion_data_release_May2014/rawData/MAGIC_InsulinSecretion_data_release_May2014.zip"
     NCASE = NULL
     N=5318
     typeCC=FALSE
 }


 if (biom.dataset == "URATE") {
     biom.fname="/sc/orga/projects/epigenAD/coloc/data/GIANT_WHR/rawData/GUGC_MetaAnalysis_Results_UA.csv"
     NCASE = NULL
     N=NULL
     typeCC=FALSE
 }



############################################################
outFile= paste(odir, biom.dataset, "_formatted", sep="")

cat("\n-----------------------------------------------------------------\nSTUDY:",biom.dataset,"file:",biom.fname,"output:", outFile,"\n-----------------------------------------------------------------\n")

#message("STUDY: ", biom.dataset)
if (typeCC) type="cc" else type="quant"

source("/hpc/users/giambc02/scripts/COLOC/pipeline/functions_coloc_pipeline.R")
      #OUT_DIR= "/sc/orga/projects/psychgen/resources/COLOC2/files/GWAS3/"
      #dir.create(file.path(OUT_DIR), showWarnings = FALSE)

if (biom.dataset %in% c("FEMORAL", "FOREARM", "SPINE")) {
   #biom.dataset = "FEMORAL"
   #fname = "/sc/orga/projects/epigenAD/coloc/data/FEMORAL/rawData/fn2stu.MAF0_.005.pos_.out_" 
   #N = 32735
   #biom.dataset = "FOREARM"
   #biom.fname="/sc/orga/projects/epigenAD/coloc/data/FOREARM/rawData/fa2stu.MAF0_.005.pos_.out_"
   #N = 8143
   #biom.dataset = "SPINE"
   #biom.fname="/sc/orga/projects/epigenAD/coloc/data/SPINE/rawData/ls2stu.MAF0_.005.pos_.out_"
   #N = 28498
   outFile= paste(odir, biom.dataset, "_formatted", sep="")
   biom.df= fread(biom.fname, header=T)
   biom.df= data.frame(biom.df)
   biom.df=biom.df[,c("rs_number", "chromosome", "position", "eaf", "beta", "se", "p.value", "reference_allele", "other_allele")]
   names(biom.df) = c("SNPID", "CHR", "POS", "F", "BETA", "SE", "PVAL", "A1", "A2")
   biom.df$N= N
   biom.df$SNPID= gsub("chr", "", biom.df$SNPID)
   biom.df$CHR= gsub("chr", "", biom.df$CHR)
   n_cutoff = quantile(biom.df$N, c(.50))
   info = data.frame(data=biom.dataset, NsnpsCompleteData_before_filter = nrow(biom.df[complete.cases(biom.df),]), PVAL_range = paste(range(biom.df$PVAL), collapse=" , "), NsamplesMax= max(biom.df$N), NsamplesMin=min(biom.df$N), Ncut=as.numeric(n_cutoff), filtered_snps=nrow(biom.df)-nrow(biom.df_cut))

   write.table(info, file=paste(odir, "_info.txt", sep=""), append = TRUE, sep="\t")

   write.table(biom.df, file=outFile, row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
   write.table(biom.df_cut, file=paste(odir, biom.dataset, "_formatted_Nfiltered", sep=""), row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
   message("File for ", biom.dataset, " saved in ", outFile)
   next()
} 

if (biom.dataset %in% c("DIAGRAM", "RA")) {
   outFile= paste(odir, biom.dataset, "_formatted", sep="")
   biom.df= fread(biom.fname, header=T)
   biom.df= data.frame(biom.df)
   if (biom.dataset=="DIAGRAM") {
     biom.df$N=biom.df$N_CASES+biom.df$N_CONTROLS
     # find SE from CI
     biom.df$SE= (biom.df$OR_95U - biom.df$OR_95L)/3.92 # used this
     #?? or SE = intervention effect estimate / Z
     biom.df=biom.df[,c("SNP", "CHROMOSOME", "POSITION", "P_VALUE", "OR", "SE", "N_CASES", "N", "RISK_ALLELE", "OTHER_ALLELE")]
   } else {
     biom.df$Ncases = NCASE
     biom.df$N = N
     biom.df$SE= (biom.df$OR_95.CIup - biom.df$OR_95.CIlow)/3.92 # used this
     biom.df = biom.df[,c("SNPID", "Chr", "Position.hg19.", "P.val", "OR.A1.", "SE", "Ncases", "N", "A1", "A2")]
   }
   names(biom.df)=c("SNPID", "CHR", "POS", "PVAL", "BETA", "SE", "Ncases", "N", "A1", "A2")
   biom.df$SNPID= gsub("chr", "", biom.df$SNPID)
   biom.df$CHR= gsub("chr", "", biom.df$CHR)
   n_cutoff = quantile(biom.df$N, c(.50))
   biom.df_cut = biom.df[biom.df$N >=n_cutoff,]
   info = data.frame(data=biom.dataset, NsnpsCompleteData_before_filter = nrow(biom.df[complete.cases(biom.df),]), PVAL_range = paste(range(biom.df$PVAL), collapse=" , "), NsamplesMax= max(biom.df$N), NsamplesMin=min(biom.df$N), Ncut=as.numeric(n_cutoff), filtered_snps=nrow(biom.df)-nrow(biom.df_cut))

   write.table(info, file=paste(odir, "_info.txt", sep=""), append = TRUE, sep="\t")

   write.table(biom.df, file=outFile, row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
   write.table(biom.df_cut, file=paste(odir, biom.dataset, "_formatted_Nfiltered", sep=""), row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
   message("File for ", biom.dataset, " saved in ", outFile)
   next()
}

if (biom.dataset %in% c("MAGIC_Manning_et_al_FastingGlucose_MainEffect", "MAGIC_Manning_et_al_lnFastingInsulin_MainEffect")) {
   biom.df= fread(biom.fname, header=T)
   biom.df= data.frame(biom.df)
   biom.df$N = N
   n_cutoff = quantile(biom.df$N, c(.50))
   biom.df_main =  biom.df[,c("Snp", "maf", "MainEffects", "MainSE", "MainP", "N", "effect_allele", "other_allele")]
   biom.df_BMIAdjusted = biom.df[,c("Snp", "maf", "BMIadjMainEffects", "BMIadjMainSE", "BMIadjMainP", "N", "effect_allele", "other_allele")]
   names(biom.df_main) = c("SNPID", "F", "BETA", "SE", "PVAL", "N", "A1", "A2")
   names(biom.df_BMIAdjusted) = c("SNPID", "F", "BETA", "SE", "PVAL", "N", "A1", "A2")
   t = addCHRPOS_fn(biom.df_main)
   outFile= paste(odir, biom.dataset, "_formatted", sep="")

   write.table(t, file=outFile, row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
   write.table(t[t$N >=n_cutoff,], file=paste(odir, biom.dataset, "_formatted_Nfiltered", sep=""), row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")

   t = addCHRPOS_fn(biom.df_BMIAdjusted)
   outFile= paste(odir, biom.dataset, "_BMIAdjusted_formatted", sep="")

   write.table(t, file=outFile, row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
   write.table(t[t$N >=n_cutoff,], file=paste(odir, biom.dataset, "_BMIAdjusted_formatted_Nfiltered", sep=""), row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
   next()
}

if (biom.dataset=="URATE") {
   biom.df= fread(biom.fname, header=T)
   biom.df= data.frame(biom.df)
   names(biom.df)=c("SNPID", "N", "A1", "A2", "BETA", "SE", "PVAL")
   biom.df = addCHRPOS_fn(biom.df)
   n_cutoff = quantile(biom.df$N, c(.50))
   biom.df_cut = biom.df[biom.df$N >=n_cutoff,]
   info = data.frame(data=biom.dataset, NsnpsCompleteData_before_filter = nrow(biom.df[complete.cases(biom.df),]), PVAL_range = paste(range(biom.df$PVAL), collapse=" , "), NsamplesMax= max(biom.df$N), NsamplesMin=min(biom.df$N), Ncut=as.numeric(n_cutoff), filtered_snps=nrow(biom.df)-nrow(biom.df_cut))

   write.table(info, file=paste(odir, "_info.txt", sep=""), append = TRUE, sep="\t")

   write.table(biom.df, file=outFile, row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
   write.table(biom.df_cut, file=paste(odir, biom.dataset, "_formatted_Nfiltered", sep=""), row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
   message("File for ", biom.dataset, " saved in ", outFile)
   next()
}

biom.df = tryCatch(formatColoc(fname = biom.fname, type=type, N=N, Ncases=NCASE, info_filter=0, maf_filter=0, fread=T, eqtl=FALSE), error=function(e) NULL )
      if (is.null(biom.df)) {
          message("******************************************!!!!! Study could not be processed")
          write.table(paste("******************************************!!!!! Study ", biom.dataset, " could not be processed", sep=""), file=paste(odir, "_info.txt", sep=""), append = TRUE, sep="\t")
          next()
       }

if (!("CHR" %in% names(biom.df)) | !("POS" %in% names(biom.df))) addCHRPOS=TRUE else addCHRPOS=FALSE
# take out "chr" from SNPID
biom.df$SNPID= gsub("chr", "", biom.df$SNPID)

if (addCHRPOS) {
   biom.df  = addCHRPOS_fn(biom.df)
}

# take out "chr" from CHR
biom.df$CHR= gsub("chr", "", biom.df$CHR)

# print some info for the study
message("Number of SNPs: ", nrow(biom.df))
message("Number of SNPs with complete data: ", nrow(biom.df[complete.cases(biom.df),]))
message("Range for PVALs ", paste(range(biom.df$PVAL), collapse=" , "))

message("Max number of samples: ", max(biom.df$N))
message("Min number of samples: ", min(biom.df$N))

n_cutoff = quantile(biom.df$N, c(.50)) 
message("The cut off for the 50th centile is: ", n_cutoff)
biom.df_cut = biom.df[biom.df$N >=n_cutoff,]
message("Would remove: ", nrow(biom.df)-nrow(biom.df_cut), " SNPs")


info = data.frame(data=biom.dataset, NsnpsCompleteData_before_filter = nrow(biom.df[complete.cases(biom.df),]), PVAL_range = paste(range(biom.df$PVAL), collapse=" , "), NsamplesMax= max(biom.df$N), NsamplesMin=min(biom.df$N), Ncut=as.numeric(n_cutoff), filtered_snps=nrow(biom.df)-nrow(biom.df_cut))
 
write.table(info, file=paste(odir, "_info.txt", sep=""), append = TRUE, sep="\t")

     write.table(biom.df, file=outFile, row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
     write.table(biom.df_cut, file=paste(odir, biom.dataset, "_formatted_Nfiltered", sep=""), row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")
     message("File for ", biom.dataset, " saved in ", outFile)


}
