sniff <- function(file="/sc/orga/projects/epigenAD/coloc/CARDIoGRAMplusC4D/rawData/cad.add.160614.website.txt", eqtl = FALSE) {
  # check if there is a header
  if (all(sapply(read.table(file, header = F, nrow = 1, stringsAsFactors=F), class)=="character")) header=TRUE else header=FALSE
  if (!header) {
     checkIflogical=sapply(read.table(file, header = F, nrow = 1, stringsAsFactors=F), class)
     if (length(checkIflogical[checkIflogical=="logical"])>0) checkIflogical[checkIflogical=="logical"]="character"
     if (all(checkIflogical=="character")) header=TRUE else header=FALSE
  }
  if (!header) stop("There is no header in the data")

      line = read.table(file, header = T, nrow = 100, stringsAsFactors=F)
      if (ncol(line)==1) {
          separator=","
          line = read.table(file, header = T, nrow = 100, stringsAsFactors=F, sep=separator)
      }
      SNPID = which(grepl("snp|SNP|MarkerName", names(line),  ignore.case = T))
      CHR = which(grepl("^chr$|^CHR$|^Chromosome$|^chr_name$", names(line),  ignore.case = T))
      POS = which(grepl("^pos$|^bp$|^bp_hg19$|^position|^chrom_start$", names(line),  ignore.case = T))
      BETA = which(grepl("^beta$|^effect$|^or$|OR.A1|BMIadjMainEffects", names(line),  ignore.case = T))
      F = which(grepl("^F$|freq|FRQ|MAF", names(line),  ignore.case = T))[1] # sometimes have 2 (one for each allele), doesn't matter whcih you take for our applications (fGWAS and coloc)
          if (is.na(F)) F=integer(0)
      PVAL = which(grepl("^p$|pval|Pvalue|P.value|P.val|BMIadjMainP", names(line),  ignore.case = T))
      SE = which(grepl("^se|^StdErr$|BMIadjMainSE", names(line),  ignore.case = T))
      if (length(PVAL)>0) {
         #Make sure this is not a heterogenous pvalue, if it is set to 0
         if (any(grepl("het", names(line)[PVAL], ignore.case =T))) {
           message("Pvalue found is a heterogenous p-value. Looking for another pvalue")
           PVAL = PVAL[!grepl("het", names(line)[PVAL], ignore.case =T)]
           }
         if (length(PVAL)==0) {
           PVAL = which(grepl("^p", names(line),  ignore.case = T))
           #print(names(line)[pval])
          }
      }

      if (eqtl) {
         ProbeID = which(grepl("ProbeID", names(line),  ignore.case = T))
         if (length(ProbeID)==0) stop("Column ProbeID is missing")
      }
      # Necessary columns
       if (all(c(length(CHR), length(POS))==0) & length(SNPID)!=0) print("No chr, pos: must find from rsid")
             if (any(c(length(SNPID), length(PVAL))==0)) stop("Column names ", c("snp", "pval")[c(length(snp), length(pval))==0], " not recognized")
      if (any(c(length(SNPID), length(CHR), length(POS), length(BETA), length(SE), length(PVAL))>1)) stop("Some values match more than one column name: change complicated column names")

      # If have a sample size for each SNP, output this
      N = which(grepl("^N$|^TotalSampleSize$|^SAMPLE_SIZE$", names(line),  ignore.case = T))
      Ncases =  which(grepl("^N_CASES$|^N_CASE$", names(line),  ignore.case = T))
      Ncontrols =  which(grepl("^N_CONTROLS$|^N_CONTROL$", names(line),  ignore.case = T))

      info = which(grepl("INFO|RSQ", names(line),  ignore.case = T))
      Z = which(grepl("^Z$|zscore", names(line),  ignore.case = T))
      ### Sanity checks
      if (length(PVAL)>0) if ( sum(line[,PVAL] <0 | line[,PVAL]>1) >0 ) message("Invalid Pvalues")
      if (length(SE)>0) if ( sum(line[,SE] <=0 | line[,SE]=="Inf" | line[,SE]>10) >0 ) message("Invalid SE")
      if (length(info)>0) if ( sum(line[,info] <0 | line[,info]>1) >0 ) message("Imputation information is less than 0 or greater than 1")

      print(line[1:2,])
      message("SNP column: ", names(line)[SNPID])
      message("CHR column: ", names(line)[CHR])
      message("POS column: ", names(line)[POS])
      message("EFFECT column: ", names(line)[BETA])
      message("PVAL column: ", names(line)[PVAL])
      message("Nsnp column: ", names(line)[N])
      message("INFO column: ", names(line)[info])
      message("SE column: ", names(line)[SE])
      message("Z column: ", names(line)[Z])
      message("FREQ column: ", names(line)[F])
      message("Ncases column: ", names(line)[Ncases])
      message("Ncontrols column: ", names(line)[Ncontrols])
if (eqtl) message("ProbeID column: ", names(line)[ProbeID])

      output_cols = c(SNPID=SNPID, CHR=CHR, POS=POS, BETA=BETA, PVAL=PVAL, N=N, info=info, SE=SE, Z=Z, F=F, Ncases=Ncases, Ncontrols=Ncontrols)
      if (eqtl) output_cols = c(SNPID=SNPID, CHR=CHR, POS=POS, BETA=BETA, PVAL=PVAL, N=N, info=info, SE=SE, Z=Z, F=F, Ncases=Ncases, Ncontrols=Ncontrols, ProbeID=ProbeID)
      # Also output a list of the classes
      output_class = c(class(line[,SNPID]), class(line[,CHR]),class(line[,POS]),class(line[,BETA]),class(line[,PVAL]),class(line[,N]),class(line[,info]),class(line[,SE]),class(line[,Z]),class(line[,F]), class(line[,Ncases]), class(line[,Ncontrols]))
      if (eqtl) output_class = c(class(line[,SNPID]), class(line[,CHR]),class(line[,POS]),class(line[,BETA]),class(line[,PVAL]),class(line[,N]),class(line[,info]),class(line[,SE]),class(line[,Z]),class(line[,F]), class(line[,Ncases]), class(line[,Ncontrols]), class(line[,ProbeID]))
      # if data.frame make it NA for missing columns
      output_class[which(output_class=="data.frame")]=rep("NA", length(which(output_class=="data.frame")))
      # remove missing
      output_class = output_class[output_class!="NA"]
      # if integer make it numeric
      output_class[which(output_class=="integer")]=rep("numeric", length(which(output_class=="integer")))

      output=list(output_cols, output_class)
      return(output)
   }


formatColoc <- function(fname = "/Users/claudiagiambartolomei/Dropbox/MEGA/coloc/example/eQTL.GENE_FDR_0p05.db", type="cc", N=NA, Ncases=NA, info_filter=0.3, maf_filter=0.001, eqtl=FALSE) {
  require(data.table)
  colsAll = sniff(file=fname, eqtl=eqtl)
  # fread select reads the columns in order so must re-order
  colsAll[[1]] = colsAll[[1]][order(colsAll[[1]])]
  data= fread(fname, select = colsAll[[1]], col.names=names(colsAll[[1]]))
  data = as.data.frame(data)
 # what to do about the sample size?
  if (type=="quant") {
  if (!("N" %in% names(colsAll[[1]]))) {
  if (is.na(N)) stop("Please specify a sample size because it is not found within the data")
     data$N = N
   }
  }
  if (type=="cc") {
  if (!all((c("Ncases", "Ncontrols") %in% names(colsAll[[1]])))) {
  if ( is.na(Ncases) ) stop("The dataset is a case-control so must specify number of cases")
     data$Ncases = Ncases
     data$N = N
    }
  }

  if (length(which(names(data)=="info"))>0) {
    data = data[data$info > info_filter,]
  } else {
  message("Cannot find info so cannot filter")
  }

  if (length(which(names(data)=="F"))>0) {
    # could be either MAF or a frequency, but to filter find MAF
    data$MAF = ifelse(data$F<0.5, data$F, 1-data$F)
    data = data[data$MAF > maf_filter,]
   } else {
   message("Cannot find frequency so cannot filter")
  }

  # you need either effect and SE or pval for coloc
  if ( all(!(c("BETA", "SE") %in% names(data))) | !("PVAL" %in% names(data)) ) stop("Need either BETA and SE, or PVAL in the data")
  return(data)
}



