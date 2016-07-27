#########################################################
# This function takes in:
# nameList -	Character vector with the names of the first layer of the datasets (e.g. Enhancers, TF, BroadPeaks...)
 ## If nameList not specified the function will take all the annotations in an RData file
# RDataFile -	.RData file path with the annotation info 
 ## Can take the complete .RData file for REMC, but keep in mind
 ## system.time(overlapGR(RDataFile="/sc/orga/projects/roussp01a/Claudia_TMP/data/REMC.BroadPeaks.Rdata", snpDf=region)) # This takes 4.6 mins
 ## system.time(overlapGR(nameList="BroadPeaks", RDataFile="/sc/orga/projects/roussp01a/Claudia_TMP/data/REMC.Rdata", snpDf=region)) # This takes 11.51 mins!
# snpDf  -	Data frame with column names "CHR", "START", "STOP", "SNP", already loaded in R
 ## ex.
 ## t=overlapGR(RDataFile="/sc/orga/projects/roussp01a/Claudia_TMP/data/SE.Rdata", snpDf=region)
 ## t1=overlapGR(nameList=c("coreHMM","BroadPeaks"), RDataFile="/sc/orga/projects/roussp01a/Claudia_TMP/data/REMC.Rdata", snpDf=region)
## For REMC (auxHMM, coreHMM) and ENCODE (chromHMM) SNPs are overlapped with each row of the Genomic Range object

overlapGR <- function(nameList=NA, RDataFile, snpDf) {
         suppressPackageStartupMessages({
         require(GenomicRanges);
         });
         #options(warn=2) # treat warnings as errors
         
         RDataName=gsub(".Rdata", "", basename(RDataFile))
         message("Loading ", RDataFile, " ...")
         listANNOT=get(load(RDataFile))
 
         # Add "chr" in seqlevels if don't have (because annotList has it)
         addChr=ifelse(any(grep("chr",snpDf$CHR))>0, FALSE,TRUE)
         if (addChr) (snpDf$CHR=paste("chr", snpDf$CHR, sep=""))

         snps=GRanges(seqnames=snpDf$CHR,
               ranges=IRanges(start=snpDf$START,
                              end=snpDf$STOP, names=snpDf$SNP),
               name=snpDf$SNP) # names of regions are the index SNPs

         ANNOT=c()
         c.names=c()

         ##if (!is.na(nameList) & all(grepl("$", nameList, fixed=T))) {
                      # subset listANNOT (for now can only subset based on the elements of one name
         ##             subnameList=names(table(as.character(lapply(strsplit(as.character(nameList), "$", fixed=TRUE), "[", 2))))
         ##             nameList=names(table(as.character(lapply(strsplit(as.character(nameList), "$", fixed=TRUE), "[", 1))))
         ##             outlist = list()
                      #for (name in nameAll) {
         ##              listANNOT = sapply(listANNOT, "[", which(names(listANNOT[[nameList]])%in%subnameList))                   
                       #outlist[[length(outlist)+1]] <- listANNOTsubset
                      #}
                      #listANNOT = outlist
         ##}

         if (all(is.na(nameList))) (nameList=names(listANNOT))

         for (name in nameList) {
         message("Making a binary overlap matrix for data ", name)
            if (!is.list(listANNOT[[name]]) & !(any(grepl("chromHMM|auxHMM|coreHMM", name)==TRUE) | any(grepl("chromHMM|auxHMM|coreHMM", RDataName)==TRUE))) { # This is for SE for ex.
             if (class(listANNOT[[name]])=="GRanges") { # This is for ATACSeq peaks
               OBJ <- countOverlaps(snps, listANNOT[[name]], type=c("within"))
               ANNOT=cbind(ANNOT, OBJ)
               c.names = c(c.names, paste(RDataName, name, colnames(OBJ), sep="."))
              } else {
               OBJ <- do.call(cbind, lapply(listANNOT[[name]], function(x) countOverlaps(snps, x, type=c("within"))), quote=TRUE)
               ANNOT=cbind(ANNOT, OBJ)
               c.names = c(c.names, paste(RDataName, name, colnames(OBJ), sep="."))
              } 
             }
            if (is.list(listANNOT[[name]]) &  !(any(grepl("chromHMM|auxHMM|coreHMM", name)==TRUE) | any(grepl("chromHMM|auxHMM|coreHMM", RDataName)==TRUE))) { # This is the case for REMC Peaks
              message(name, " contains additional list layers") 
              for (i in names(listANNOT[[name]])) {
                   if (!is.list(listANNOT[[name]][[i]])) {
                    OBJ <- do.call(cbind, lapply(listANNOT[[name]][[i]], function(x) countOverlaps(snps, x, type=c("within"))), quote=TRUE)
                    ANNOT=cbind(ANNOT, OBJ)
                    c.names = c(c.names, paste(RDataName, name, i, colnames(OBJ), sep="."))
                    }
                   }
               }
            if (any(grepl("chromHMM|auxHMM|coreHMM", name)==TRUE) | any(grepl("chromHMM|auxHMM|coreHMM", RDataName)==TRUE)) { 
              if (any(grepl("auxHMM|coreHMM", RDataName)==TRUE)) {
                list_state = split(listANNOT[[name]], listANNOT[[name]]$id)
                OBJ <- do.call(cbind, lapply(list_state, function(x) countOverlaps(snps, x, type=c("within"))), quote=TRUE)
                ANNOT=cbind(ANNOT, OBJ)
                c.names = c(c.names, paste(RDataName, name, colnames(OBJ), sep="."))
                }
               if (RDataName=="REMC" & any(grepl("coreHMM", name)==TRUE) | (RDataName=="ENCODE" & any(grepl("chromHMM", nameList)==TRUE))) {
                 for (i in names(listANNOT[[name]])) {
                  list_state = split(listANNOT[[name]][[i]], listANNOT[[name]][[i]]$id)
                  OBJ <- do.call(cbind, lapply(list_state, function(x) countOverlaps(snps, x, type=c("within"))), quote=TRUE)
                  ANNOT=cbind(ANNOT, OBJ)
                  c.names = c(c.names, paste(RDataName, name, i, colnames(OBJ), sep="."))
                  }
                }
            }
          colnames(ANNOT) = c.names
          #if(max(ANNOT)>1) message("Matrix is not binary") # Meaning the SNP overlaps several regions of this annotation
          ANNOT[ANNOT >0] <- 1

          }
          return(ANNOT)
          options(warn=1)  # restore default for warnings, i.e. print warnings as they occur
          }


# This function makes snpDf providing a character vector of SNP IDs with format chr:pos
MAKEsnpDf <- function(snpid) {
        # check that it is in the correct format
        if (length(grep("[[:digit:]]+:[[:digit:]]", snpid)) != length(snpid)) stop("SNP name is not in correct format")
        addChr=ifelse(any(grep("chr",snpid))>0, FALSE,TRUE)
        chr=lapply(strsplit(as.character(snpid), ":", fixed=TRUE), "[", 1)
        if (addChr) (chr=paste("chr", chr, sep="")) else (chr=as.character(chr))
        snpDf=data.frame(CHR=chr, 
        START=as.numeric(lapply(strsplit(as.character(snpid), ":", fixed=TRUE), "[", 2)),
        STOP=as.numeric(lapply(strsplit(as.character(snpid), ":", fixed=TRUE), "[", 2))+1, 
        SNP=snpid)
        return(snpDf)
        }


# FOR SNPs:        
# if it is 1 base only, need to add 1 base pair otherwise overlaps won't work
GRtoBED <- function(gr=gr) {
  df <- data.frame(seqnames=seqnames(gr),
   starts=start(gr)-1,
   ends=end(gr),
   id = gr$id)
#   names=c(rep(".", length(gr))),
#   scores=c(rep(".", length(gr))),
#   strands=strand(gr))

#write.table(df, file="foo.bed", quote=F, sep="\t", row.names=F, col.names=F)
   return(df)
   }

GRListtoBED <- function(grList) {
   BEDList= lapply(grList, GRtoBED)
   return(BEDList)
   }

# FOR INTERVALS
# if it is 1 base only, need to add 1 base pair otherwise overlaps won't work
GRtoBEDforIntervals <- function(gr=gr) {
  df <- data.frame(seqnames=seqnames(gr),
   starts=start(gr),
   ends=end(gr),
   id = gr$id)
   return(df)
   }

GRListtoBEDforIntervals <- function(grList) {
   BEDList= lapply(grList, GRtoBEDforIntervals)
   return(BEDList)
   }

# Function to read files and make Granges while using IDs directly from the BED files
# For REMC (auxHMM, coreHMM); ENCODE (chromHMM) we want to keep the same IDs as the BED file
BEDtoGR = function(file_path) {
  dat0 <- read.table(gzfile(as.character(file_path)), header = FALSE, colClasses=c(NA, NA, NA, NA, "NULL", "NULL"))
  colnames(dat0) = c("chr", "start", "end", "id")
  bed0 <- with(dat0, GRanges(chr, IRanges(start, end), id=id))
  return(bed0)
}

# Function to make GrangesList from a set of BED files
# e.g. GRList= BEDtoGRList(list.files("/sc/orga/projects/roussp01a/Claudia_TMP/ATACseq_MULTIREGIONAL/", full.names=T))
BEDtoGRList = function(list.files) {
   BEDList= lapply(list.files, BEDtoGR)
   names(BEDList) = basename(list.files)
   return(BEDList)
   }
   

# Function to apply to all combinations of two lists and produce a list
comb_apply <- function(f,..., MoreArgs=list()){
      exp <- unname(as.list(expand.grid(...,stringsAsFactors = FALSE)))
      do.call(mapply, c(list(FUN=f, SIMPLIFY=FALSE, MoreArgs=MoreArgs), exp))
}


overlapGRwithListDF <- function(nameList=NA, RDataFile, snpDfList) {
         suppressPackageStartupMessages({
         require(GenomicRanges);
         });

         RDataName=gsub(".Rdata", "", basename(RDataFile))
         message("Loading ", RDataFile, " ...")
         listANNOT=GRangesList(get(load(RDataFile)))

         snpsGRList = GRangesList(lapply(snpDfList,makeGRangesFromDataFrame))
         ANNOT=comb_apply(function(i, j) countOverlaps(snpsGRList[[i]], listANNOT[[j]]), seq_along(snpsGRList), seq_along(listANNOT))
         names(ANNOT)=apply(expand.grid(names(snpsGRList), names(listANNOT)), 1, paste, collapse="_")

         d <- do.call(rbind, ANNOT)
         rownames(d)=as.character(lapply(strsplit(as.character(rownames(d)), "_", fixed=TRUE), "[", 3))
         d=as.data.frame(d)
         #if(max(ANNOT)>1) message("Matrix is not binary") # Meaning the SNP overlaps several regions of this annotation
         d[d >0] <- 1
         ANNOTreformat=split(d, rownames(d))
         ANNOTreformat=lapply(ANNOTreformat, t)

          #names(ANNOTreformat) = names(listANNOT)
          if (length(ANNOTreformat)!=length(names(listANNOT))) stop("List does not have correct elements")
          return(ANNOTreformat)
}




# Function written by Panagiotis Roussos in PAINTOR pipeline
# This function gets the output of PAINTOR and estimates Log2RR and P value
# GWAS where the file GAMMA_LIKELIHOOD_CAUSAL_SNPS_CAUSAL_SNPS file is
# ESTIMATE_STAT_PAINTOR(GWAS = "PGC_BIP32b_log8cutoff", CAUSAL_SNPS = 2)

ESTIMATE_STAT_PAINTOR = function(
  DIR="/sc/orga/projects/roussp01a/PAINTOR/",
  CAUSAL_SNPS = NULL) {

    DATA0 = read.table(paste(DIR, "/results/GAMMA_LIKELIHOOD_CAUSAL_SNPS_", CAUSAL_SNPS, sep = ""), header = TRUE, fill = TRUE)
    MISSING_ANNOT = nrow(DATA0[is.na(DATA0$Baseline_Gamma),])

    DATA = na.omit(DATA0[DATA0$Annotation!="NO_ANNOTATION", ])
    BASELINE_LIKELIHOOD = DATA0[DATA0$Annotation=="NO_ANNOTATION", ]$Likelihood

    paste("Number of missing annotations: ", MISSING_ANNOT, sep = "")
    paste("Number of included annotations: ", nrow(DATA), sep = "")

    RESULTS = NULL
    for (i in 1:nrow(DATA)) {
       annotation = as.character(DATA$Annotation[i])
       dat0 = DATA[DATA$Annotation == annotation, ]
       LRT = -2 * (BASELINE_LIKELIHOOD - dat0$Likelihood)
       P = pchisq(LRT, df = 1, lower.tail=FALSE)
       BASELINE_PP = 1 / (1 + exp(dat0$Baseline_Gamma))
       ANNOT_PP = 1 / (1 + exp(dat0$Baseline_Gamma + dat0$Annotation_Gamma))
       REL_PROB = ANNOT_PP / BASELINE_PP
       LOG2_REL_PROB = log2(REL_PROB)
       out0 = data.frame(Annotation = annotation, RELATIVE_PROB = REL_PROB, LOG2_RELATIVE_PROB = LOG2_REL_PROB, LRT = LRT, Pvalue = P)
       if (is.null(RESULTS)) RESULTS = out0 else RESULTS = rbind(RESULTS, out0)
    }
    RESULTS$FDR = p.adjust(RESULTS$Pvalue, method = "fdr")
    write.table(RESULTS, file = paste(DIR, "/results/RELATIVE_RISK_PVALUE_CAUSAL_SNPS_", CAUSAL_SNPS, sep = ""), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
}

#########################################
## This function is used in reformatting data
## It takes in the file with a header and outputs the most likely position of useful columns in the data

sniff <- function(file="/sc/orga/projects/epigenAD/coloc/CARDIoGRAMplusC4D/rawData/cad.add.160614.website.txt", eqtl = FALSE) {

  # check if there is a header
  if (all(sapply(read.table(file, header = F, nrow = 1, stringsAsFactors=F), class)=="character")) header=TRUE else header=FALSE
  # Column called "F" interpreted as a logical, make sure thi sis not happening:
  if (!header) {
     checkIflogical=sapply(read.table(file, header = F, nrow = 1, stringsAsFactors=F), class)
     if (length(checkIflogical[checkIflogical=="logical"])>0) checkIflogical[checkIflogical=="logical"]="character"
     if (all(checkIflogical=="character")) header=TRUE else header=FALSE
  }
  # For cc need: 
  # SNPID, CHR, POS, F, Z, NCASE, NCONTROL or
  # SNPID, CHR, POS, SE, Z
  # but ro compute Z need either pvalue and beta, or beta and SE so mucst import those too
  # Try to find correct columns to extract
  if (!header) stop("There is no header in the data")

      line = read.table(file, header = T, nrow = 100, stringsAsFactors=F)
      if (ncol(line)==1) {
          separator=","
          line = read.table(file, header = T, nrow = 100, stringsAsFactors=F, sep=separator)
      }
      SNPID = which(grepl("snp|SNP|MarkerName", names(line),  ignore.case = T))
      CHR = which(grepl("^chr$|^CHR$|^Chromosome$|^chr_name$", names(line),  ignore.case = T))
      POS = which(grepl("^pos$|^bp$|^bp_hg19$|^position|^chrom_start$", names(line),  ignore.case = T))
      BETA = which(grepl("^beta$|^b$|^effect$|^or$|OR.A1|BMIadjMainEffects", names(line),  ignore.case = T))
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
      # If there is no "chr" column see if I can extract it from the SNP name (only if all SNP names contain chr info)
      # if (length(chr)==0) & all(grepl("chr",line[,snp])) {
      #  message("Try to extract chr info from file") # sapply(strsplit(as.character(basename(files)), "_",fixed = TRUE), "[", 1)

      # Necessary columns
      # if only chr, pos are missing could match with biomaRt to find rsid
      if (all(c(length(CHR), length(POS))==0) & length(SNPID)!=0) print("No chr, pos: must find from rsid")
      #if (any(c(length(snp), length(chr), length(pos), length(effect), length(pval))==0)) stop("Column names ", c("snp", "chr", "pos", "effect", "pval")[c(length(snp), length(chr), length(pos), length(effect), length(pval))==0], " not recognized")
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


# This function takes in a named vector with columns to extract and a character vector with corresponding class for each name (from function "sniff") and outputs the bash command to import the columns and reads the data frame (quicker than reading in R)
#  addColsName=c("GENE_SYMBOL", "ENSEMBL")
# addColsClasses = c("character", "character")
quickRead <- function(file="/sc/orga/projects/epigenAD/coloc/CARDIoGRAMplusC4D/rawData/cad.add.160614.website.txt", importCols = colsAll[[1]], colClasses = colsAll[[2]], addFilter=TRUE, colToBeFiltered=13, valueToBeFiltered=44, addColsName=NA, addColsClasses=NA) {

  # This is so the imported IDs are not in scientific format
  options(scipen=999)

  if (all(sapply(read.table(file, header = F, nrow = 1, stringsAsFactors=F), class)=="character")) header=TRUE else header=FALSE
  if (!header) {
     checkIflogical=sapply(read.table(file, header = F, nrow = 1, stringsAsFactors=F), class)
     if (length(checkIflogical[checkIflogical=="logical"])>0) checkIflogical[checkIflogical=="logical"]="character"
     if (all(checkIflogical=="character")) header=TRUE else header=FALSE
  }

  iszip = grepl(".gz|zip", file)
  if (iszip) {
      nrec = as.numeric(system(paste('zcat ', file, ' | wc -l', sep=''), intern=TRUE))
      } else {
      nrec = as.numeric(system(paste('cat ', file, ' | wc -l', sep=''), intern=TRUE))
  }
      if (header) nrec=nrec -1

   if (all(!is.na(addColsName))) {
      line = read.table(file, header = T, nrow = 1, stringsAsFactors=F)
      #if (is.na(addColsClasses)) addColsClasses = rep("character", length(addColsName)) # not convenient, should specify correct columns to make it quicker
      addColsName0 = c()
      for (i in addColsName) {
            col = which(grepl(paste(i, "$", sep=""), names(line),  ignore.case = F))
            if (length(col)>1) stop("Name ", i, " matches more than one column")
            addColsName0 = c(addColsName0, col)
            }
      names(addColsName0) = addColsName
      importCols = c(importCols, addColsName0)
      colClasses = c(colClasses, addColsClasses)
   }

   names = names(importCols)
   importCols=paste(importCols, collapse=",$")

    if (grepl("csv", file)) seperator="," else seperator=" "

# could add this to read in characters but it is just too slow... | awk '{gsub(/\"/, \"\")}1 
    if (iszip) {
            cmd = paste("zcat ", file, " | awk -F'", seperator, "' '{print $", importCols, "}' ", sep="")
            } else {
            cmd = paste("awk -F'", seperator, "' '{print $", importCols, "}' ", file, sep="")
    }

    if (addFilter) {
       if (iszip) {
            cmd = paste("zcat ", file, " | awk -F'", seperator, "' -v filter=", valueToBeFiltered, " '{if ($", colToBeFiltered, "> filter) print $", importCols, "}' ", sep="")
       } else {
            cmd = paste("awk -F'", seperator, "' -v filter=", valueToBeFiltered, " '{if ($", colToBeFiltered, "> filter) print $", importCols, "}' ", file, sep="")
       }
    }

 
     print(cmd)
     if (seperator==",") {
             data <- tryCatch(read.table(pipe(cmd), stringsAsFactors=FALSE, colClasses =colClasses, col.names=names, nrows= nrec, comment.char="", header=header, sep=" "), error=function(e) NULL )
     } else {
             data <- tryCatch(read.table(pipe(cmd), stringsAsFactors=FALSE, colClasses =colClasses, col.names=names, nrows= nrec, comment.char="", header=header), error=function(e) NULL )
    }
     # If can't read it is probably because of colClasses, so try it without
     # This makes it A LOT slower
     if ( is.null(data ) ) {
         message("Problems reading the data: trying without colClasses. NOTE: This is MUCH slower!")
         data <- tryCatch(read.table(pipe(cmd), stringsAsFactors=FALSE, col.names=names, nrows= nrec, comment.char="", header=header), error=function(e) NULL )
     }
     if ( is.null(data ) )  stop("Enter correct columns")
     return(data)
}


# Match rs to 1000 genomes
# First column must be RSID!
#  data= read.table("/sc/orga/projects/epigenAD/coloc/liver_Schadt/eQTLs/fgwas/fgwas_Liver_individual_files//chr1_Liver_1:1011087:CG_C_Probe10025907601_C1orf159.tab", header=T, sep="\t")
from_rs_to_chrpos = function(data =data) {
       message("First column must be RSID!")
       removed_list <- data.frame(Marker_removed = character(), reason = character())
       rsid.data = data[grep("^rs", data[,1]),1]
       require(biomaRt)
       max <- 100000
       x <- seq_along(rsid.data)
       rsid_splits <- split(rsid.data, ceiling(x/max))
       pos<- data.frame()
       for (i in 1:length(rsid_splits)) {
          print(paste('Split rsids number ', i, sep=""))


       # You can get the list of the available GRCh37 marts by running the following command:
       # listMarts(host='grch37.ensembl.org')
       # The following command will give you the list of datasets available for the GRCh37 snp mart:
       #grch37_snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org")
       #listDatasets(grch37_snp)

          # To query human GRCh37 snp mart, you can just use one of the following command:
          grch37_snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", dataset="hsapiens_snp")
          #snp.db = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
          snp.db <- useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", dataset="hsapiens_snp")

          #pos = getBM(c('refsnp_id','chr_name','chrom_start'), filters = c('snp_filter'), values = rsid.data, mart = snp.db )
          pos_split = getBM(c('refsnp_id','chr_name','chrom_start'), filters = c('snp_filter'), values = rsid_splits[i], mart = snp.db )

          # pos = getBM(c('refsnp_id','allele', 'chr_name','chrom_start'), filters = c('snp_filter'), values = rsid.data, mart = snp.db )
          # pos$alleles <- gsub("/", "_", pos$allele)
          # Add the alleles to chr:pos:al1_al2   ex. 1:111864465:TC_T
          # pos$chr_pos_alleles <- paste(pos$chr_name, pos$chrom_start, pos$alleles, sep=":")
          pos <- rbind(pos, pos_split)
          print(dim(pos))
       }

       pos.na <- pos$refsnp_id[which(is.na(pos$chrom_start))]  # Some results from bomart come out as 'NA'
       pos <- subset(pos, !pos$refsnp_id %in% pos.na)
       no.match.genome<- as.character(rsid.data[!rsid.data %in% pos$refsnp_id])
       # Write the removed SNPs in a .txt file
       if (length(no.match.genome)>0) {
          removed_list <- rbind(removed_list, data.frame(Marker_removed = c(no.match.genome,pos.na), reason = "Does not map to genome"))
       }
       # For now remove SNPs that map to >1 region:
       multiple.match <- names(which(table(pos$refsnp_id) > 1))
       # Write the removed SNPs in .txt file
       if (length(multiple.match)>0) {
           removed_list <- rbind(removed_list, data.frame(Marker_removed = multiple.match, reason = "Map to >1 pos in genome"))
       }
       # Take out the multiple matches

       pos <- subset(pos, !pos$refsnp_id %in% multiple.match)
       pos$chr_pos <- paste(pos$chr_name, pos$chrom_start, sep=":")
       # Remove these SNPs from original dataset
       to.remove <-  data[,1] %in% (c(multiple.match,no.match.genome))
       data <- data[!to.remove,]
       data <- merge(data, pos[,c(1,4)], by.x=1, by.y="refsnp_id", sort = FALSE)
       # biomarker.data$chr_pos <- ifelse(is.na(biomarker.data$chr_pos), biomarker.data[,1], biomarker.data$chr_pos) 
       # biomarker.data$chr_pos <- ifelse(to.replace, pos$chr_pos_alleles, biomarker.data[,1]) 

       # Reorder so you have chr_pos in first column and p-value in second column
       # find.col.df1 = c(which(names(data)=="chr_pos"), 2, which(names(data)==input.marker.name))
       # data <- data[,find.col.df1]
       #This is if want to keep all columns in dataframes: Must change function if want this!
       chrpos.col = ncol(data)
       ncol = 1:ncol(data[-chrpos.col])
       data <- data[,c(chrpos.col, ncol)]

       # Must either change the name of original biomarker data or find a very unique name for marker--not 'SNP' b/c it messes up the merging
       names(data)[2] <- "input_name"
   
       return(list(data, removed_list))
 
}


# Take a data frame with SNPID, CHR, POS columns and give SNPID as chr:pos for SNPs and chr:pos:alleles for indels (chr as numeric)
# can also filter to exclude sex and indels
# e.g.
# data <- tryCatch(read.table(pipe("awk '{print $6,$1,$2}' /sc/orga/projects/epigenAD/coloc/CMC//rawData/eQTL.GENE_FDR_0p05.db"), stringsAsFactors=FALSE, colClasses=c("character", "character", "numeric"), col.names=c("SNPID", "CHR", "POS"), header=T), error=function(e) NULL )
format_indels <- function(data=data, no_sex=FALSE, no_indels=FALSE, onlyChrPos=TRUE, changeFormat=FALSE) {

# if changeFormat = TRUE, it changes format of INDELS from chr:pos:alleles to che:pos:I or chr:pos:D, but this will depend on which allele is conaidered as A1 and which as A2
# if onlyChrPos = TRUE, only chr:pos, so indels are treated like other SNPs and if there are duplicates between a SNP and an indel, indel is removed and reported
# length(grep("^[0-9]{1,2}[:][0-9]+[:I]", df[,1]))
# length(grep("^[0-9]{1,2}[:][0-9]+[:][ACTG]+", data[,1]))

###########################################
### USE CHR:POS or CHR:POS:ALLELES as SNPID:
  removed_list <- data.frame(Marker_removed = character(), reason = character())

   
# CHR must be numeric
  removeChr=ifelse(any(grep("chr",data$CHR))>0, TRUE,FALSE)
  if (removeChr) (data$CHR=gsub("chr", "", data$CHR))
  # Make sex SNPs into 23 and 24?
  data$CHR=as.factor(data$CHR)
  levels(data$CHR)[levels(data$CHR)=="X"]="23"
  levels(data$CHR)[levels(data$CHR)=="Y"]="24"
  # Add a numeric chromosome column to facilitate ordering  
  data$CHR=as.numeric(as.character(data$CHR))
  # Remove sex SNPs
  #if (removeSex) (data=data[data$chrNum %in% 1:22,])


  # If have a "chr" in front in marker column: take out
  if (length(grep("^chr", data[,1]))>0) ( data$SNPID = gsub("^chr", "", data$SNPID) )
  if (length(grep("^X", data[,1]))>0) ( data$SNPID = gsub("^X", "23", data$SNPID) )
  if (length(grep("^Y", data[,1]))>0) ( data$SNPID = gsub("^Y", "23", data$SNPID) )

  if (length(grep(":SNP$", data[,1]))>0) ( data$SNPID = gsub(":SNP", "", data$SNPID) )

  data$SNPID = as.character(data$SNPID)
  if (!(any(grep("^rs", data$SNPID)) | any(grep("^[0-9]{1,2}[:][1-9][0-9]*$", data$SNPID)))) message("Enter correct column with names of marker")

  # Check if indels are all of this format "5:137454915:GAC_G", then I can use this:
  data$indels = ifelse(data$SNPID %in% grep("[^0-9]$", data$SNPID, value = TRUE, perl=TRUE), TRUE, FALSE)
  # Could be that it is only rsids, check if you find indels from alleles 1 and alleles 2
  if (all(data$indels=="FALSE")) {
      allele2 = which(grepl("^allele1$|^a1$|^al1$|^Effect_allele$|^Allele1$|^A1_effect|ALT$", names(data),  ignore.case = T)) # If this is the output from METAL, Allele1 refers to the allele found in the first study used in the meta-analysis
      allele1 = which(grepl("^allele2$|^a2$|^al2|^noneffect_allele$|^Non_Effect_allele$|^Allele2$|^A2_other|REF$", names(data),  ignore.case = T))
      if (length(allele1)>0 & length(allele2)>0) {
          data$indels = ifelse(nchar(data[,allele1])>2|nchar(data[,allele2])>2, TRUE, FALSE)
      }
  }
  # if already have SNPID all in chr:pos format don't do anything
  if (length(grep("^[0-9]{1,2}[:][1-9][0-9]*$", data$SNPID)) != nrow(data)) {

  # if either the expression data or the biomarker data end with a letter (INDEL) or start with a letter (rsid), format data:
  # if (length(grep("^[0-9]{1,2}[:][0-9]+[:A-Za-z]", biomarker.data[,1])) >0 | length(grep("^rs", biomarker.data[,1])) >0) {

  if (changeFormat) {
     if (length(grep("^[0-9]{1,2}[:][0-9]+[:][ID]$", data$SNPID)) == length(data[data$indels==TRUE,1])) stop()
     tochange = grep("^[0-9]{1,2}[:][0-9]+[:][ACTG]+", data$SNPID)
     t = data[tochange,]
     t$alleles = as.character(lapply(strsplit(as.character(t$SNPID), ":", fixed=TRUE),"[", 3))
     t$A1 = as.character(lapply(strsplit(as.character(t$alleles), "_", fixed=TRUE),"[", 1))
     t$A2 = as.character(lapply(strsplit(as.character(t$alleles), "_", fixed=TRUE),"[", 2))
     t$recode = paste(t$CHR, t$POS, sep=":")
     t$recode = ifelse(nchar(t$A1)>nchar(t$A2), paste(t$recode, "D", sep=":"), paste(t$recode, "I", sep=":"))
     data$SNPID[tochange] = t$recode
  }
# test some 
# library(biomaRt)
# snp.db <- useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", dataset="hsapiens_snp")
# values=list(1, 42722952, 42722953)
# getBM(attributes=c("refsnp_id", "allele"), filters=c("chr_name","start","end"), values=list(1, 42722952, 42722953), mart=snp.db)    

  names(data)[which(names(data)=="SNPID")] = "input_name"
  if (!onlyChrPos) {
   # use SNPID as chr:pos or chr:pos:ALLELES for indels
   data$SNPID = ifelse(data$indel, as.character(data$input_name), paste(data$CHR, data$POS, sep=":"))
  }

  if (onlyChrPos) {
   # use SNPID as chr:pos and remove indels with same positions as SNPs
   data$SNPID = paste(data$CHR, data$POS, sep=":")
   data = data[order(data$SNPID, data$indels),] 
   message("Number of indels removed because in same position as SNPs : ", nrow(data[duplicated(data$SNPID),]))
   dupl = data[duplicated(data$SNPID),]
   if (nrow(dupl)>0) {
     removed_list <- rbind(removed_list, data.frame(Marker_removed = dupl$input_name, reason = "INDEL map to same position as a SNP"))
     data=data[!duplicated(data$SNPID),]
   }
  }

  chrpos.col = ncol(data)
  ncol = 1:ncol(data[-chrpos.col])
  data <- data[,c(chrpos.col, ncol)]
  if (length(table(is.na(data$SNPID)))==2) stop("Problem: missing IDs created")
  # If all the values in marker column are chr:pos.   
  # if (length(grep("^[0-9]{1,2}[:][1-9][0-9]*$", data[,1])) == length(data[data$indels==FALSE,1]) ) ( format_chr_pos <- TRUE ) else ( format_chr_pos <- FALSE   )


  # If already have position, check that markers are given according to the GRC37/hg19 human genome assembly... How can I do this if I don't have both rsid and chr:pos??

  }

     # if (no_sex) ( sex.chr<- data$SNPID[grep("^[XxYy][:]", data$SNPID)] ) else ( sex.chr<-as.character())  ## For now remove sex chromosomes (since don't have in expression data)
     if (no_sex) ( sex.chr<- data$SNPID[data$CHR %in% c(23,24)] ) else ( sex.chr<-as.character())  ## For now remove sex chromosomes (since don't have in expression data)
     if (no_indels) ( indels<- data$SNPID[data$indels==TRUE] ) else ( indels<-as.character()) ## For now remove indels
     #if (length(sex.chr)>0) {
     #    removed_list <- rbind(removed_list, data.frame(Marker_removed = sex.chr, reason = "Map to sex chromosome: no match for this at the moment"))
     #  }
     other <- data$SNPID[grep("^[A-Za-z]", data$SNPID)] ## For now remove sex chromosomes (since don't have in expression data)
     to.remove <- data$SNPID %in% c(sex.chr,other,indels)
     if (length(other)>0) {
         removed_list <- rbind(removed_list, data.frame(Marker_removed = other, reason = "Map to MT chromosome?"))
       }
     if (length(indels)>0) {
         removed_list <- rbind(removed_list, data.frame(Marker_removed = indels, reason = "Looks like an indel, SNP with multiple alleles, or position is NA"))
       }

     if (nrow(removed_list)>0) {
         message("Removed ", nrow(removed_list), " SNPs because are sex SNPs, indels or others")
         data <- data[!to.remove,] 
     } else {
         removed_list = NA 
     }
     return(list(data, removed_list))
}

combineCHR <- function(base.folder="/sc/orga/projects/epigenAD/coloc/results/SCZ_CMC/eRNA", prefix=NULL, suffix=NULL, OutFile="summary.tab", Xchrom=FALSE) {
   # must make sure that everything is bounded in the same order for all chrs
   #allFiles= paste(temp, "/temp/chr", 1:22, "_withANNOT", sep="")
   if (is.null(prefix)) {
      prefix = as.character(lapply(strsplit(grep("chr", list.files(base.folder), fixed=TRUE, value=T), "[0-9]"), "[", 1))
      prefix = unique(prefix[!is.na(prefix)&prefix!=""])
   }
   if (is.null(suffix)) {
     suffix = as.character(lapply(strsplit(grep("chr", list.files(base.folder), fixed=TRUE, value=T), "[0-9]"), "[", 2))
     suffix = unique(suffix[!is.na(suffix)&suffix!=""])
   }
   allFiles= paste(base.folder, "/", prefix, 1:22, suffix, sep="")
 if (Xchrom) allFiles = c(allFiles, paste(base.folder, "/", prefix, 'X', suffix, sep=""))
 if (all(file.exists(allFiles))==TRUE) {
    final = as.character(read.table(allFiles[1], nrows = 1, header = FALSE, sep =' ', stringsAsFactors = FALSE))
    system(paste("head -1 ", allFiles[1], " > ", OutFile, sep=""))
    system(paste("tail -n +2 -q ", paste(allFiles, collapse=" "), " >> ", OutFile, sep=""))
    
    message("Final input file in ", OutFile)
    } else {
   stop("Files ", paste(allFiles[!file.exists(allFiles)], collapse=" ; "), " do not exist")
   }
}
