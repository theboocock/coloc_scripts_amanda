# Written by Claudia Giambartolomei on 03/31/2016 
#########################################
  # eQTL:
  # mandatory (in any order): SNPID  CHR  POS  ([BETA  SE] or [PVAL])  N ProbeID # optional: Gene.name/ensemblID
  # if output from matrixEQTL, get SE from: eQTL.data$se.beta  = eQTL.data$beta / eQTL.data$t.stat
  # biom quant
  # mandatory (in any order): SNPID  CHR  POS  ([BETA  SE] or [PVAL])  N
  # biom case control
  # mandatory (in any order): SNPID  CHR  POS  ([BETA  SE] or [PVAL])  N Ncases
  # MAF in at least one of the datasets
# Please make sure that the betas refer to the same allele!
#########################################
# read sniff and quickRead functions from here: 
#  source("/hpc/users/giambc02/scripts/GENERAL/functions_pipeline.R")
# Function "sniff" is used in reformatting data
## It takes in the file with a header and outputs the most likely position of useful columns in the data
# Function "quickRead" takes in a named vector with columns to extract and a character vector with corresponding class for each name (from function "sniff") and outputs the bash command to import the columns and reads the data frame (quicker than reading in R)
##  addColsName=c("GENE_SYMBOL", "ENSEMBL")
##  addColsClasses = c("character", "character")
#########################################
# This function takes in eQTL and biomarker files and outputs columns needed for coloc
# Must specify whether it is a quantitative or a case-control study, and the sample size only if it is not found within the file name (it will give a warning if it is not found)
  # eQTL:
  # mandatory (in any order): SNPID  CHR  POS  [BETA SE, or  PVAL]  N MAF ProbeID # optional: Gene.name/ensemblID
  # biom quant
  # mandatory (in any order): SNPID  CHR  POS  [BETA SE, or PVAL]  N
  # biom case control
  # mandatory (in any order): SNPID  CHR  POS  [BETA SE, or PVAL]  N Ncases

formatColoc <- function(fname = "/Users/claudiagiambartolomei/Dropbox/MEGA/coloc/example/eQTL.GENE_FDR_0p05.db", type="cc", N=NA, Ncases=NA, info_filter=0.3, maf_filter=0.001, fread=FALSE, eqtl=FALSE) { 
# read sniff and quickRead functions from here: 
 source("/hpc/users/giambc02/scripts/GENERAL/functions_pipeline.R")

 colsAll = sniff(file=fname, eqtl=eqtl)

 # filter by HetDf > 3 before importing?
 # if ("HetDf" %in% names(colsAll[[1]])) {
 #    colToBeFiltered = colsAll[[1]][which(names(colsAll[[1]])=="HetDf")]
 #    data = quickRead(file = fname, importCols = colsAll[[1]], colClasses = colsAll[[2]], addFilter=TRUE, colToBeFiltered=colToBeFiltered, valueToBeFiltered=3)

  if (!fread) {
  data = quickRead(file = fname, importCols = colsAll[[1]], colClasses = colsAll[[2]], addFilter=FALSE)
  }
  if (fread) {
  require(data.table)
  #data = fread(fname, select = as.numeric(colsAll[[1]]), col.names=names(colsAll[[1]]), colClasses = colsAll[[2]])
  # fread select reads the columns in order so must re-order
  colsAll[[1]] = colsAll[[1]][order(colsAll[[1]])]
  #sel <- as.numeric(colsAll[[1]])
  #types of all columns
  data= fread(fname, select = colsAll[[1]], col.names=names(colsAll[[1]]), stringsAsFactors=FALSE)
  data = as.data.frame(data)
  }
  if ( is.null(data ) )  {
     cat("Enter correct column with pvalues.\n")
  }

  # Find rsid chr pos info
  #if (all(!c("chr", "pos") %in% names(data))) {
  #   biomart=FALSE
  #   hapmap = FALSE
  #   if (biomart) data=find_chrpos(data=data)[[1]] # biomaRt -- takes a while
     #if (hapmap) {
     # map = /sc/orga/projects/roussp01a/resources/imputation/hapmap3_r2_b36/hapmap3_r2_b36_chr18.legend
     #}
  #}
  
  #if ("ProbeID" %in% names(colsAll[[1]])) {
  #  names(data)[which(names(data)=="ProbeID")] = "ProbeID" # If there is no ProbeID use the ensemblID as the ProbeID
  #}
  #if ((!"ProbeID" %in% names(colsAll[[1]])) & ("ensemblID" %in% names(colsAll[[1]]))) {
  #  names(data)[which(names(data)=="ensemblID")] = "ProbeID" # If there is no ProbeID use the ensemblID as the ProbeID
  #}

 # what to do about the sample size?
  if (type=="quant") {
  if (!("N" %in% names(colsAll[[1]]))) {
     #if (is.na(N)) stop("Please specify a sample size because it is not found within the data")
     if (!is.na(N)) {
     message("Sample size not found within the data: using ", N)
     data$N = N
     }
   }
  }
  if (type=="cc") {
  if (!all((c("Ncases", "Ncontrols") %in% names(colsAll[[1]])))) {
     #if ( is.na(Ncases) ) stop("The dataset is a case-control so must specify number of cases")
     if (!is.na(Ncases)) {
     message("Sample size not found within the data: using ", Ncases, " ", N)
     data$Ncases = Ncases
     data$N = N
     #data$Ncontrols = Ncontrols
     #data$N = data$Ncases + data$Ncontrols
     }
    }
  }

  if (length(which(names(data)=="info"))>0 & info_filter!=0 & !is.na(info_filter)) {
    data = data[data$info > info_filter,]
  } else {
  message("Cannot find info so cannot filter")
  }

  if (length(which(names(data)=="F"))>0 & maf_filter!=0 & !is.na(maf_filter)) {
    # could be either MAF or a frequency, but to filter find MAF
    data$F = as.numeric(as.character(data$F))
    data$MAF = ifelse(data$F<0.5, data$F, 1-data$F)
    data = data[data$MAF > maf_filter,]
   } else {
   message("Cannot find frequency so cannot filter")
  }
  
  # you need either effect and SE or pval for coloc
  if ( all(!(c("BETA", "SE") %in% names(data))) | !("PVAL" %in% names(data)) ) stop("Need either BETA and SE, or PVAL in the data")
  #outfile=paste(basename(fname), ".formatted.txt", sep="")
  #write.table(x = data, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
  return(data)
}


# data is a data frame with SNPID column
find_chrpos <- function(data=data) {

####################################
# IF DON'T HAVE IT, FIND CHR AND POS IN HG19 POSITION FOR THE GIVEN RSIDS 
#################################### 
  # Match rs to 1000 genomes
  removed_list <- data.frame(Marker_removed = character(), reason = character())
       #input.marker.name = names(biomarker.data) [1]
       rsid.data = data[grep("^rs", data$SNPID),1]
       #biomarker.data[which(!biomarker.data$SNPID %in% rsid.data),] # SNPID is SNP_A-2097957

       require(biomaRt)
       max <- 100000
       x <- seq_along(rsid.data)
       rsid_splits <- split(rsid.data, ceiling(x/max))
       pos<- data.frame()
       for (i in 1:length(rsid_splits)) {
          print(paste('Split rsids number ', i, sep=""))
          grch37_snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org", dataset="hsapiens_snp")
          #snp.db <- useMart("snp", dataset="hsapiens_snp")
          #pos = getBM(c('refsnp_id','chr_name','chrom_start'), filters = c('snp_filter'), values = rsid.data, mart = snp.db )
          pos_split = getBM(c('refsnp_id','chr_name','chrom_start'), filters = c('snp_filter'), values = rsid_splits[i], mart = grch37_snp )

          # pos = getBM(c('refsnp_id','allele', 'chr_name','chrom_start'), filters = c('snp_filter'), values = rsid.data, mart = snp.db )
          # pos$alleles <- gsub("/", "_", pos$allele)
          # Add the alleles to chr:pos:al1_al2   ex. 1:111864465:TC_T
          # pos$chr_pos_alleles <- paste(pos$chr_name, pos$chrom_start, pos$alleles, sep=":")
          pos <- rbind(pos, pos_split)
          print(dim(pos))
       }

       pos.na <- pos$refsnp_id[which(is.na(pos$chrom_start))]  # Some results from bomart come out as 'NA'
       pos <- subset(pos, !pos$refsnp_id %in% pos.na)
       no.match.genome<- rsid.data[!rsid.data %in% pos$refsnp_id]
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
       to.remove <- data[,1] %in% (c(multiple.match,no.match.genome))
       data <- data[!to.remove,]
       data <- merge(data, pos, by.x=1, by.y="refsnp_id", sort = FALSE)
       names(data)[which(names(data)=="chr_name")] = "CHR"
       names(data)[which(names(data)=="chrom_start")] = "POS"

       #SNPID, CHR, POS, then BETA.biom, SE.biom, PVAL.biom, N.biom in biom_chrXX.RData where biom is a code (ht), df
       #biomarker.data=biomarker.data[,c("SNPID", "CHR", "POS", "BETA", "SE", "PVAL", "N")]
       #names(biomarker.data)[-c(1:3)]=paste(names(biomarker.data)[-c(1:3)], ".height", sep="")
       output=list(data, removed_list)
       return(output)

   #write.table(x = removed_list, file = paste(fn, "_chrpos.txt", sep=""), row.names = FALSE, quote = FALSE, sep = '\t')
}

#####################################
# This function takes in:
#  eqtl.df -  eQTL data frame
#  biom.df -  biomarker data frame
  # eQTL:
  # mandatory (in any order): SNPID  CHR  POS  ([BETA  SE] or [PVAL])  N ProbeID # optional: Gene.name/ensemblID
  # if output from matrixEQTL, get SE from: eQTL.data$se.beta  = eQTL.data$beta / eQTL.data$t.stat
  # biom quant
  # mandatory (in any order): SNPID  CHR  POS  ([BETA  SE] or [PVAL])  N
  # biom case control
  # mandatory (in any order): SNPID  CHR  POS  ([BETA  SE] or [PVAL])  N Ncases
  # MAF in at least one of the datasets
  #!! Please make sure that the betas refer to the same allele!
#  p12        -  probability of trait 1 and trait 2 
#  useBETA  - should BETA and SE be used for coloc, or only PVAL?
#  outfile   -  name of output file

coloc.eqtl.biom <- function(eqtl.df, biom.df, p12=1e-6, useBETA=TRUE, plot=TRUE, outfolder, prefix= "pref") {

# Use the function "format_indels" from here /hpc/users/giambc02/scripts/GENERAL/functions_pipeline.R
source("/hpc/users/giambc02/scripts/GENERAL/functions_pipeline.R")

suppressPackageStartupMessages({
  require(coloc);
});

if (!file.exists(outfolder)) dir.create(outfolder)
if (plot) {
   plot.fld = paste(outfolder, "plot/", sep="")
   pval.fld= paste(outfolder, "pval/", sep="")
   dir.create(file.path(plot.fld), showWarnings = FALSE)
   dir.create(file.path(pval.fld), showWarnings = FALSE)
}

###
  # Set Variables 
  maf_filter = 0.001 # 0.05  #MAF filter applied to datasets
  rsq_filter = 0.3 #Imputation quality filter applied to datasets

###
if ("Ncases" %in% names(biom.df)) cc=TRUE else cc=FALSE
maf.eqtl = ifelse("MAF" %in% names(eqtl.df), TRUE, FALSE)
maf.biom = ifelse("MAF" %in% names(biom.df), TRUE, FALSE)
if (!maf.eqtl & !maf.biom) message("There is no MAF information in neither datasets, looking for frequency column in eQTL data")

## check all columns exist
if (useBETA) cols.eqtl = c("SNPID", "CHR", "POS", "BETA", "SE", "PVAL", "ProbeID", "N") else cols.eqtl = c("SNPID", "CHR", "POS", "PVAL", "ProbeID", "N")
if (!all(  cols.eqtl %in% names(eqtl.df))) stop("These columns are missing from the eQTL data: ", cols.eqtl[!cols.eqtl %in% names(eqtl.df)])
if (useBETA) cols.biom = c("SNPID", "CHR", "POS", "BETA", "SE", "PVAL", "N") else cols.biom = c("SNPID", "CHR", "POS", "PVAL", "N")
if (cc) cols.biom = c(cols.biom, "Ncases")
if (!all(  cols.biom %in% names(biom.df))) stop("These columns are missing from the biomarker data: ", cols.biom[!cols.biom %in% names(biom.df)])

#####################
# Filter by imputation quality if column exists
info.columns <- grep( names(biom.df), pattern = 'info', value=TRUE)
if (length(info.columns) > 0)        {
    biom.df = subset(biom.df, biom.df[,info.columns] > rsq_filter)
    }
info.columns <- grep( names(eqtl.df), pattern = 'info', value=TRUE)
if (length(info.columns) > 0)        {
    eqtl.df = subset(eqtl.df, eqtl.df[,info.columns] > rsq_filter)
    }

# use only one of the MAFs from the two datasets
# First check if there is a MAF in eQTL data and use this, if not take the one in biom data
# Filter by MAF
if (maf.eqtl) {
   cols.eqtl = c(cols.eqtl, "MAF")
   eqtl.df = subset(eqtl.df, eqtl.df$MAF > maf_filter)
   }
if (!maf.eqtl & maf.biom) {
   cols.biom = c(cols.biom, "MAF")
   biom.df = subset(biom.df, biom.df$MAF > maf_filter)
   }

# Otherwise look for a frequency column
if (!maf.eqtl & !maf.biom) {
   freq.eqtl = ifelse("F" %in% names(eqtl.df), TRUE, FALSE)
   if (freq.eqtl) {
   #"^F$|freq|FRQ|MAF"
   eqtl.df$MAF = ifelse(eqtl.df$F<0.5, eqtl.df$F, 1-eqtl.df$F)
   eqtl.df = subset(eqtl.df, eqtl.df$MAF > maf_filter)
   maf.eqtl = TRUE
   cols.eqtl = c(cols.eqtl, "MAF")
   }
}

if (!maf.eqtl & !maf.biom) stop("There is no MAF information in neither datasets")


################## SNPID MATCHING # The reasoning here is that the user can give either only "SNPID", or "SNPID" and "input_names", and if there is no "input_names" column thenthe software will find the chr:pos and see if it matches better with the eQTL data than the SNPID provided
#eqtl.df = eqtl.df[,cols.eqtl]
#biom.df = biom.df[,cols.biom]
if (!"input_name" %in% colnames(biom.df)) {
      biom.df = biom.df[,cols.biom]
      message("Number of rows in original biomarker data ", nrow(biom.df))
      biom.df.form = format_indels(data=biom.df, no_sex=FALSE, no_indels=FALSE, onlyChrPos=FALSE, changeFormat=FALSE)
      biom.df = biom.df.form[[1]]
      #write.table(biom.df.form[[2]], file=paste(OUT_DIR, "log_removed_snps", sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE, sep="\t")
      message("Number of rows in formatted biomarker data ", nrow(biom.df))
}

if (!"input_name" %in% colnames(eqtl.df)) {
      eqtl.df = eqtl.df[,cols.eqtl]
      message("Number of rows in original eQTL data ", nrow(eqtl.df))
      eqtl.df.form = format_indels(data=eqtl.df, no_sex=FALSE, no_indels=FALSE, onlyChrPos=FALSE, changeFormat=FALSE)
      eqtl.df = eqtl.df.form[[1]]
      #write.table(biom.df.form[[2]], file=paste(OUT_DIR, "log_removed_snps", sep=""), row.names = FALSE, col.names=TRUE, quote = FALSE, sep="\t")
      message("Number of rows in formatted eQTL data ", nrow(eqtl.df))
}

# Find the combinations of SNPID that matches the most SNPs between the two datasets
# if there is a column called "input_name" and it matches more SNPs than the SNPID, take this one for either biom or eQTL
##match_snpid = max(length(biom.df$SNPID[biom.df$SNPID %in% eqtl.df$SNPID]), length(eqtl.df$SNPID[eqtl.df$SNPID %in% biom.df$SNPID])) # eQTL df has a different number of matches just because data contains the same SNP for different ProbeIDs -- doesn't matter which column to change
##match_input_snpid = max(length(biom.df$input_name[biom.df$input_name %in% eqtl.df$SNPID]), length(biom.df$input_name[biom.df$input_name %in% eqtl.df$SNPID]))
##match_input = max(length(biom.df$input_name[biom.df$input_name %in% eqtl.df$input_name]), length(biom.df$input_name[biom.df$input_name %in% eqtl.df$input_name]))
##find_best_column = which.max(c(match_snpid, match_input_snpid, match_input))

# Find the combinations of SNPID that matches the most SNPs between the two datasets
# if there is a column called "input_name" and it matches more SNPs than the SNPID, take this one for either biom or eQTL
biomSNPID = unique(biom.df$SNPID)
eqtlSNPID = unique(eqtl.df$SNPID)
match_snpid = max(length(biomSNPID[biomSNPID %in% eqtlSNPID]), length(eqtlSNPID[eqtlSNPID %in% biomSNPID])) # eQTL df has a different number of matches just because data contains the same SNP for different ProbeIDs -- doesn't matter which column to change

biomInput_name = unique(biom.df$input_name)
eqtlInput_name = unique(eqtl.df$input_name)
match_input_snpid = max(length(biomInput_name[biomInput_name %in% eqtlSNPID]), length(biomInput_name[biomInput_name %in% eqtlSNPID]))
match_input = max(length(biomInput_name[biomInput_name %in% eqtlInput_name]), length(biomInput_name[biomInput_name %in% eqtlInput_name]))
find_best_column = which.max(c(match_snpid, match_input_snpid, match_input))

if (find_best_column==1) message("Best combination is SNPID: do not change column names")
if (find_best_column==2) {
    message("Best combination is SNPID in one dataset and input_name in the other: change one of the column names from input_name to SNPID")
    names(biom.df)[names(biom.df)=="SNPID"] <- "SNPID2"
    names(biom.df)[names(biom.df)=="input_name"] <- "SNPID"
}
if (find_best_column==3) {
    message("Best combination is input_name in both the datasets: change both of the column names from input_name to SNPID")
    names(biom.df)[names(biom.df)=="SNPID"] <- "SNPID2"
    names(biom.df)[names(biom.df)=="input_name"] <- "SNPID"
    names(eqtl.df)[names(eqtl.df)=="SNPID"] <- "SNPID2"
    names(eqtl.df)[names(eqtl.df)=="input_name"] <- "SNPID"
}

message("There are in total ", length(biom.df$SNPID[(biom.df$SNPID%in%eqtl.df$SNPID )]), " SNPs overlapping biom and eqtl and ", length(biom.df$SNPID[!(biom.df$SNPID%in%eqtl.df$SNPID )]), " not overlapping")

#    biom_input_match_eqtl_snpid = length(biom.df$input_name[biom.df$input_name %in% eqtl.df$SNPID])
#    if (biom_input_match_eqtl_snpid>match_snpid) {
#        names(biom.df)[names(biom.df)=="SNPID"] <- "SNPID2"
#        names(biom.df)[names(biom.df)=="input_name"] <- "SNPID"
#        }

#    eqtl_input_match_biom_snpid = length(biom.df$input_name[biom.df$input_name %in% eqtl.df$SNPID])
#        if (eqtl_input_match_biom_snpid>match_snpid) {
#        names(eqtl.df)[names(eqtl.df)=="SNPID"] <- "SNPID2"
#        names(eqtl.df)[names(eqtl.df)=="input_name"] <- "SNPID"
#        }

#####################
eqtl.df = eqtl.df[,cols.eqtl]
biom.df = biom.df[,cols.biom]

# Remove missing data
eqtl.df = eqtl.df[complete.cases(eqtl.df),]
biom.df = biom.df[complete.cases(biom.df),]

# if there is a "chr" in front of CHR column
hasChr=ifelse(any(grep("chr",eqtl.df$CHR))>0, TRUE,FALSE)
  if (!hasChr) (eqtl.df$CHR=paste("chr", eqtl.df$CHR, sep=""))
hasChr=ifelse(any(grep("chr",biom.df$CHR))>0, TRUE,FALSE)
  if (!hasChr) (biom.df$CHR=paste("chr", biom.df$CHR, sep=""))

# Try matching simply by chr:pos for both datasets and take this if everything else fails
#match_snpid = max(length(biom.df$SNPID[biom.df$SNPID %in% eqtl.df$SNPID]), length(eqtl.df$SNPID[eqtl.df$SNPID %in% biom.df$SNPID]))
#if ((match_snpid)==0) { # change the format of the data with rsids
#     source("/hpc/users/giambc02/scripts/GENERAL/functions_pipeline.R")     

#####################

   res.all <- data.frame()

  ############################### now start the loop
  list.probes <- unique(eqtl.df$ProbeID)

  # There can be more than one probe per gene 
  for (i in 1:length(list.probes)) {  ### find each gene with a cis-eQTL

    ProbeID = as.character(list.probes[i]) ##the character bit is important for probe names that are numbers
    region.eqtl <- subset(eqtl.df, ProbeID == as.character(list.probes[i]))
    pos.start <- min(region.eqtl$POS)
    pos.end   <- max(region.eqtl$POS)
    my.chr = unique(region.eqtl$CHR)
    
      matches <- which(biom.df$CHR==my.chr & biom.df$POS > pos.start & biom.df$POS < pos.end )
      region.biom <- biom.df[matches, ]
      # matches <- which(my_split_list[[as.character(my.chr)]]
      # region.biom <- subset(my_split_list[[as.character(my.chr)]], biom.df[matches, ])

      # Loop over each biomarker 
      # message(ProbeID, ": ", length(matches), " snps in biomarkers. From: ", pos.start, " To: ", pos.end)

      if (cc) {
          type= "cc"
          #  s = proportion of individuals that are cases (cases / N)
         region.biom$s1 = region.biom$Ncases/region.biom$N
         } else {  
         type = "quant"
         region.biom$s1=rep(0.5, length(region.biom$N)) ## This will be ignored since the type is "quant"
         } 


         # merged.data <- merge(region.biom[, c("SNPID", colname.pval, colname.N, "s1", colname.beta, colname.se)], region.eqtl, by = "SNPID")
         merged.data <- merge(region.biom, region.eqtl, by = "SNPID",  suffixes=c(".biom", ".eqtl"))
         # Remove the pvalues at zero, otherwise it gives an error!
         merged.data = merged.data[merged.data$PVAL.biom>0 & merged.data$PVAL.eqtl>0,]
                  
         nsnps = nrow(merged.data)

         n_occur <- data.frame(table(merged.data$SNPID))
         dupl = merged.data[merged.data$SNPID %in% n_occur$Var1[n_occur$Freq > 1],]
         message("There are ", nrow(dupl)/2, " duplicated SNP names in the data")
         if (nrow(dupl)>0) {
          #removed_list <- rbind(removed_list, data.frame(Marker_removed = dupl$SNPID, reason = "Duplicated SNPs"))
          dupl=dupl[order(dupl$MAF, decreasing=T),]
          toremove = rownames(dupl[ !duplicated(dupl$SNPID), ])
          merged.data = merged.data[!(rownames(merged.data) %in% toremove),]
         }

         message(ProbeID, ": ", nsnps, " snps in both biomarker and eQTL data. From: ", pos.start, " To: ", pos.end)

         if (nsnps <= 2 ) ("There are not enough common snps in the region")
         if (nsnps > 2 ) {

         # For now run with p-values (better for cc data)
         # dataset.biom = list(snp = merged.data$SNPID, beta = merged.data[, colname.beta],varbeta = merged.data[,colname.se]^2,
         if (!useBETA) {
         dataset.biom = list(snp = merged.data$SNPID, pvalues = merged.data$PVAL.biom,
                         N = merged.data$N.biom, s=merged.data$s1, type = type, MAF=merged.data$MAF)
         dataset.eqtl = list(snp = merged.data$SNPID, pvalues = merged.data$PVAL.eqtl,
                           N = merged.data$N.eqtl, type = "quant", MAF=merged.data$MAF)
         } else {
         # ‘beta’ and ‘varbeta’
         dataset.biom = list(snp = merged.data$SNPID, beta = merged.data$BETA.biom, varbeta= (merged.data$SE.biom)^2,
                         N = merged.data$N.biom, s=merged.data$s1, type = type, MAF=merged.data$MAF)
         dataset.eqtl = list(snp = merged.data$SNPID, beta = merged.data$BETA.eqtl, varbeta= (merged.data$SE.eqtl)^2,
                           N = merged.data$N.eqtl, type = "quant", MAF=merged.data$MAF)
         #dataset.eqtl$MAF <-  maf.eqtl[match(merged.data$SNPID, maf.eqtl$snp ) ,"maf"]
         }
         suppressMessages(capture.output(coloc.res <- coloc.abf(dataset.biom, dataset.eqtl, p12 = p12)))
         pp3       <- as.numeric(coloc.res$summary[5])
         pp4       <- as.numeric(coloc.res$summary[6])
         snp.biom <- merged.data[which.min(merged.data$PVAL.biom), "SNPID"]
         snp.eqtl <- merged.data[which.min(merged.data$PVAL.eqtl), "SNPID"]
         min.pval.biom <- min(merged.data$PVAL.biom)
         min.pval.eqtl <- min(merged.data$PVAL.eqtl)
         best.causal = as.character(coloc.res$results$snp[which.max(coloc.res$results$SNP.PP.H4)])

         res.temp = data.frame(ProbeID = ProbeID, Chr = my.chr, pos.start=pos.start, pos.end=pos.end, snp.biom=snp.biom, snp.eqtl=snp.eqtl, min.pval.biom=min.pval.biom, min.pval.eqtl=min.pval.eqtl, best.causal=best.causal, pp3, pp4, files=NA)


         ############# PLOT
         if (plot & (pp4 > 0.2 | pp3 >=0.2) & nsnps > 2) {
         #if (plot & pp4 > 0.5 & nsnps > 2) {
         # For plotting, input called chr has to be numeric!
         # Make sure this is the case otherwise no plot is produced (because the locusZoom script coerces this value to numeric and NA is produced if not numeric
         my.chr = gsub("chr","", my.chr)

                pvalue_BF_df = as.data.frame(coloc.res[2])
                #region_name <- paste(ProbeID,'.', biom.names[j], ".chr", chr.name, "_", pos.start, "_", pos.end, sep= '')
                region_name <- paste(prefix, '.', ProbeID,'.', ".chr", my.chr, "_", pos.start, "_", pos.end, sep= '')
                pvalue_BF_file <- paste(pval.fld, 'pval_', region_name, '.txt', sep="")

                ### LocusZoom arguments:
                pvalue_BF_df$chr = my.chr
                pvalue_BF_df$pos = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "POS.biom"]
                pvalue_BF_df$locus_snp <- paste(pvalue_BF_df$chr, pvalue_BF_df$pos, sep=":")
                pvalue_BF_df$locus_snp <- paste("chr", pvalue_BF_df$locus_snp, sep="")
                # INDELS FORMAT FOR LOCUSZOOM: chr1:117930794:AAG_A (not rsid)
                pvalue_BF_df$locus_snp <- ifelse(grepl("*[:][:]*", pvalue_BF_df$results.snp), paste("chr", as.character(pvalue_BF_df$results.snp), sep=""), as.character(pvalue_BF_df$locus_snp))
                pvalue_BF_df$results.pvalues.df1 = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "PVAL.biom"]
                pvalue_BF_df$results.pvalues.df2 = merged.data[match(pvalue_BF_df$results.snp, merged.data$SNPID), "PVAL.eqtl"]
                image.biom = paste(plot.fld, "/", region_name, '_df1.pdf', sep='')
                image.eqtl = paste(plot.fld, "/", region_name, '_df2.pdf', sep='')

                write.table(x =  pvalue_BF_df , file = pvalue_BF_file, row.names = FALSE, quote = FALSE, sep = '\t')

                message('Output pdf for biomarker: ', image.biom)
                pdf(image.biom, width = 9, height = 9)
                #output of region_ld.ld is in /SAN/biomed/biomed14/vyp-scratch/vincent/eQTLs/ ?
                # If INDEL ALLELES do not match exactly (for ex. are reversed from the reference EUR files in here /cluster/project8/vyp/vincent/toolsVarious/locuszoom/EUR/), skip for now:
                plotted = tryCatch(locuszoom.ms(metal = pvalue_BF_file,
                  refSnp = pvalue_BF_df[pvalue_BF_df$results.snp==best.causal,"locus_snp"] , #rs10877835
                  title = 'A',
                  pvalCol='results.pvalues.df1',
                  legend = 'left',
                  markerCol='locus_snp',
                  ylab = '-log10( biomarker P-value )',
                  chrCol= 'chr',
                  posCol = 'pos',
                  chr = my.chr,
                  showGenes = TRUE,
                  show_xlab=FALSE,
                  temp.file.code = region_name,
                  start= pos.start ,
                  end = pos.end
                  ), error=function(e) NULL )
                dev.off()
                # If the plot is empty unlink:
                if (is.null(plotted)) unlink(image.biom)

                message('Output pdf for eQTL: ', image.eqtl)
                pdf(image.eqtl, width = 9, height = 9)
                plotted = tryCatch(locuszoom.ms(metal = pvalue_BF_file,
                  refSnp = pvalue_BF_df[pvalue_BF_df$results.snp==best.causal,"locus_snp"] , #rs10877835
                  title = 'B',
                  pvalCol='results.pvalues.df2',
                  legend = 'left',
                  markerCol='locus_snp',
                  ylab = '-log10( expression P-value )',
                  chrCol= 'chr',
                  posCol = 'pos',
                  chr = my.chr,
                  showGenes = TRUE,
                  show_xlab=FALSE,
                  temp.file.code = region_name,
                  start= pos.start ,
                  end = pos.end
                  ), error=function(e) NULL )
                dev.off()
                # If the plot is empty unlink:
                if (is.null(plotted)) unlink(image.biom)

                #do.call(file.remove,list(list.files(wd, pattern="region_ld")))
                #file.remove(paste(wd, region_name, "_df2", ".log", sep=""))

                # Now merge the two pdfs??
                #system(paste("convert ", region_name, "_df1.pdf", " ", region_name, "_df2.pdf", " ", region_name, ".pdf", sep=""))
                #file.remove(paste(wd, region_name, "_df1.pdf", sep=""))
                #file.remove(paste(wd, region_name, "_df2.pdf", sep=""))

                 res.temp$files= paste(as.character(paste(plot.fld, region_name, "_df1.pdf", sep="")), as.character(paste(plot.fld, region_name, "_df2.pdf", sep="")))
         }


         res.all <- rbind(res.all, res.temp)

     }
   }

   res.all <- data.frame(res.all)
   res.all$ProbeID <- as.character(res.all$ProbeID)
   res.all$snp.eqtl <- as.character(res.all$snp.eqtl)
   res.all$best.causal <- as.character(res.all$best.causal)
   res.all$files <- as.character(res.all$files)

   res.all <- res.all[with(res.all, order(pp4, decreasing=T)),]
   outfname = paste(outfolder, prefix, '_summary.tab', sep='')
   write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')

   # If Gene.name is missing, use ensemblID instead, then try to retrieve name from biomaRt. 
   # if (nrow(eqtl.data[eqtl.data$Gene.name=="",]) > 0) {
   # eqtl.data$Gene.name = ifelse(eqtl.data$Gene.name=="", eqtl.data$ensemblID, eqtl.data$Gene.name)
   if (length(res.all$ProbeID[grep("ENSG", res.all$ProbeID)]) >0  & !("Gene.name" %in% names(res.all))) {
   # if (!("Gene.name" %in% names(res.all))) {
   res.all$Gene.name = res.all$ProbeID
   biomart=FALSE # it doesn't work sometimes -- cannot connect etc
      if (biomart) {
      library(biomaRt)
      #if (length(res.all$Gene.name[grep("ENSG", res.all$Gene.name)]) >0 ) {
        library(biomaRt)
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

   return(res.all)
}


