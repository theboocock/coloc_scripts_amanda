coloc.eqtl.biom <- function(eqtl.df, biom.df, p12=1e-6, useBETA=TRUE, plot=FALSE, outfolder, prefix= "pref", save.coloc.output=FALSE) {

  source("~/psychgen/resources/COLOC2/COLOC_scripts/scripts/claudia.R")
  source("~/psychgen/resources/COLOC2/COLOC_scripts/scripts/coloc.eqtl.biom.compare.R")

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
   message(names(biom.df))
   freq.biom = ifelse("F" %in% names(biom.df), TRUE, FALSE)
   if (freq.eqtl) {
       #"^F$|freq|FRQ|MAF"
       eqtl.df$MAF = ifelse(eqtl.df$F<0.5, eqtl.df$F, 1-eqtl.df$F)
       eqtl.df = subset(eqtl.df, eqtl.df$MAF > maf_filter)
       maf.eqtl = TRUE
       cols.eqtl = c(cols.eqtl, "MAF")
   }else if(freq.biom){
       biom.df$MAF = ifelse(biom.df$F<0.5, biom.df$F, 1-biom.df$F)
       biom.df = subset(biom.df, biom.df$MAF > maf_filter)
       maf.biom= TRUE
       cols.biom = c(cols.biom, "MAF")
   }
}

if (!maf.eqtl & !maf.biom) stop("There is no MAF information in neither datasets")
################## SNPID MATCHING
  # if there is a "chr" in front of CHR column
  hasChr=ifelse(any(grep("chr",eqtl.df$CHR[1:2]))>0, TRUE,FALSE)
    if (hasChr) (eqtl.df$CHR=gsub("chr", "", eqtl.df$CHR))
  hasChr=ifelse(any(grep("chr",biom.df$CHR[1:2]))>0, TRUE,FALSE)
    if (hasChr) (biom.df$CHR=gsub("chr", "", biom.df$CHR))

  #if (!"input_name" %in% colnames(biom.df) && (length(grep("^[0-9]{1,2}[:][1-9][0-9]*$", biom.df$SNPID))!=nrow(biom.df))) {
  if (length(grep("^[0-9]{1,2}[:][1-9][0-9]*$", biom.df$SNPID))!=nrow(biom.df)) addChrposBiom = TRUE
  if (length(grep("^[0-9]{1,2}[:][1-9][0-9]*$", eqtl.df$SNPID))!=nrow(eqtl.df)) addChrposEQTL = TRUE
  if (addChrposBiom) {
    biom.df$chrpos = paste(biom.df$CHR, biom.df$POS, sep=":")
  }

  #if (!"input_name" %in% colnames(eqtl.df)) {
  if (addChrposEQTL) {
    eqtl.df$chrpos = paste(eqtl.df$CHR, eqtl.df$POS, sep=":")
  }

# Find the combinations of SNPID that matches the most SNPs between the two datasets
biomSNPID = unique(biom.df$SNPID)
eqtlSNPID = unique(eqtl.df$SNPID)
match_snpid = max(length(biomSNPID[biomSNPID %in% eqtlSNPID]), length(eqtlSNPID[eqtlSNPID %in% biomSNPID]))
match_chrpos_snpid = 0
match_chrpos = 0

if (addChrposBiom) {
  biomchrpos = unique(biom.df$chrpos)
  match_chrpos_snpid = max(length(biomchrpos[biomchrpos %in% eqtlSNPID]), length(biomchrpos[biomchrpos %in% eqtlSNPID]))
}

#if (addChrposEQTL & !addChrposBiom) {
if (addChrposEQTL) {
   eqtlchrpos = unique(eqtl.df$chrpos)
   if (!addChrposBiom) {
   match_chrpos_snpid = max(length(eqtlchrpos[eqtlchrpos %in% biomSNPID]), length(eqtlchrpos[eqtlchrpos %in% biomSNPID]))
   }
}

if (addChrposBiom & addChrposEQTL) match_chrpos = max(length(biomchrpos[biomchrpos %in% eqtlchrpos]), length(biomchrpos[biomchrpos %in% eqtlchrpos]))
find_best_column = which.max(c(match_snpid, match_chrpos_snpid, match_chrpos))

if (find_best_column==1) message("Best combination is SNPID: do not change column names")
if (find_best_column==2) {
    message("Best combination is SNPID in one dataset and input_name in the other: change one of the column names from input_name to SNPID")
    names(biom.df)[names(biom.df)=="SNPID"] <- "SNPID2"
    names(biom.df)[names(biom.df)=="chrpos"] <- "SNPID"
}
if (find_best_column==3) {
    message("Best combination is chrpos in both the datasets: change both of the column names from chrpos to SNPID")
    names(biom.df)[names(biom.df)=="SNPID"] <- "SNPID2"
    names(biom.df)[names(biom.df)=="chrpos"] <- "SNPID"
    names(eqtl.df)[names(eqtl.df)=="SNPID"] <- "SNPID2"
    names(eqtl.df)[names(eqtl.df)=="chrpos"] <- "SNPID"
}

#####################
  biom.df = biom.df[,cols.biom]
  eqtl.df = eqtl.df[,cols.eqtl]

  # Remove missing data
  eqtl.df = eqtl.df[complete.cases(eqtl.df),]
  biom.df = biom.df[complete.cases(biom.df),]

  res.all <- data.frame()

###################
# Find ProbeIDs overlapping with biom.df
gap = 2000000 # cushion from start and end of probe/gene
commonChr = intersect(unique(eqtl.df$CHR), unique(biom.df$CHR))
message('For maximum speed, split the data by chromosome first')
biom.dfByChr = split(biom.df, f=as.factor(biom.df$CHR))
eqtl.dfByChr = split(eqtl.df, f=as.factor(eqtl.df$CHR))
##################################################### now start the loop
# Now go over all regions that overlap between eQTL table and input.data
for (chr in commonChr){
      message("Looping through chr ", chr)
      my.chr = chr
      biom.df.chr = biom.dfByChr[[as.character(chr)]]

      eqtl.df.chr = eqtl.dfByChr[[as.character(chr)]]
      eqtl.dfByProbe = split(eqtl.df.chr, f=as.factor(eqtl.df.chr$ProbeID))

      list.probes <- unique(eqtl.df.chr$ProbeID)

      # There can be more than one probe per gene 
      for (i in 1:length(list.probes)) {  

       ProbeID = as.character(list.probes[i]) ##the character bit is important for probe names that are numbers
       region.eqtl <- eqtl.dfByProbe[[as.character(ProbeID)]]
       pos.start <- min(region.eqtl$POS)
       pos.end   <- max(region.eqtl$POS)

       #matches <- which(biom.df.chr$CHR==my.chr & biom.df.chr$POS > pos.start & biom.df.chr$POS < pos.end )
       matches <- which(biom.df.chr$POS > pos.start & biom.df.chr$POS < pos.end )
       region.biom <- biom.df.chr[matches, ]

      if (cc) {
          type= "cc"
          #  s = proportion of individuals that are cases (cases / N)
         region.biom$s1 = region.biom$Ncases/region.biom$N
         } else {
         type = "quant"
         region.biom$s1=rep(0.5, length(region.biom$N)) ## This will be ignored since the type is "quant"
         }

         merged.data <- merge(region.biom, region.eqtl, by = "SNPID",  suffixes=c(".biom", ".eqtl"))
         # Remove the pvalues at zero, otherwise it gives an error!
         merged.data = merged.data[merged.data$PVAL.biom>0 & merged.data$PVAL.eqtl>0,]


         n_occur <- data.frame(table(merged.data$SNPID))
         dupl = merged.data[merged.data$SNPID %in% n_occur$Var1[n_occur$Freq > 1],]
         message("There are ", nrow(dupl)/2, " duplicated SNP names in the data")
         if (nrow(dupl)>0) {
          #removed_list <- rbind(removed_list, data.frame(Marker_removed = dupl$SNPID, reason = "Duplicated SNPs"))
          dupl=dupl[order(dupl$MAF, decreasing=T),]
          toremove = rownames(dupl[ !duplicated(dupl$SNPID), ])
          merged.data = merged.data[!(rownames(merged.data) %in% toremove),]
         }

         nsnps = nrow(merged.data)

         message(ProbeID, ": ", nsnps, " snps in both biomarker and eQTL data. From: ", pos.start, " To: ", pos.end)

         if (nsnps <= 2 ) ("There are not enough common snps in the region")
         if (nsnps > 2 ) {
           if (!useBETA) {
                dataset.biom = list(snp = merged.data$SNPID, pvalues = merged.data$PVAL.biom,
                           N = merged.data$N.biom, s=merged.data$s1, type = type, MAF=merged.data$MAF)
                dataset.eqtl = list(snp = merged.data$SNPID, pvalues = merged.data$PVAL.eqtl,
                           N = merged.data$N.eqtl, type = "quant", MAF=merged.data$MAF)
           } else {
                dataset.biom = list(snp = merged.data$SNPID, beta = merged.data$BETA.biom, varbeta= (merged.data$SE.biom)^2,
                           N = merged.data$N.biom, s=merged.data$s1, type = type, MAF=merged.data$MAF)
                dataset.eqtl = list(snp = merged.data$SNPID, beta = merged.data$BETA.eqtl, varbeta= (merged.data$SE.eqtl)^2,
                           N = merged.data$N.eqtl, type = "quant", MAF=merged.data$MAF)
            }
         suppressMessages(capture.output(coloc.res <- coloc.abf(dataset.biom, dataset.eqtl, p12 = p12)))
         pp0       <- as.numeric(coloc.res$summary[2])
         pp1       <- as.numeric(coloc.res$summary[3])
         pp2       <- as.numeric(coloc.res$summary[4])
         pp3       <- as.numeric(coloc.res$summary[5])
         pp4       <- as.numeric(coloc.res$summary[6])
         snp.biom <- merged.data[which.min(merged.data$PVAL.biom), "SNPID"]
         snp.eqtl <- merged.data[which.min(merged.data$PVAL.eqtl), "SNPID"]
         min.pval.biom <- min(merged.data$PVAL.biom)
         min.pval.eqtl <- min(merged.data$PVAL.eqtl)
         best.causal = as.character(coloc.res$results$snp[which.max(coloc.res$results$SNP.PP.H4)])
         # Take the logsum of the 4 models
         l1 = coloc.res$results$lABF.df1
         l2 = coloc.res$results$lABF.df2
         lsum <- coloc.res$results$internal.sum.lABF # lsum = l1 + l2

         lH0.abf <- 0
         lH1.abf <-  logsum(l1) - log(nsnps)
         lH2.abf <-  logsum(l2) - log(nsnps)
         lH3.abf <- logdiff(logsum(l1) + logsum(l2), logsum(lsum))  - log(nsnps^2)
         lH4.abf <- logsum(lsum) -log(nsnps)

         all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
                 res.temp = data.frame(ProbeID = ProbeID, Chr = my.chr, pos.start=pos.start, pos.end=pos.end, nsnps = nsnps, snp.biom=snp.biom, snp.eqtl=snp.eqtl, min.pval.biom=min.pval.biom, min.pval.eqtl=min.pval.eqtl, best.causal=best.causal, PP0.coloc.priors=pp0, PP1.coloc.priors=pp1, PP2.coloc.priors=pp2, PP3.coloc.priors = pp3, PP4.coloc.priors=pp4, lH0.abf=lH0.abf, lH1.abf=lH1.abf, lH2.abf=lH2.abf, lH3.abf=lH3.abf, lH4.abf=lH4.abf, plotFiles=NA, files=NA)

         if (save.coloc.output) {
           coloc.out = paste(outfolder, "/coloc.output.perSNP/", sep="")
           if (!file.exists(coloc.out)) dir.create(coloc.out)
           write.table(x=coloc.res$results, file=paste(coloc.out, ProbeID, '_results.tab', sep=''),row.names = FALSE, quote = FALSE, sep = '\t')
           res.temp$files= as.character(coloc.out)
         }

         res.all <- rbind(res.all, res.temp)

     }
   }
   } # chr loop
   lkl.frame = res.all[,c("lH0.abf", "lH1.abf", "lH2.abf", "lH3.abf", "lH4.abf")]
   alphas = optim(c(2, -2, -2, -2, -2), fn.pw.gwas, data=lkl.frame, method = "Nelder-Mead", control=list(fnscale=-1))
   optim.alphas = exp(alphas$par)/ sum(exp(alphas$par))
   new.coloc = apply(res.all[,c("lH0.abf", "lH1.abf", "lH2.abf", "lH3.abf", "lH4.abf")], 1, function(x) combine.abf.locus(x[1],x[2], x[3], x[4], x[5], a0 = optim.alphas[1], a1 = optim.alphas[2], a2 = optim.alphas[3], a3 = optim.alphas[4], a4 = optim.alphas[5]))
   new.coloc=t(new.coloc)

   res.all = cbind.data.frame(res.all, new.coloc)
   outfname = paste(outfolder, prefix, '_summary.tab', sep='')
   write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')
   write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')

   return(res.all)
}

