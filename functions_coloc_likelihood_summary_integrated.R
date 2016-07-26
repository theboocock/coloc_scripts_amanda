
coloc.eqtl.biom <- function(eqtl.df, biom.df, p12=1e-6, useBETA=TRUE, plot=FALSE, outfolder, prefix= "pref", save.coloc.output=FALSE, match_snpid=FALSE) { # match_snpid match by SNPID or match by best combination?
  if (class(eqtl.df$ProbeID)!="character") stop("When reading the data frame, make sure class of ProbeID in eQTL data is a character")

  source("/hpc/users/giambc02/scripts/COLOC/claudia.R")
  source("/hpc/users/giambc02/scripts/COLOC/optim_function.R")

if (!file.exists(outfolder)) dir.create(outfolder)
if (plot) {
   plot.fld = paste(outfolder, "plot/", sep="")
   pval.fld= paste(outfolder, "pval/", sep="")
   dir.create(file.path(plot.fld), showWarnings = FALSE)
   dir.create(file.path(pval.fld), showWarnings = FALSE)

   # For plotting with locuszoom
   refFlat_path = "/hpc/users/giambc02/scripts/locuszoom/refFlat.RData"
   source("/hpc/users/giambc02/scripts/locuszoom/call_locuszoom3_temp2.R")
   load(refFlat_path)
   refFlatRaw <- refFlatRaw.VP
}
###
  # Set Variables 
  maf_filter = 0.001 # 0.05  #MAF filter applied to datasets
  rsq_filter = 0.3 #Imputation quality filter applied to datasets

###
if ("Ncases" %in% names(biom.df)) cc=TRUE else cc=FALSE
#if (all(c("CHR", "POS") %in% names(biom.df))) haveCHRPOS.biom=TRUE else haveCHRPOS.biom=FALSE
#if (all(c("CHR", "POS") %in% names(eqtl.df))) haveCHRPOS.eqtl=TRUE else haveCHRPOS.eqtl=FALSE
maf.eqtl = ifelse("MAF" %in% names(eqtl.df), TRUE, FALSE)
maf.biom = ifelse("MAF" %in% names(biom.df), TRUE, FALSE)
if (!maf.eqtl & !maf.biom) message("There is no MAF information in neither datasets, looking for frequency column in eQTL data")

## check all columns exist
#if (useBETA) cols.eqtl = c("SNPID", "CHR", "POS", "BETA", "SE", "PVAL", "ProbeID", "N") else cols.eqtl = c("SNPID", "CHR", "POS", "PVAL", "ProbeID", "N")
#if (!all(  cols.eqtl %in% names(eqtl.df))) stop("These columns are missing from the eQTL data: ", cols.eqtl[!cols.eqtl %in% names(eqtl.df)])
#if (useBETA) cols.biom = c("SNPID", "CHR", "POS", "BETA", "SE", "PVAL", "N") else cols.biom = c("SNPID", "CHR", "POS", "PVAL", "N")
#if (cc) cols.biom = c(cols.biom, "Ncases")
#if (!all(  cols.biom %in% names(biom.df))) stop("These columns are missing from the biomarker data: ", cols.biom[!cols.biom %in% names(biom.df)])

#if ("PVAL" %in% names(biom.df))
if (useBETA) {
   cols.eqtl = c("SNPID", "CHR", "POS", "PVAL", "BETA", "SE", "ProbeID") # We need the N only if we do the sdYest step...
   cols.biom = c("SNPID", "CHR", "POS", "PVAL", "BETA", "SE") # We need the N only if we do the sdYest step...
   }
if (!useBETA) {
   cols.eqtl = c("SNPID", "CHR", "POS", "PVAL", "ProbeID", "N")
   cols.biom = c("SNPID", "CHR", "POS", "PVAL", "N")
   }
if (!all(  cols.eqtl %in% names(eqtl.df))) stop("These columns are missing from the eQTL data: ", cols.eqtl[!cols.eqtl %in% names(eqtl.df)])
if (cc & !useBETA) cols.biom = c(cols.biom, "Ncases")
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

################## SNPID MATCHING
# The reasoning here is that the user can give either only "SNPID", or "SNPID" and "input_names" (or we can find it as rsid or chr:pos and add it in the data as input_names)
# if there is no "input_names" column and the SNPID is not already chr:pos format, then the software will find the chr:pos and see if it matches better with the eQTL data than the SNPID provided

# if there is a "chr" in front of SNPID take out:
  #hasChrSNPID=ifelse(any(grep("chr", biom.df$SNPID[1:100]))>0, TRUE,FALSE)
  #  if (hasChrSNPID) biom.df$SNPID = gsub("chr", "", biom.df$SNPID)
  #hasCharSNPID=ifelse(any(grep("_|-|.", biom.df$SNPID[1:100]))>0, TRUE,FALSE)
  #  if (hasCharSNPID) biom.df$SNPID = gsub("_|-|.", ":", biom.df$SNPID)

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

if (!match_snpid) {
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

# match_snpid = max(length(intersect(biomSNPID, eqtlSNPID)), length(intersect(eqtlSNPID, biomSNPID)))
# Match is faster than intersect
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

} # if !match_snpid
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

# To create summary table: this allows for duplicated ProbeIDs:
# t=eqtl.df[eqtl.df$POS == ave(eqtl.df$POS, eqtl.df$ProbeID, FUN=min),]
# table(duplicated(t$ProbeID)); while this does not:
# t2 = eqtl.df[sapply(split(1:nrow(eqtl.df),eqtl.df$ProbeID),function(x) x[which.min(eqtl.df$POS[x])]),]
#t3 = do.call("cbind", list((eqtl.df[eqtl.df$POS == ave(eqtl.df$POS, eqtl.df$ProbeID, FUN=min),]), eqtl.df[eqtl.df$POS == ave(eqtl.df$POS, eqtl.df$ProbeID, FUN=max),]))

eqtlTable = FALSE
if (eqtlTable) {
library(data.table)
DT <- as.data.table(eqtl.df)
expr_table <- DT[, list(CHR= unique(CHR), START = min(POS), STOP = max(POS), minP = min(PVAL)), by = ProbeID]
message("There are ", nrow(expr_table), " ProbeIDs in the eQTL data")
}

onlyOverlap = FALSE
if (onlyOverlap) {
  library(GenomicRanges)
  my.BM.GRanges <- GRanges(seqnames = biom.df$CHR,
                              IRanges(start = biom.df$POS, end= biom.df$POS))  ##some start positions seem negative

  my.QTL.GRanges <- GRanges(seqnames = expr_table$CHR,
                              IRanges(start= ifelse(expr_table$START-gap<0, 0, expr_table$START-gap), end= expr_table$STOP+gap))

  my.overlap <- findOverlaps(query = my.BM.GRanges, subject = my.QTL.GRanges) 
  regions.of.interest <- unique(my.overlap@subjectHits)
  expr_table_interesting <- expr_table[ regions.of.interest, ]
  message("There are ", nrow(expr_table_interesting), " ProbeIDs (+/- gap) that overlap with at least 1 SNP in the biom data")
}
#  expr_table_interesting = unique(expr_table_interesting$ProbeID)

#probes = unique(eqtl.df$ProbeID)
commonChr = intersect(unique(eqtl.df$CHR), unique(biom.df$CHR))
### split biomarker data by chromosome
message('For maximum speed, split the data by chromosome first')
# if don't have chr pos info?
biom.dfByChr = split(biom.df, f=as.factor(biom.df$CHR))
eqtl.dfByChr = split(eqtl.df, f=as.factor(eqtl.df$CHR))
# eqtl.dfByProbe = split(eqtl.dfByChr, f=as.factor(eqtl.dfByChr$ProbeID))

# lapply(names(dataByChr), function(x){write.table(dataByChr[[x]], row.names = FALSE, quote = FALSE, col.names = TRUE, sep=" ", file = paste(DIR, "/data/", x, sep = ""))})
##################################################### now start the loop
# Now go over all regions that overlap between eQTL table and input.data
for (chr in commonChr){
      message("Looping through chr ", chr)
      my.chr = chr
      biom.df.chr = biom.dfByChr[[as.character(chr)]]

      eqtl.df.chr = eqtl.dfByChr[[as.character(chr)]]
      eqtl.dfByProbe = split(eqtl.df.chr, f=as.factor(eqtl.df.chr$ProbeID))
 
      #for (eqtl.region in 1:nrow(expr_table)) {
      # biom.df.chr <- subset(my_split_list_biom[[as.character(my.chr)]], POS >=  expr_table_interesting$start[ eqtl.region ] & POS < expr_table_interesting$stop[ eqtl.region ])biom.df.chr <- subset(my_split_list_biom[[as.character(my.chr)]])

      list.probes <- unique(eqtl.df.chr$ProbeID[eqtl.df.chr$SNPID %in% biom.df.chr$SNPID])

      # There can be more than one probe per gene 
      for (i in 1:length(list.probes)) {  ### find each gene with a cis-eQTL

       ProbeID = as.character(list.probes[i]) ##the character bit is important for probe names that are numbers
       #region.eqtl <- subset(eqtl.df.chr, ProbeID == as.character(list.probes[i]))
       region.eqtl <- eqtl.dfByProbe[[as.character(ProbeID)]]
       pos.start <- min(region.eqtl$POS)
       pos.end   <- max(region.eqtl$POS)
       #my.chr = unique(region.eqtl$CHR)

       #matches <- which(biom.df.chr$CHR==my.chr & biom.df.chr$POS > pos.start & biom.df.chr$POS < pos.end )
       matches <- which(biom.df.chr$POS >= pos.start & biom.df.chr$POS <= pos.end )
       region.biom <- biom.df.chr[matches, ]

       # matches <- which(my_split_list[[as.character(my.chr)]]
       # region.biom <- subset(my_split_list[[as.character(my.chr)]], biom.df[matches, ])

       # Loop over each biomarker 
       # message(ProbeID, ": ", length(matches), " snps in biomarkers. From: ", pos.start, " To: ", pos.end)

      if (cc & !useBETA) {
          type= "cc"
          #  s = proportion of individuals that are cases (cases / N)
          region.biom$s1 = region.biom$Ncases/region.biom$N
         }
      if (!cc & !useBETA) {
          type = "quant"
          region.biom$s1=rep(0.5, length(region.biom$N)) ## This will be ignored since the type is "quant"
         }

         merged.data <- merge(region.biom, region.eqtl, by = "SNPID",  suffixes=c(".biom", ".eqtl"))
         # Remove the pvalues at zero, otherwise it gives an error!
         if (!useBETA) merged.data = merged.data[merged.data$PVAL.biom>0 & merged.data$PVAL.eqtl>0,]


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

           # For now run with p-values (better for cc data)
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
                #dataset.eqtl$MAF <-  maf.eqtl[match(merged.data$SNPID, maf.eqtl$snp ) ,"maf"]
         }

         coloc.res <- coloc.abf(dataset.biom, dataset.eqtl, p12 = p12)
         suppressMessages(capture.output(coloc.res <- coloc.abf(dataset.biom, dataset.eqtl, p12 = p12)))
         #suppressMessages(capture.output(coloc.res <- coloc.abf2(dataset.biom, dataset.eqtl, p12 = p12)))
         #nsnps <- as.numeric(coloc.res$summary[1]) # I think this is wrong -- computed before the merging of dataests
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

         ## Per locus likelihood
         # Take the logsum of the 4 models
         l1 = coloc.res$results$lABF.df1
         l2 = coloc.res$results$lABF.df2
         lsum <- coloc.res$results$internal.sum.lABF # lsum = l1 + l2

         lH0.abf <- 0
         lH1.abf <-  logsum(l1) - log(nsnps)
         lH2.abf <-  logsum(l2) - log(nsnps)

         #lH4.abf <-  logsum(lsum)-log(nsnps)
         #lH3.abf <-  logdiff(lH1.abf + lH2.abf, lH4.abf)
         lH3.abf <- logdiff(logsum(l1) + logsum(l2), logsum(lsum))  - log(nsnps^2)
         lH4.abf <- logsum(lsum) -log(nsnps)

         all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)

         # Optimize to find the best parameters: if by region do here
         #lkl.frame = cbind(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
         #alphas = optim(c(2, -2, -2, -2, -2), fn, data=lkl.frame, method = "Nelder-Mead", control=list(fnscale=-1))
         #optim.alphas = exp(alphas$par)/ sum(exp(alphas$par))

         #res.temp = data.frame(ProbeID = ProbeID, Chr = my.chr, pos.start=pos.start, pos.end=pos.end, snp.biom=snp.biom, snp.eqtl=snp.eqtl, min.pval.biom=min.pval.biom, min.pval.eqtl=min.pval.eqtl, best.causal=best.causal, PP0.coloc.priors=pp0, PP1.coloc.priors=pp1, PP2.coloc.priors=pp2, PP3.coloc.priors = pp3, PP4.coloc.priors=pp4, lH0.abf=lH0.abf, lH1.abf=lH1.abf, lH2.abf=lH2.abf, lH3.abf=lH3.abf, lH4.abf=lH4.abf, p0.optim.priors=optim.alphas[1], p1.optim.priors=optim.alphas[2], p2.optim.priors=optim.alphas[3], p3.optim.priors=optim.alphas[4], p4.optim.priors=optim.alphas[5], files=NA)

         res.temp = data.frame(ProbeID = ProbeID, Chr = my.chr, pos.start=pos.start, pos.end=pos.end, nsnps = nsnps, snp.biom=snp.biom, snp.eqtl=snp.eqtl, min.pval.biom=min.pval.biom, min.pval.eqtl=min.pval.eqtl, best.causal=best.causal, PP0.coloc.priors=pp0, PP1.coloc.priors=pp1, PP2.coloc.priors=pp2, PP3.coloc.priors = pp3, PP4.coloc.priors=pp4, lH0.abf=lH0.abf, lH1.abf=lH1.abf, lH2.abf=lH2.abf, lH3.abf=lH3.abf, lH4.abf=lH4.abf, plotFiles=NA, files=NA)

         if (save.coloc.output) {
           coloc.out = paste(outfolder, "/coloc.output.perSNP/", sep="")
           if (!file.exists(coloc.out)) dir.create(coloc.out)
           write.table(x=coloc.res$results, file=paste(coloc.out, ProbeID, '_results.tab', sep=''),row.names = FALSE, quote = FALSE, sep = '\t')
           res.temp$files= as.character(coloc.out)
         }


         ############# PLOT
         if (plot & (pp4 > 0.2 | pp3 >=0.2) & nsnps > 2) {
         #if (plot & pp4 > 0.5 & nsnps > 2) {
         # For plotting, input called chr has to be numeric!
         # Make sure this is the case otherwise no plot is produced (because the locusZoom script coerces this value to numeric and NA is produced if not numeric
         my.chr = gsub("chr","", my.chr)

                pvalue_BF_df = as.data.frame(coloc.res[2])
                #region_name <- paste(ProbeID,'.', biom.names[j], ".chr", chr.name, "_", pos.start, "_", pos.end, sep= '')
                region_name <- paste(prefix, ".", ProbeID, ".chr", my.chr, "_", pos.start, "_", pos.end, sep= '')
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
                if (is.null(plotted)) unlink(image.eqtl)

                 res.temp$plotFiles= paste(as.character(paste(plot.fld, region_name, "_df1.pdf", sep="")), as.character(paste(plot.fld, region_name, "_df2.pdf", sep="")))
         }

         res.all <- rbind(res.all, res.temp)
         #res.all <- res.all[with(res.all, order(pp4, decreasing=T)),]
         outfname = paste(outfolder, prefix, '_summary.tab', sep='')
         write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')

     }
   }
   } # chr loop

   res.all <- data.frame(res.all)
   res.all$ProbeID <- as.character(res.all$ProbeID)
   res.all$snp.eqtl <- as.character(res.all$snp.eqtl)
   res.all$best.causal <- as.character(res.all$best.causal)
   res.all$plotFiles <- as.character(res.all$plotFiles)
   res.all$files <- as.character(res.all$files)

   optim.res =  paste(outfolder, 'maximization_results.txt', sep='') 
   # Optimize to find the best parameters
   lkl.frame = res.all[,c("lH0.abf", "lH1.abf", "lH2.abf", "lH3.abf", "lH4.abf")]
   alphas = optim(c(2,-2,-2,-2), fn, data=lkl.frame, method = "Nelder-Mead", control=list(fnscale=-1))
   optim.alphas = exp(alphas$par)/ sum(exp(c(alphas$par,alphas$par[2] + alphas$par[3])))
   write(paste("Model with 4 parameters: ", prefix, ": ", paste(optim.alphas, collapse =" , "), sep=""), file = optim.res, append=TRUE)

   lkl.frame = res.all[,c("lH0.abf", "lH1.abf", "lH2.abf", "lH3.abf", "lH4.abf")]
   alphas = optim(c(2, -2, -2, -2, -2), fn.pw.gwas, data=lkl.frame, method = "Nelder-Mead", control=list(fnscale=-1))
   optim.alphas = exp(alphas$par)/ sum(exp(alphas$par))
   write(paste("Model with 5 parameters: ", prefix, ": ", paste(optim.alphas, collapse =" , "), sep=""), file = optim.res, append=TRUE)

  # compute posteriors using the already computed likelihoods per locus (lH1.abf etc) and the optimized parameters
  # l1 = res.all$lH1.abf[1]; l2 = res.all$lH2.abf[1]; # nsnp = res.all$nsnp[1]
  # first = combine.abf.locus(l0.locus=res.all$lH0.abf[1], l1.locus=res.all$lH1.abf[1], l2.locus=res.all$lH2.abf[1], l3.locus=res.all$lH3.abf[1], l4.locus=res.all$lH4.abf[1], a0, a1, a2, a3, a4) 
   
   new.coloc = apply(res.all[,c("lH0.abf", "lH1.abf", "lH2.abf", "lH3.abf", "lH4.abf")], 1, function(x) combine.abf.locus(x[1],x[2], x[3], x[4], x[5], a0 = optim.alphas[1], a1 = optim.alphas[2], a2 = optim.alphas[3], a3 = optim.alphas[4], a4 = optim.alphas[5]))
   new.coloc=t(new.coloc)

   res.all = cbind.data.frame(res.all, new.coloc)

   #res.all <- res.all[with(res.all, order(pp4, decreasing=T)),]
   #outfname = paste(outfolder, prefix, '_summary.tab', sep='')
   write.table(x =  res.all , file = outfname, row.names = FALSE, quote = FALSE, sep = '\t')

   # If Gene.name is missing, use ensemblID instead, then try to retrieve name from biomaRt. 
   if (length(res.all$ProbeID[grep("ENSG", res.all$ProbeID)]) >0  & !("Gene.name" %in% names(res.all))) addGeneName = TRUE
   addGeneName= FALSE
   if (addGeneName) {
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

