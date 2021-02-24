
require(data.table)
require(tidyr)
require(dplyr)
require(ggplot2)

option_list = list(
   make_option(c("-w", "--workingDir"), type="character", default="./", 
               help="Working directory [default= %default]", metavar="character"),
   make_option(c("-f", "--bedpe"), type="character", default=NULL, 
               help="intercation file as bedpe", metavar="character"),
   make_option(c("-o", "--outputfile"), type="character", default=NULL, 
               help="name of the output file ", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt
if (is.null(opt$geneList)){
   print_help(opt_parser)
   stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

########## Variables######
wkdir=opt$workingDir
wkdir=paste0(wkdir,"/")
cis=opt$bedpe
out=opt$outputfile
head(cis,n=5)

###########
setwd(wkdir)
cis=fread(cis)

### Error checks to make sure cis file has required columns
if (any (is.na (cis$chrom1)))
   stop ("Please specify chromsome names for all bin1 using the 'chrom1' column in peaks file...")
if (any (is.na (cis$chrom2)))
   stop ("Please specify chromsome names for all bins using the 'chrom2' column in peaks file...")
if (any (is.na (cis$start1)))
   stop ("Please specify start for all bin1 using the 'start1' column in peaks file...")
if (any (is.na (cis$start2)))
   stop ("Please specify start for all bin2 using the 'start2' column in peaks file...")
if (any (is.na (cis$end1)))
   stop ("Please specify end for all bin1 using the 'end1' column in peaks file...")
if (any (is.na (cis$end2)))
   stop ("Please specify end for all bin2 using the 'end2' column in peaks file...")
if (any (is.na (cis$coverage_bin1)))
   stop ("Please specify coverage for all bin1 using the 'coverage_bin1' column in peaks file...")
if (any (is.na (cis$coverage_bin2)))
   stop ("Please specify coverage for all bin2 using the 'coverage_bin2' column in peaks file...")
if (any (is.na (cis$rawcounts)))
   stop ("Please specify raw intercation counts for all intercations using the 'rawcounts' column in peaks file...")


### calculate distance
cis$dist=abs(as.numeric(as.character(cis$start2))-as.numeric(as.character(cis$start1)))
cis=cis[which(cis$d<=2000000 & cis$d >5000),]
saveRDS(cis,"all_cis_2MB.rds")

### Fit NB ######### 
s=sample(1:nrow(cis),nrow(cis)*0.2)
mod_negB=MASS::glm.nb(rawcounts ~ log(coverage_bin1)+log(coverage_bin2)+log(dist),data=cis[s,])
summary(mod_negB)
pred_negB=predict(mod_negB,dispersion=1/mod_negB$theta,cis_pp,type="response")
outliers=(qnbinom(n, mu=pred_negB,size=mod_negB$theta))
cis_refit=cis[cis$rawCounts<outliers,]

s=sample(1:nrow(cis_refit),nrow(cis_refit)*0.2)
mod_negB_reFit=MASS::glm.nb(rawcounts ~ log(coverage_bin1)+log(coverage_bin2)+log(dist),data=cis[s,])
summary(mod_negB_reFit)
  
cis$pred_negB=predict(mod_negB_reFit,cis,type="response",dispersion=1/mod_negB_reFit$theta)
cis$pval_NB=pnbinom((cis$rawcounts),size=mod_negB_reFit$theta,mu=cis$pred_negB,lower.tail = F)
cis$padj_NB=p.adjust(cis$pval_NB,method = "fdr")
cis$OE_NB=cis$rawcounts/cis$pred_negB
cis$residuals=cis$rawcounts-cis$pred_negB
cis$logpadj=-log10(round(cis$padj_NB,digits=15))
colnames(cis)[1]="#chrom_2"
saveRDS(cis,paste0(out,"_all_cis_2MB_NB.rds"))

  
