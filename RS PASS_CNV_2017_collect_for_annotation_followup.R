# This is awkward but it had to be done because h1/h2 cannot currently be run automatically but only manually, by hard setting Batch="h1"/"h2" in PASS_CNV_2017_collect_for_annotation.R

# Load the content of the manual run
rm(list=ls())
setwd("~/Documents/shared/proj/R/RStudio/PASS_CNV_follow_up")
load("~/Documents/shared/proj/PASS_CNV/Analysis/Code/PASS_CNV_2017_collect_for_annotation.2018_06_18.RData")

# Clean up
rm(list=setdiff(ls(), save_object_list))

# Modify/introdduce running parameters for local execution  
Results_dir_local = "Results/"

# Overlap subs
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
PASS_CNV.overlap_sum_cov_02 <- function(xx1, xx2, db_levl = 2) {
# Initial, buggy version, not with breaks.  Use direct overlaps.  Double counts multiple overlaps
ov_size = ov_frac = rep(0, nrow(xx1))
for(ic1 in 1:nrow(xx1)) {
  s_sc    = (xx1[ic1, "Sample_Group"] == xx2[, "Sample_Group"]) & (xx1[ic1, "Chr"] == xx2[, "Chr"])
  
  sI      = findInterval(xx2[,"Start"], xx1[ic1, c("Start","End")])
  eI      = findInterval(xx2[,  "End"], xx1[ic1, c("Start","End")])
  s_ov    = !((sI+eI == 0) | (sI+eI == 4)) ## some overlap
  if(sum(s_sc & s_ov) == 0) next
  if(db_levl > 0) cat("------------------------------------------------------ ic1=", ic1, "\n")
  s = sI==1 & eI==1 & s_sc;               ov_size[ic1] = ov_size[ic1] +  sum(xx2[  s,"End"] - xx2[  s,"Start"], na.rm=T)
  if(db_levl > 1) cat("all in =", sum(s), "\n")
  s = sI==0 & eI==1 & s_sc;               ov_size[ic1] = ov_size[ic1] +  sum(xx2[  s,"End"] - xx1[ic1,"Start"], na.rm=T)
  if(db_levl > 1) cat("left   =", sum(s), "\n")
  s = sI==1 & eI==2 & s_sc;               ov_size[ic1] = ov_size[ic1] +  sum(xx1[ic1,"End"] - xx2[  s,"Start"], na.rm=T)
  if(db_levl > 1) cat("right  =", sum(s), "\n")
  s = sI==0 & eI==2 & s_sc; if(sum(s))  { ov_size[ic1] = xx1[ic1,"End"] - xx1[ic1,"Start"] }
  if(db_levl > 1) cat("all out=", sum(s), "\n")
  
  ov_size[ic1] = ov_size[ic1]/1e6
  ov_frac[ic1] = ov_size[ic1]/xx1[ic1,"Size"]
  
  if(db_levl  > 0) {
    ip = c(1,3:6, 8, 11:12)
    print(xx1[ic1,         ip]) 
    cat("------------- \n")
    print(xx2[s_sc & s_ov, ip]) 
    cat("size, fraction =", ov_size[ic1]/1e6, ov_frac[ic1], "\n")
  }
#   if(ic1>1) stop(" -- some overlap -- ")  
}
return(cbind(ov_size_02=ov_size, ov_frac_02=ov_frac))
}
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
PASS_CNV.overlap_sum_cov    <- function(xx1, xx2_in, db_levl = 2) {
  ov_size    = ov_frac = ov_n = rep(0, nrow(xx1))
  ov_ID      = array()
  ov_L_array = list()

  xx2 = xx2_in
  for(ic1 in 1:nrow(xx1)) {
    if(is.null(xx2_in)) { xx2 = xx1[setdiff(1:nrow(xx1), ic1),] }
    s_sc                         = xx1[ic1, "Chr"] == xx2[, "Chr"]
    if(!is.null(xx2_in))    s_sc = (xx1[ic1, "Sample_Group"] == xx2[, "Sample_Group"]) & s_sc  
  
    x                 = unlist(xx1[ic1, c("Start","End")])
    ov_L              = Find_overlap_regions(x, xx2[s_sc, c("Start","End")], db_levl=db_levl)
#     cat(" ---- ic1 = ", ic1, "\n")
#     cat("is.null(ov_L$ov_size)", is.null(ov_L$ov_size), "ov_L$ov_size =", ov_L$ov_size, "\n")
    ov_ID[ic1] = ""
    if(!is.null(ov_L)) { 
      ov_size[ic1]      = ov_L$ov_size/1e6
      ov_frac[ic1]      = ov_L$ov_frac
      xx                = xx2[s_sc,]
      x                 = unique(xx[ov_L$ss_table_flag, "Sample_Group"])
      if(length(x) > 0) {
        ov_n[ic1]       = length(x)
        z = NULL; for(i in 1:length(x)) z = paste(z, x[i])
        ov_ID[ic1]      = z
      }
    }
    ov_L_array[[ic1]] = ov_L
  }

  return( list(xx=cbind(ov_size, ov_frac), ov_n=ov_n, ov_ID=ov_ID, ov_L_array=ov_L_array) )
}
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
Find_overlap_regions        <- function(ss_interval, ss_table, db_levl=0) {
  # ss_interval   start stop interval
  # ss_table      start stop table
  if(nrow(ss_table) == 0) return(NULL)
  breaks  = unique(sort(c(ss_interval, ss_table[,1], ss_table[,2])))
  breaks  = breaks[ss_interval[1] <= breaks  &  breaks <= ss_interval[2]]
  mpoints = breaks[1:(length(breaks)-1)] + diff(breaks)/2
  ss_table_flag    = array(F, nrow(ss_table))
  ss_interval_flag = array()
  ss_interval_indx = list()
  for(ip in 1:length(mpoints)) {
#     ip = 1; cbind(mpoints[ip], ss_table)
    xxt  = cbind(mpoints[ip], ss_table)
    ov_r = apply(xxt, MARGIN=1, function(x) findInterval(x[1], x[2:3])); if(db_levl >= 2) { cat("---- ip=", ip, "\n"); print(cbind(xxt, ov_r)) }
    ov_r                   = ov_r == 1
    ss_table_flag          = ss_table_flag | ov_r
    ss_interval_flag[ip]   = sum(ov_r) > 0
    ss_interval_indx[[ip]] = which(ov_r)
  }
  ov_size = sum(diff(breaks)*as.integer(ss_interval_flag))
  ov_frac = ov_size/diff(ss_interval)
  
  return(ov_L=list(ov_size=ov_size, ov_frac=ov_frac,
  ss_table_flag=ss_table_flag, ss_interval_flag=ss_interval_flag, ss_interval_indx=ss_interval_indx))
}
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

# Read annotated spreadsheet from Klaus
fn = paste(RESULTS_dir, "2018_07_Klaus_Annotations/v00_2018_06_15_Cuts.I..KSA.Annotated.csv", sep="")
CNV_I = read.csv(fn,    stringsAsFactors=F, header=T)
CNV_I = CNV_I[,3:ncol(CNV_I)]
fn = paste(RESULTS_dir, "2018_07_Klaus_Annotations/v00_2018_06_15_Cuts.P..KSA.Annotated.csv", sep="")
CNV_P = read.csv(fn,    stringsAsFactors=F, header=T)
CNV_P = CNV_P[,3:ncol(CNV_P)]

# Attach overlap
attach_version_02_overlap = F

if(attach_version_02_overlap) {
xx = PASS_CNV.overlap_sum_cov_02(CNV_I, CNV_P, db_levl = 0)
CNV_I = cbind(CNV_I, xx)
} ## ---- attach version 02
xx_L = PASS_CNV.overlap_sum_cov(CNV_I, CNV_P, db_levl = 0)
CNV_I = cbind(CNV_I, xx_L$xx)
xx_L = PASS_CNV.overlap_sum_cov(CNV_I, NULL, db_levl = 0)
CNV_I = cbind(CNV_I, ov_size_w=xx_L$xx[,1], ov_frac_w=xx_L$xx[,2], ov_n=xx_L$ov_n, ov_ID=xx_L$ov_ID)

if(attach_version_02_overlap) {
xx = PASS_CNV.overlap_sum_cov_02(CNV_P, CNV_I, db_levl = 0)
CNV_P = cbind(CNV_P, xx)
} ## ---- attach version 02
xx_L = PASS_CNV.overlap_sum_cov(CNV_P, CNV_I, db_levl = 0)
CNV_P = cbind(CNV_P, xx_L$xx)
xx_L = PASS_CNV.overlap_sum_cov(CNV_P, NULL, db_levl = 0)
CNV_P = cbind(CNV_P, ov_size_w=xx_L$xx[,1], ov_frac_w=xx_L$xx[,2], ov_n=xx_L$ov_n, ov_ID=xx_L$ov_ID)

## Output result
# fn = paste(Results_dir_local, "2018_07_Klaus_Annotations/v00_2018_06_15_Cuts.I..KSA.Annotated.ov.csv", sep="")
fn = paste(Results_dir_local, "2018_07_Klaus_Annotations/v00_2018_08_27_Cuts.I..KSA.Annotated.ov.csv", sep="")
write.csv(CNV_I, fn, row.names=F);  system(paste("open ", fn)) 
# fn = paste(Results_dir_local, "2018_07_Klaus_Annotations/v00_2018_06_15_Cuts.P..KSA.Annotated.ov.csv", sep="")
fn = paste(Results_dir_local, "2018_07_Klaus_Annotations/v00_2018_08_27_Cuts.P..KSA.Annotated.ov.csv", sep="")
write.csv(CNV_P, fn, row.names=F);  system(paste("open ", fn)) 

## order & Output result
uPos = 3e8*CNV_I[,"Chr"] + CNV_I[,"Start"]
o = order(uPos);    # CNV_I[o,"Chr"]
fn = paste(Results_dir_local, "2018_07_Klaus_Annotations/v00_2018_08_27_Cuts.I..KSA.Annotated.ov.Chr_order.csv", sep="")
write.csv(CNV_I[o,], fn, row.names=F);  system(paste("open ", fn)) 
uPos = 3e8*CNV_P[,"Chr"] + CNV_P[,"Start"]
o = order(uPos);    # CNV_P[o,"Chr"]
fn = paste(Results_dir_local, "2018_07_Klaus_Annotations/v00_2018_08_27_Cuts.P..KSA.Annotated.ov.Chr_order.csv", sep="")
write.csv(CNV_P[o,], fn, row.names=F);  system(paste("open ", fn)) 

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
# Continuation from      PASS_CNV_2017_collect_for_annotation_followup.R  small change
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
# ---- ---- output only confirmed
if(1) {
  ALG = "I"  #  I P
  if(ALG == "I") xx = CNV_I 
  if(ALG == "P") xx = CNV_P 
  
  uPos = 3e8*xx[,"Chr"] + xx[,"Start"]
  o    = order(uPos);    # xx[o,"Chr"]
  xx   = xx[o,]
  xxn  = xx
  ic   = c("C.Frec.Hap", "C.Sam.Hap", "P.Per.Hap", "P.Frec.Hap", "P.Sam.Hap", "C.Frec.All", "C.Sam.All", "P.Per.All", "P.Frec.All", "P.Sam.All")
  for(i in ic) {xxn[xxn[,i]==".",i] = 0;  xxn[,i] = as.numeric(xxn[,i])}
  s0      =  xxn[,"ov_frac"]   > 0; tNA(s0)
  xxn     =  xxn[s0,]
  xx      =  xx[ s0,]
  s.P.Hap = (xxn[,"P.Per.Hap"] < 30) | (xxn[,"P.Frec.Hap"] < 10)
  s.C.Hap =  xxn[,"C.Sam.Hap"] < 2
  s1      = s.P.Hap & s.C.Hap; tNA(s1)
  ip      = c(1,3,6, 8, 11:12, 24:38, 39:44)
  cbind(xx[,ip], s1, s.P.Hap, s.C.Hap)
  cbind(xx[s1,ip])
  fn = paste(Results_dir_local, "2018_07_Klaus_Annotations/v00_2018_08_27_Cuts.", ALG, "..KSA.Annotated.ov.Chr_order.rare.csv", sep="")
  write.csv(xx[s1,], fn, row.names=F);  system(paste("open ", fn)) 
} ## ---- used to produce output spreadsheets
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
# ---- ---- output all with selection flags
if(1) {
  ALG = "I"  #  I P
  if(ALG == "I") xx = CNV_I 
  if(ALG == "P") xx = CNV_P 
  
  uPos = 3e8*xx[,"Chr"] + xx[,"Start"]
  o    = order(uPos);    # xx[o,"Chr"]
  xx   = xx[o,]
  xxn  = xx
  ic   = c("C.Frec.Hap", "C.Sam.Hap", "P.Per.Hap", "P.Frec.Hap", "P.Sam.Hap", "C.Frec.All", "C.Sam.All", "P.Per.All", "P.Frec.All", "P.Sam.All")
  for(i in ic) {xxn[xxn[,i]==".",i] = 0;  xxn[,i] = as.numeric(xxn[,i])}
  s0      =  xxn[,"ov_frac"]   > 0; tNA(s0)
  s.P.Hap = (xxn[,"P.Per.Hap"] < 30) | (xxn[,"P.Frec.Hap"] < 10)
  s.C.Hap =  xxn[,"C.Sam.Hap"] < 2
  s1      = s0 & s.P.Hap & s.C.Hap; tNA(s1)
  ip      = c(1,3,6, 8, 11:12, 24:38, 39:44)
  cbind(xx[, ip], s1, s0, s.P.Hap, s.C.Hap)[1:5,]
  xx      = cbind(xx, s1, s0, s.P.Hap, s.C.Hap)
  if(0) {
    fn = paste(Results_dir_local, "2018_07_Klaus_Annotations/v00_2018_08_27_Cuts.", ALG, "..KSA.Annotated.ov.Chr_order.rare_flags.csv", sep="")
    write.csv(xx, fn, row.names=F);  system(paste("open ", fn)) 
  } ## ---- output spreadsheet
}
CNV_flg = xx ## CNVI_I/P with confirmation and common cnv flags
## ----------------------------------------------------------------------  make CNV vs exposure comparisons
if(0) {
  # this is table with CNVs found but not individually verified, checked with multiple software, checked for pop. frequency, or manually curated
  # ---- get names of chips with CNVs from CNV_s
  cuts               = CNV_cuts[[2]]
  cuts$Confidence[1] = 1200
  
  s                = PASS_CNV.select_CNV(CNV_s, cuts)
  x                = PASS_CNV.remove_AB(CNV_s[s, "Sample_Group"]);  length(x); 
  sample_names     = unique(x); length(sample_names);                                     # ---- names of individuals (chips) with CNVs
  # ---- create two tables, _all & _CNV
  s                = u_Summary_report[,"pass_QC"]  &  u_Summary_report[,"Collaborator.Participant"] != "NA12878"    # check for duplicates, none ## sort(tNA(u_Summary_report[s,"Collaborator.Participant"]))
  xx_all           = PASS_CNV.pheno_table_u_Summary_report(u_Summary_report, s, fn="")
  s                = s & u_Summary_report[,"Collaborator.Participant"] %in% sample_names
  xx_CNV           = PASS_CNV.pheno_table_u_Summary_report(u_Summary_report, s, fn="")
  # ---- align tables, _all & _CNV
  m = match(unlist(dimnames(xx_CNV)[1]), unlist(dimnames (xx_all)[1]))
  x = rep(0, dim(xx_all)[1])
  x[m] = xx_CNV[,"Total"]
  xx = cbind(xx_all, CNV=x)
  xx = cbind(xx, ratio=(xx[,"CNV"])/(xx[,"Total"]))
  xx
} ## ---- CNV discovery rate per phenotype  --> copy of old

if(1) {
  if(0) {
    cuts             = CNV_cuts[[1]]
    s                = PASS_CNV.select_CNV(CNV_flg, cuts)
    x                = PASS_CNV.remove_AB(CNV_flg[s, "Sample_Group"]);  length(x); 
    sample_names     = unique(x); length(sample_names);                                     # ---- names of individuals (chips) with CNVs
  } ## ---- no longer relevant
  if(0) {
    # ---- create two tables, _all & _CNV
    sample_names     = unique(CNV_flg[, "Sample_Group"]); length(sample_names);                                     # ---- names of individuals (chips) with CNVs
    s                = u_Summary_report[,"pass_QC"]  &  u_Summary_report[,"Collaborator.Participant"] != "NA12878"    # check for duplicates, none ## x = u_Summary_report[s,"Collaborator.Participant"]; length(x);  length(unique(x))
    xx_all           = PASS_CNV.pheno_table_u_Summary_report(u_Summary_report, s, fn="")
    s                = s & u_Summary_report[,"Collaborator.Participant"] %in% sample_names
    xx_CNV           = PASS_CNV.pheno_table_u_Summary_report(u_Summary_report, s, fn="")
  } ## ---- tables for presentation
  if(0) {
    m = match(unlist(dimnames(xx_CNV)[1]), unlist(dimnames (xx_all)[1]))
    x = rep(0, dim(xx_all)[1])
    x[m] = xx_CNV[,"Total"]
    xx = cbind(xx_all, CNV=x)
    xx = cbind(xx, ratio=(xx[,"CNV"])/(xx[,"Total"]))
    xx
  }  # ---- align tables, _all & _CNV and show discovery ratio per phenotype
  
  #    Collaborator.Participant             Demise.Type                          Specimen.Type Call.Rate Batch pass_QC E_Alc E_Smk
  # ---- select names of samples that carry CNVs
  s                = T
  s                = CNV_flg[,"s1"]
  sample_names     = unique(CNV_flg[s, "Sample_Group"]); length(sample_names);                                     # ---- names of individuals (chips) with CNVs
  # ---- select unique samples with Exposure
  u_T              = u_Summary_report
  s_u              = u_T[,"pass_QC"]  &  u_T[,"Collaborator.Participant"] != "NA12878";    # check for duplicates, none ## x = u_T[s,"Collaborator.Participant"]; length(x);  length(unique(x))
  u_T              = u_T[s_u,]
  if(0) {
    for(iPheno in c("E_Alc","E_Smk")) {
      s                = !is.na(u_T[,iPheno])
      print( tNA(u_T[  ,iPheno]) )
      print( tNA(u_T[s ,iPheno]) )
    }
  } ## ---- check exposure distribution
  ## ---- build table
  u_T             = u_T[!is.na(u_T[,"E_Alc"]),]        # both phenotypes are missing always together
  dAlc_EXP        = u_T[,"E_Alc"] > 2;  tNA(dAlc_EXP)
  dSmk_EXP        = u_T[,"E_Smk"] > 2;  tNA(dSmk_EXP)
  dual_EXP        = 2*as.integer(dAlc_EXP) + as.integer(dSmk_EXP)
  dEXP_m          = cbind(dAlc_EXP, dSmk_EXP, dual_EXP)
  dEXP_m_names    = c("Alchohol","Smoking","Dual")
  s_CNV            = u_T[,"Collaborator.Participant"] %in% sample_names;  tNA(s_CNV)
  for(i in 1:3) {
    cat(dEXP_m_names[i], "-------------------------------","\n")
    s_EXP            = dEXP_m[,i]
    print(table(s_CNV, s_EXP))
    if(i<3) print(fisher.test(table(s_CNV, s_EXP)))
  }
  
  if(0) {
    s_EXP            = dAlc_EXP
    cbind(s_CNV, s_EXP, u_T[,c("E_Alc","E_Smk")], x)[1:25,]
    table(s_CNV, s_EXP)
    fisher.test(table(s_CNV, s_EXP))
  } # -- different way to do table
  if(0) {
    s                = T
    s                = CNV_flg[,"s0"]
    sample_names     = unique(CNV_flg[s, "Sample_Group"]); length(sample_names);                                     # ---- names of individuals (chips) with CNVs
    # ---- select unique samples with Exposure
    u_T              = u_Summary_report
    s_u              = u_T[,"pass_QC"]  &  u_T[,"Collaborator.Participant"] != "NA12878";    # check for duplicates, none ## x = u_T[s,"Collaborator.Participant"]; length(x);  length(unique(x))
    u_T              = u_T[s_u,]
    s_CNV            = u_T[,"Collaborator.Participant"] %in% sample_names;  tNA(s_CNV)
    s                = is.na(u_T[,"E_Alc"])
    tNA(s_CNV)
    tNA(s_CNV[s])
  } ## ---- table usin s0 (only confirmed, not rare) instead of s1
  
} ## ---- CNV discovery rate per phenotype

# 1247 HapMap samples, 575/.4611 = 1247.018
# Looks like I should rerun tables with an AND for partial coverage CNVs

## ----------------------------------------------------------------------  2019_04_25 notes: cuts used for analysis
if(0) {
  # ---- ---- from PASS_CNV_2017_collect_for_annotation.R
  # ---- Illumina
  s = T
  s = s & xx[,       "Chr"] %in%   1:22;  tNA(s)
  s = s & xx[,     "Value"] !=    2;  tNA(s)
  s = s & xx[,"Confidence"] >   800;  tNA(s)
  # ---- PennCNV
  s = T
  s = s & xx[,       "Chr"] %in%   1:22;  tNA(s)
  s = s & xx[,     "Value"] !=    2;  tNA(s)
  s = s & xx[, "Confidence"] >   200;  tNA(s)
  # ---- 
  # For PennCNV the variable is conf=
  #   
  #   I have separae annotation from Klaus for CNV_I, CNV_P (Illumina, PennCNV).  
  # Overlap is determined by PASS_CNV_2017_collect_for_annotation_followup.R
  # 
  # CNV_I[1:10, c("Sample_Group", "Chr", "Start", "End","Size", "ov_size", "ov_frac", "ov_size_w", "ov_frac_w", "ov_n", "ov_ID")]
  # 
  # Annotation from Klaus:
  #   ic   = c("C.Frec.Hap", "C.Sam.Hap", "P.Per.Hap", "P.Frec.Hap", "P.Sam.Hap", "C.Frec.All", "C.Sam.All", "P.Per.All", "P.Frec.All", "P.Sam.All")
  # C/P          --> complete / partial
  # Freq/Sam/Per --> Frequency / # of samples / Percentage coverage
  #   Hap/All      --> Hapmap 1247 / All 10k samples
}
## ----------------------------------------------------------------------  


# ------------------------------------------------ Analysis

if(0) {

if(0) {
fn = paste(RESULTS_dir, "2018_07_Klaus_Annotations/v00_2018_06_15_Cuts.I..KSA.Annotated.csv", sep="")
xx = read.csv(fn,    stringsAsFactors=F, header=T)
s  = xx[,"Chr"] == xx[,"chr"]  &  xx[,"Start"] == xx[,"start"]  &  xx[,"End"] == xx[,"end"] 
tNA(s)
fn = paste(RESULTS_dir, "2018_07_Klaus_Annotations/v00_2018_06_15_Cuts.P..KSA.Annotated.csv", sep="")
xx = read.csv(fn,    stringsAsFactors=F, header=T)
s  = xx[,"Chr"] == xx[,"chr"]  &  xx[,"Start"] == xx[,"start"]  &  xx[,"End"] == xx[,"end"] 
tNA(s)
} ## ---- confirm alignment

s = CNV_I[,"Batch"] == "h2"
for(i in 24:26) { print(names(CNV_I)[i]);  print(tNA(CNV_I[s,i]))}

# cbind(1:ncol(CNV_I), names(CNV_I))
xx = CNV_I
s  = xx[,"ov_frac"] > 0
summary(xx[s,"ov_frac"])
quantile(xx[s,"ov_frac"], prob=0.1*0:10)



}
