example.pfb
Name	Chr	Position	PFB
rs1000061	11	88115425	0.388376383763838
rs1000071	11	4557878	0.870609981515712
rs1000079	3	152108465	0.877551020408163
rs1000121	20	16171957	0.172897196261682
rs1000122	20	61192853	0.144194756554307
rs1000131	3	1551741	0.109778597785978
rs1000198	3	120596510	0.475092250922509

-------------------------------------------------------------------
head(bed_coord)
  Chr pos_b0 Position       Name
1   1  82153    82154  rs4477212
2   1 752565   752566  rs3094315
3   1 752720   752721  rs3131972
4   1 762319   762320 exm2268640
5   1 768447   768448 rs12562034
6   1 776545   776546 rs12124819

> tNA(bed_coord[,"Chr"])
x
    1    10    11    12    13    14    15    16    17    18    19     2    20    21    22     3     4     5     6     7     8     9     M     X     Y  <NA> 
83020 48003 51672 47332 31771 30669 29579 32782 33010 25183 29676 74080 24594 12910 15411 61096 50050 52854 62629 48538 45316 42528   208 23500  1387     0 

-------------------------------------------------------------------
head(pfb_Yale)
         Name Chr  Position   PFB
1 kgp23671479   7 110336274 0.972
2  kgp4072255  10 130192301 0.397
3   rs4935417  10  55055218 0.211
4  kgp5626194  16  56341037 0.328
5  rs10806671   6 170270028 0.553
6 kgp10678636   6 105777953 0.154
> tNA(pfb_Yale[,"Chr"])
x
     0      1     10     11     12     13     14     15     16     17     18     19      2     20     21     22      3      4      5      6      7      8      9 
 17468 371562 226039 222463 225609 162299 151756 145202 157042 140641 133825 104213 373936 115668  64706  67661 326591 293615 279109 320907 250962 240336 208243 
    MT      X     XY      Y   <NA> 
   288 129344   1039   3071      0 

-------------------------------------------------------------------






# ------------------------------------------------------------------- This is for `2018_05` PPT prep

Consolidated in "PASS_CNV_2017_both_analysis.R"

# ------------------------------------------------------------------- Develop code

# ---- ---- 

Pheno_fn          = "~/Documents/shared/proj/PASS_CNV/Analysis/Genetics List 1st and 2nd Batch 9.26.17.csv"  ## Phenotypes from Jacob, 1st & 2nd batch
Pheno          = read.csv(Pheno_fn,          stringsAsFactors=F, header=T, skip= 0)

# ---- ---- 
h1$CNV_s          = PASS_CNV.attach_pheno_info(h1$CNV_s,          h2$Pheno, "Sample_Group",             "Participant.ID", c("Demise.Type", "Specimen.Type"))
# ---- ---- Create unified tables
u_Summary_report = PASS_CNV.create_u_Summary_report(h1, h2)
CNV_s            = rbind(h1$CNV_s, h2$CNV_s); # dim(CNV_s); dim(h1$CNV_s); dim(h2$CNV_s); 


# ------------------------------------------------------------------- This is for `2018_05` PPT prep
# Comput PennCNV in a small subset

if(0) {
module load dev/perl/5.8.9
cd /groups/markianos/programs/PennCNV/penncnv/example

DATA=/groups/markianos/Data/Chip_Data_Human/2016_PASS_CNV/
ls $DATA
PFB=$DATA/PennCNV_files/KM_InfiniumOmniExpressExome-8v1-3_A_hg19.pfb
GCM=$DATA/PennCNV_files/KM_InfiniumOmniExpressExome-8v1-3_A_hg19.gcmodel
GEN=$DATA/Broad_PASS_processed/200557030023_R02C01.gtc.txt
OUTD=/groups/markianos/km/Chips_human/PASS_CNV_2016/zz_OUT
LOG=$OUTD/test.log
OUT=$OUTD/test.rawcnv

# detect_cnv.pl -test -hmm example.hmm -pfb $PFB -gcmodel $GCM -conf -log $LOG -out $OUT $GEN
# detect_cnv.pl -test -hmm example.hmm -pfb $PFB               -conf -log $LOG -out $OUT $GEN
# *********** -----> this worked, use their gcmodel file:   detect_cnv.pl -test -hmm example.hmm -pfb $PFB -gcmodel ../lib/hh550.hg18.gcmodel -conf -log $LOG -out $OUT $GEN

# -----------------------------------> actual run
cd $OUTD/../
LOG=$OUTD/KM_2017_09_10_with_gcmodel.log
OUT=$OUTD/KM_2017_09_10_with_gcmodel.rawcnv
HMM=/groups/markianos/programs/PennCNV/penncnv/example/example.hmm

ls $DATA/Broad_PASS_processed/*.gtc.txt > inputlist_2017_09_10.txt
#    detect_cnv.pl -test -hmm  example.hmm -pfb $PFB -gcmodel $GCM -conf -log $LOG -out $OUT $GEN
time detect_cnv.pl -test -hmm  $HMM        -pfb $PFB -gcmodel $GCM -conf -log $LOG -out $OUT -list inputlist_2017_09_10.txt &

# detect_cnv.pl -test -hmm example.hmm -pfb example.pfb -log ex3.log -out ex3.rawcnv -gcmodel ../lib/hh550.hg18.gcmodel -conf -list inputlist
} ## ---- original

cd /n/groups/markianos/km/Chips_human/PASS_CNV_2016

I repeated p1st using sbatch at O2. 
sbatch_00.sh --> test
sbatch_01.sh --> run again


 
if(0) {
#!/bin/bash
#SBATCH -c 1                               # 1 core
#SBATCH -t 0-00:05                         # Runtime of 5 minutes, in D-HH:MM format
#SBATCH -p short                           # Run in short partition
#SBATCH -o hostname_sinfo_%j.out           # File to which STDOUT + STDERR will be written, including job ID in filename
#SBATCH --mail-type=ALL                    # ALL email notification type
#SBATCH --mail-user=Kyriacos.Markianos@childrens.harvard.edu # Email to which notifications will be sent
 
srun hostname
srun sinfo
} ## ---- O2 example, the basis for sbatch_00.sh
if(0) {
sbatch sbatch_00.sh
squeue -u km111

sbatch sbatch_01.sh

ls -lh zz_OUT/2018_05_03/;  ls -lh
rm -f *.out;  rm -f *.err
cat zz_OUT/KM_2017_09_10_with_gcmodel.log | grep 'LRR values for all autosome markers from' | wc -l
cat zz_OUT/2018_05_04_all_2nd/test.log    | grep 'LRR values for all autosome markers from' | wc -l

# -----------------------------------> actual run
cd $OUTD/../
LOG=$OUTD/KM_2017_09_10_with_gcmodel.log
OUT=$OUTD/KM_2017_09_10_with_gcmodel.rawcnv
HMM=/n/groups/markianos/programs/PennCNV/penncnv/example/example.hmm

ls $DATA/Broad_PASS_processed/*.gtc.txt > inputlist_2017_09_10.txt
#    detect_cnv.pl -test -hmm  example.hmm -pfb $PFB -gcmodel $GCM -conf -log $LOG -out $OUT $GEN
time detect_cnv.pl -test -hmm  $HMM        -pfb $PFB -gcmodel $GCM -conf -log $LOG -out $OUT -list inputlist_2017_09_10.txt &

# detect_cnv.pl -test -hmm example.hmm -pfb example.pfb -log ex3.log -out ex3.rawcnv -gcmodel ../lib/hh550.hg18.gcmodel -conf -list inputlist

} ## ---- sbatch_01.sh
 
Updated documentation in "PASS_CNV_2016_Analysis_Notes_and_project_notes.tex".  
Preperation steps for pfb, gcmodel, gtx.txt in "tex" file.

PASS_CNV_2017.PennCnv_pfb_prep.R (pfb & gcmodel) runs local and it is chip specific.
Make _p1st & _p2nd files, p1st identical to original, and process locally?  No, just double in/out files and make process_batch = c(1,2) to direct processing.

PASS_CNV_2017.PennCnv_gen_prep.R (gtx.txt, file local but run in O2) needs updating too. Use process_batch = c(1,2) too.

Make new directories for all outputs.


# ------------------------------------------------------------------- `2018_05_16`

missing pfb, rsids, bed files

revisit bad assays file (InfiniumOmniExpressExome-8v1-4_A1.1.2.bad_assays.csv)

/Volumes/hd0/Users/kmarkianos/Documents/shared/proj/PASS_CNV/Data/Broad_2nd_hyb/Holm_PASS_OmniExpressExomeArrays_50_samples_8-2017/Holm_PASS_OmniExpressExomeArrays_50_samples_8-2017.vcf

fn_vcf = "~/Documents/shared/proj/PASS_CNV/Data/Broad_2nd_hyb/Holm_PASS_OmniExpressExomeArrays_50_samples_8-2017/Holm_PASS_OmniExpressExomeArrays_50_samples_8-2017.vcf"

vcf = as.data.frame(fread(fn,       header=T, skip=11, nrows=-1))

Google: Converting Illumina Raw Genotype Data Into ascii Format
Export final report with BeadStudio --> https://www.biostars.org/p/2240/

my_gen_1b = my_gen
dim(my_gen_1b); dim(my_gen);
m = my_gen[,"SNP Name"] %in% my_gen_1b[,"SNP Name"] 
x = my_gen[!m,"SNP Name"]
x = my_gen[!m,"Chromosome"];  tNA(x)
bed_coord_1b = as.data.frame(fread("~/Documents/shared/proj/PASS_CNV/Analysis/Infinium_PASS_CNV_PennCNV_related_downloads/InfiniumOmniExpressExome-8v1-3_A/InfiniumOmniExpressExome-8v1-3_A.bed",   header=F, skip=1));
names(bed_coord_1b) = c("Chr", "pos_b0", "Position", "Name")
m = match(bed_coord_1b[, "Name"], bed_coord[, "Name"])
xx = cbind(bed_coord_1b, bed_coord[m,])

# ------------------------------------------------------------------- `2018_05_25`, `2018_06_01`

- finished pfb & GC
- move to O2
- prepare gtc files


cd /n/groups/markianos/km/Chips_human/
cd /n/groups/markianos/km/Chips_human/PASS_CNV_2016

process once more 1st and put in tmp
process 2nd
update local sbatch_02_2nd.sh and run ... where is the slurm command?

# DATA_dir = "/Volumes/hd0/Users/kmarkianos/Documents/shared/proj/PASS_CNV/Data/Broad_2nd_hyb/export_gtc/gtc"
# ReadStatT_L_names   = list.files(path=DATA_dir, recursive = T, full.names = F)

cd /n/groups/markianos/km/Chips_human/PASS_CNV_2016
sbatch sbatch_02_2nd.sh

Checked sample Holm_OmniEx_Exome_RW_01_A04_S002-BISH-06078_201486520003_R01C01 (Sample 7 in GenomeStudio), exm-rs1183201 and it is NaN for BalleleF & logR in GenomeStudio despite being a called SNP in 49/51 samples.

fn           = RUN_dir_L[1]; 
fn_gen_out   = paste(gen_out_DIR, ReadStatT_L_names[ifn], sep="/")
cat(ifn, fn_gen_out, "\n")
my_gen       = as.data.frame(fread(fn,           header=T, skip=11, nrows=-1)); ## colClasses = c("character", "character", "integer", "character", "numeric", "numeric")
x            = paste(my_gen[, "Allele1 - AB"], my_gen[, "Allele2 - AB"], sep=""); x[x=="--"] = "NC"; ## tNA(x);  tNA(example_gen[,"99HI0700A.GType"])
xx           = cbind(my_gen[, c("SNP Name", "Chromosome", "Position")], x.GType=x, my_gen[,c("Log R Ratio Illumina","bAllele Freq")])
names(xx)    = c("Name", "Chr", "Position", "x.GType", "x.Log R Ratio", "x.B Allele Freq")

xx = xx[xx[,"Chr"]==1,]
o = order(xx[,"Position"])
xx = xx[o,]
s = which(is.na(xx[,"x.B Allele Freq"]))
x = c(s, s-1, s+1)
x = sort(x)
xx[x[1:18],]

# ------------------------------------------------------------------- `2018_06_14`

Finished running PennCNV in O2


CNV_table_for_Klaus_fn = "~/Documents/shared/proj/PASS_CNV/PASS_CNV_tmp/Results/v00_all_CNV_table_select_CNVa_2017_02_10_v01_for_Kyriacos.csv"
CNV_table_for_Klaus    = read.csv(CNV_table_for_Klaus_fn,      stringsAsFactors=F, header=T, skip= 0)
names(CNV_table_for_Klaus)[3] = "Sample_Group"

s = CNV_table_s[,"Confidence"] > 1000

x.I = paste(CNV_table_s[s       ,"Sample_Group"], CNV_table_s[s,        "Size"])
x.K = paste(CNV_table_for_Klaus[,"Sample_Group"], CNV_table_for_Klaus[, "Size"])

x = match(x.K, x.I)
unique(diff(x)) 

So selection of Illumina only with a strict (CNV_table_s[,"Confidence"] > 1000) selection.

# ------------------------------------------------------------------- `2018_06_15`
h1$CNV_s          = PASS_CNV.attach_pheno_info(h1$CNV_s,          h2$Pheno, "Sample_Group",             "Participant.ID", c("Demise.Type", "Specimen.Type"))
xx                = PASS_CNV.attach_pheno_info(CNV_table_Pcnv,    h2$Pheno, "Sample_Group",             "Participant.ID", c("Demise.Type", "Specimen.Type"))
xx                = PASS_CNV.attach_pheno_info(CNV_table_Pcnv,    h2$Pheno, "Sample_Group",             "Participant.ID", c("Demise.Type", "Specimen.Type"))

o  = order(CNV_table_Pcnv[,"Size"], decreasing=T)
CNV_table_Pcnv = CNV_table_Pcnv[o,]

xh1 = xx[, "filename"]
xh2 = xx[, "filename"]

xh1[1:2]
xh2[1:2]

sub("\\/groups\\/markianos\\/Data\\/Chip_Data_Human\\/2016_PASS_CNV\\/\\/Broad_PASS_processed\\/(.*)\\.gtc\\.txt", "\\1", xh1[1:2])
sub("\\/groups\\/markianos\\/Data\\/Chip_Data_Human\\/2016_PASS_CNV\\/\\/Broad_PASS_processed\\/(.*)\\.gtc\\.txt", "\\1", xh2[1:2])

sub(".*\\/Broad_PASS_processed.*\\/(.*)\\.gtc\\.txt", "\\1", xh1[1:2])
sub(".*\\/Broad_PASS_processed.*\\/(.*)\\.gtc\\.txt", "\\1", xh2[1:2])

# ------------------------------------------------------------------- `2018_06_17 & 18`

Sent CNVs to Klaus.  Need QuantiSNP and and detailed review of current result.

See comments in "PASS_CNV_2017_collect_for_annotation.R", I run twice and save catalogs 
in h1 & h2$CNV_GStudio_s, h2$CNV_PennCNV_s 

Used save.image() to save session to PASS_CNV_2017_collect_for_annotation.2018_06_18.RData.
Just load("PASS_CNV_2017_collect_for_annotation.2018_06_18.RData")

# ------------------------------------------------------------------- `2018_06_27`

from PennCNV documentation:   Note: By default only autosome CNVs will be detected, but the --chrx argument can be used in the above command to generate CNV calls for chromosome X only. The CNV calling for chrX is slightly different from that of autosomes. It is highly recommended to use the --sexfile argument to supply gender annotation for all genotyped samples.

# ------------------------------------------------------------------- `2018_08_06`

cbind(names(h1$CNV_GStudio_s), names(h2$CNV_GStudio_s))

x = PASS_CNV.overlap_sum(xx1, xx2)

xx1         = h2$CNV_GStudio_s
xx2         = h2$CNV_PennCNV_s
overlap_thr = 0.1

PASS_CNV.overlap_sum <- function(xx1, xx2, overlap_thr) {
ov = rep(0, nrow(xx1))
for(ic1 in 1:nrow(xx1)) {
  xx = xx2[xx1[ic1, "Chr"] == xx2[, "Chr"], ]
  if(nrow(xx) == 0) next
  for(ic in 1:nrow(xx)) {
    min_length = min(xx1[ic1,"End"] - xx1[ic1,"Start"], xx[ic,"End"] - xx[ic,"Start"])
    if((xx1[ic1,"End"] - xx [ic, "Start"]) / min_length > overlap_thr) { ov[ic1] = 1; next }
    if((xx[ ic ,"End"] - xx1[ic1,"Start"]) / min_length > overlap_thr) { ov[ic1] = 1; next }
  }
}
return(ov)
}
# ------------------------------------------------------------------- overlap I, P

xx1         = h2$CNV_GStudio_s
xx2         = h2$CNV_PennCNV_s
PASS_CNV.overlap_sum_cov <- function(xx1, xx2, db_levl = 2) {
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
return(cbind(ov_size=ov_size, ov_frac))
}
xx = PASS_CNV.overlap_sum_cov(xx1, xx2, db_levl = 2)
xx

ip = c(1,3:6, 8, 11:12)
cbind(xx1[,ip], xx)

s = xx[,"ov_frac"] == 0.
s = paste(xx2[,"Sample_Group"], xx2[,"Chr"], sep="_") %in% paste(xx1[s,"Sample_Group"], xx1[s,"Chr"], sep="_") 
cbind(xx1[,ip], xx)
xx2[s,ip]

# ------------------------------------------------------------------- overlap P, I

xx = PASS_CNV.overlap_sum_cov(xx2, xx1, db_levl = 2)
xx

# ------------------------------------------------------------------- overlap final
xx = PASS_CNV.overlap_sum_cov(h2$CNV_GStudio_s, h2$CNV_PennCNV_s, db_levl = 0)
h2$CNV_GStudio_s_ov = cbind(h2$CNV_GStudio_s, xx)
xx = PASS_CNV.overlap_sum_cov(h2$CNV_PennCNV_s, h2$CNV_GStudio_s, db_levl = 0)
h2$CNV_PennCNV_s_ov = cbind(h2$CNV_PennCNV_s, xx)

ip = c(1,3:6, 8, 11:12)
ip = c(1,3:6, 8, 11:19)
h2$CNV_GStudio_s_ov[,ip]


# ------------------------------------------------------------------- `2018_08_10` final iteration

# ------------------------------------------------------------------- final iteration
read joined h1/h2 from Klaus
comput overlaps and produce CNV_GStudio_s_ov, CNV_PennCNV_s_ov

Alternative, join h1/h2 and process.  Then read annotations from Klaus.... not much sense

~/Documents/shared/proj/PASS_CNV/PASS_CNV_tmp/Results/2018_07_Klaus_Annotations/v00_2018_06_15_Cuts.I..KSA.Annotated.csv

fn = paste(RESULTS_dir, "2018_07_Klaus_Annotations/v00_2018_06_15_Cuts.I..KSA.Annotated.csv", sep="")
CNV_I = read.csv(fn,    stringsAsFactors=F, header=T)
CNV_I = CNV_I[,3:ncol(CNV_I)]
fn = paste(RESULTS_dir, "2018_07_Klaus_Annotations/v00_2018_06_15_Cuts.P..KSA.Annotated.csv", sep="")
CNV_P = read.csv(fn,    stringsAsFactors=F, header=T)
CNV_P = CNV_P[,3:ncol(CNV_P)]

xx = PASS_CNV.overlap_sum_cov(CNV_I, CNV_P, db_levl = 0)
CNV_I = cbind(CNV_I, xx)

# ------------------------------------------------------------------- post included in _followup
if(0) {
s = CNV_I[,"Batch"] == "h2"
ip = c(1,3:6, 8, 11:12)
ip = c(1,3:6, 8, 11:12, 22, 24:40)
ip = c(1,3:6, 8, 11:12,     24:40)
CNV_I[s,ip]
} ## ---- visualize h2 overlap

if(0) {
cbind(1:ncol(CNV_I), names(CNV_I))
      [,1] [,2]           
 [1,] "1"  "Sample_Group" 
 [2,] "2"  "SampleID"     
 [3,] "3"  "Chr"          
 [4,] "4"  "Start"        
 [5,] "5"  "End"          
 [6,] "6"  "Size"         
 [7,] "7"  "numsnp"       
 [8,] "8"  "Value"        
 [9,] "9"  "startsnp"     
[10,] "10" "endsnp"       
[11,] "11" "Confidence"   
[12,] "12" "Batch"        
[13,] "13" "Algr"         
[14,] "14" "Demise.Type"  
[15,] "15" "Specimen.Type"
[16,] "16" "E_Alc"        
[17,] "17" "E_Smk"        
[18,] "18" "X.."          
[19,] "19" "chr"          
[20,] "20" "start"        
[21,] "21" "end"          
[22,] "22" "UCSC"         
[23,] "23" "X..Hap"       
[24,] "24" "C.Frec.Hap"   
[25,] "25" "C.Sam.Hap"    
[26,] "26" "C.cnvs.Hap"   
[27,] "27" "P.Per.Hap"    
[28,] "28" "P.Frec.Hap"   
[29,] "29" "P.Sam.Hap"    
[30,] "30" "P.cnvs.Hap"   
[31,] "31" "X..All"       
[32,] "32" "C.Frec.All"   
[33,] "33" "C.Sam.All"    
[34,] "34" "C.cnvs.All"   
[35,] "35" "P.Per.All"    
[36,] "36" "P.Frec.All"   
[37,] "37" "P.Sam.All"    
[38,] "38" "P.cnvs.All"   
[39,] "39" "ov_size"      
[40,] "40" "ov_frac"      
} ## ---- save column names

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

PASS_CNV.overlap_sum_cov_00 <- function(xx1, xx2, db_levl = 2) {
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
return(cbind(ov_size=ov_size, ov_frac))
}

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

PASS_CNV_Interval_index_list <- function(Sample_eq, x, xx) {
  # Sample_eq flag, if T require Sample == xx[,"Sample"] corresponding to compare against external database
  # x  -> 4 element vector    [Sample, Chr, Start, Stop]
  # xx -> 4 column data.frame [Sample, Chr, Start, Stop]
  IL      = list()
  s_sc    =  x["Chr"] == xx[, "Chr"]
  if(Sample_eq)  s_sc = (x[1] == xx[, 1]) & s_sc
  IL$s_sc    = s_sc
  IL$sI      = findInterval(xx[,"Start"], x[c("Start","End")])
  IL$eI      = findInterval(xx[,  "End"], x[c("Start","End")])
  IL$s_ov    = !((sI+eI == 0) | (sI+eI == 4)) ## some overlap
  IL$ov_flag = sum(IL$s_sc & IL$s_ov) > 0
  return(IL)
}

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

PASS_CNV.overlap_sum_cov_01 <- function(Sample_eq, xx1, xx2, db_levl = 2) {
ov_size = ov_frac = rep(0, nrow(xx1))
for(ic1 in 1:nrow(xx1)) {
  cnvI    = xx1[ic1, c("Sample_Group", "Chr", "Start","End")]
  IL      = PASS_CNV_Interval_index_list(Sample_eq, cnvI, xx2)
  if(!IL$ov_flag) next
  if(db_levl > 0) cat("------------------------------------------------------ ic1=", ic1, "\n")

  s = sum(IL$sI==0 & IL$eI==2 & IL$s_sc);  if(sum(s)) { cnvI["End"] == cnvI["Start"];   next;} ## ---- full overlap
  if(db_levl > 1) cat("all out=", sum(s), "\n")

  s = sum(IL$sI==0 & IL$eI==1 & IL$s_sc);  if(s) { cnvI["Start"] = max(xx2[  s,"End"]) }
  if(db_levl > 1) cat("left   =", sum(s), "\n")
  
  s = sum(IL$sI==1 & IL$eI==2 & IL$s_sc);  if(s) { cnvI[  "End"] = max(c(cnvI["Start"], min(xx2[s,"Start"]))) }
  if(db_levl > 1) cat("right  =", sum(s), "\n")


  if(s) IL = PASS_CNV_Interval_index_list(Sample_eq, cnvI, xx2)
  s = sI==1 & eI==1 & s_sc;               ov_size[ic1] = ov_size[ic1] +  sum(xx2[  s,"End"] - xx2[  s,"Start"], na.rm=T)
  if(db_levl > 1) cat("all in =", sum(s), "\n")
  
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
return(cbind(ov_size=ov_size, ov_frac))
}

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
xx1         = h2$CNV_GStudio_s
xx2         = h2$CNV_PennCNV_s
ic1         = 1
s_sc        = xx1[ic1, "Chr"] == xx2[, "Chr"]
s_sc        = (xx1[ic1, "Sample_Group"] == xx2[, "Sample_Group"]) & s_sc
x           = unlist(xx1[ic1,  c("Start","End")])
xx          =        xx2[s_sc, c("Start","End")]


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
ss_interval = c(1,4)
xx = cbind(-1:1, 0:2)
xx = rbind(xx[1:3,], xx[3,])
xx = rbind(xx      , cbind(3:5, 4:6))
# xx = rbind(xx,cbind(1,3))
# xx = cbind(1,4)
xx = cbind(-1:0, 0:1)
ss_table = xx  


# ss_table = cbind(-5,2)
ov_L = Find_overlap_regions(ss_interval,ss_table, db_levl=2)
ov_L$ov_size; ov_L$ov_frac; ov_L$ss_interval_flag ;cbind(ss_table, ov_L$ss_table_flag)

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

s = CNV_I[,"ov_size_02"] != CNV_I[,"ov_size"]
CNV_I[s,]

s = CNV_P[,"ov_size_02"] != CNV_P[,"ov_size"]
CNV_P[s,]

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

uPos = 3e8*CNV_I[,"Chr"] + CNV_I[,"Start"]
o    = order(uPos);    # CNV_I[o,"Chr"]
ip   = c(1,3,6, 8, 11:12, 24:38, 39:44)
xx = CNV_I  # [o,]

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


if(0) {
s = CNV_I[,"ov_frac"] > 0; tNA(s)
CNV_I[s,]
ip = c(1,3:6, 8, 11:12)
ip = c(1,3:6, 8, 11:12, 24:38)
ip = c(1,3:6, 8, 11:12, 24:38, 39:44)
ip = c(1,3,6, 8, 11:12, 24:38, 39:44)
CNV_I[s,ip]
} ## ---- select for display

# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
# Continuation from      PASS_CNV_2017_collect_for_annotation_followup.R
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

1247 HapMap samples, 575/.4611 = 1247.018
Looks like I should rerun tables with an AND for partial coverage CNVs

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
For PennCNV the variable is conf=

I have separae annotation from Klaus for CNV_I, CNV_P (Illumina, PennCNV).  
Overlap is determined by PASS_CNV_2017_collect_for_annotation_followup.R

CNV_I[1:10, c("Sample_Group", "Chr", "Start", "End","Size", "ov_size", "ov_frac", "ov_size_w", "ov_frac_w", "ov_n", "ov_ID")]

Annotation from Klaus:
ic   = c("C.Frec.Hap", "C.Sam.Hap", "P.Per.Hap", "P.Frec.Hap", "P.Sam.Hap", "C.Frec.All", "C.Sam.All", "P.Per.All", "P.Frec.All", "P.Sam.All")
C/P          --> complete / partial
Freq/Sam/Per --> Frequency / # of samples / Percentage coverage
Hap/All      --> Hapmap 1247 / All 10k samples
}
## ----------------------------------------------------------------------  











