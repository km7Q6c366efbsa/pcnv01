# Create tables for analysis
## ----
RESULTS_dir = "~/Documents/shared/proj/PASS_CNV/PASS_CNV_tmp/Results/"
PLOT_dir    = paste("","zzPDF/", sep="");
save_object_list = c("save_object_list", "HOME_dir", "RESULTS_dir", "PLOT_dir", "h1", "h2")
save_object_list = c(save_object_list, "tNA", "PASS_CNV.remove_AB",        "PASS_CNV.attach_pheno_info")
save_object_list = c(save_object_list, "PASS_CNV.create_u_Summary_report", "PASS_CNV.pheno_table_u_Summary_report")
save_object_list = c(save_object_list, "PASS_CNV.select_CNV",              "PASS_CNV.define_CNV_cuts",            "CNV_cuts")
save_object_list = c(save_object_list, "u_Summary_report", "CNV_s")
## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- Joint analysis
source("PASS_CNV_2017_p1st.R", echo=F)
rm(list=setdiff(ls(), save_object_list))
source("PASS_CNV_2017_p2nd.R", echo=F)
rm(list=setdiff(ls(), save_object_list))

# ls(); names(h1)
# [1] "SSheet"         "Summary_report" "Pheno"          "CNV_table"      "s"            "CNV_s"            

## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- Joint analysis
# ---- ---- output the Summary_report for samples shared between runs
if(0) {
x1 = PASS_CNV.remove_AB( h1$Summary_report[, "Collaborator.Participant"] )
x2 = PASS_CNV.remove_AB( h2$Summary_report[, "Collaborator.Participant"] )
x  = intersect(x1, x2)
# [1] "S002-BISH-01139" "S002-BISH-01604" "S002-BISH-02035" "S002-BISH-00633" "S002-BISH-02228" "S002-BISH-03722" "NA12878"        

s = x2 %in% x;    fn = "~/Documents/shared/proj/PASS_CNV/PASS_CNV_tmp/Results/xx_Summary_report_repeats.csv";   write.csv(h2$Summary_report[s, ], fn);  system(paste("open ", fn))
s = x1 %in% x;    fn = "~/Documents/shared/proj/PASS_CNV/PASS_CNV_tmp/Results/xx_Summary_report_repeats_h1_performance.csv";   write.csv(h1$Summary_report[s, ], fn);  system(paste("open ", fn))
# saved repeats in xx_Summary_report_repeats_2018_03_05.xlsx
}
# ---- ---- identify failed chips and remove from CNV_s
s = h1$Summary_report[, "Call.Rate"] < 98.;  h1$failed_chips = h1$Summary_report[s, "Collaborator.Participant"]
s = h2$Summary_report[, "AutoCall"]  < 98.;  h2$failed_chips = h2$Summary_report[s, "Collaborator.Participant"]

s = h1$CNV_s[,"Sample_Group"] %in% h1$failed_chips; tNA(s);  h1$CNV_s = h1$CNV_s[!s, ]
s = h2$CNV_s[,"Sample_Group"] %in% h2$failed_chips; tNA(s);  h2$CNV_s = h2$CNV_s[!s, ]

# ---- ---- attach phenotype to CNV_s and other tables
# use only h2$Pheno that contains all cases.  Repeats appear twice in h2$Pheno
if(0) {
i = which(duplicated(h2$Pheno[,"Participant.ID"]))
print( h2$Pheno[i, 1:4] )
} ## ---- show that the 6 repeats appear twice in the second pheno file

if(1) {
if(0) {
x = PASS_CNV.remove_AB( h1$CNV_s[, "Sample_Group"] )
i  = match(x, PASS_CNV.remove_AB(h2$Pheno[,"Participant.ID"]))
cbind(h1$CNV_s[,"Sample_Group"], x, h2$Pheno[i,])[1:10, 1:5]
h1$CNV_s[,"Batch"] = "h1"
h1$CNV_s[,"pheno"] = h2$Pheno[i, "Demise.Type"]

x = PASS_CNV.remove_AB( h2$CNV_s[, "Sample_Group"] )
i  = match(x, PASS_CNV.remove_AB(h2$Pheno[,"Participant.ID"]))
cbind(h2$CNV_s[,"Sample_Group"], x, h2$Pheno[i,])[1:10, 1:5]
h2$CNV_s[,"Batch"] = "h2"
h2$CNV_s[,"pheno"] = h2$Pheno[i, "Demise.Type"]
} ## ---- previous version of code
h1$CNV_s[,"Batch"] = "h1"
h2$CNV_s[,"Batch"] = "h2"

# ---- Rename erroneously entered chip in manifest
i = grep("S002-SANF-01782", h2$CNV_s[ ,                    "Sample_Group"]);  h2$CNV_s[i,                     "Sample_Group"] = "S002-BISH-01782"
i = grep("S002-SANF-01782", h2$SSheet[,                    "Sample_Group"]);  h2$SSheet[i,                    "Sample_Group"] = "S002-BISH-01782"
i = grep("S002-SANF-01782", h2$Summary_report[,"Collaborator.Participant"]);  h2$Summary_report[i,"Collaborator.Participant"] = "S002-BISH-01782"

# ---- Attach phenotype
h1$CNV_s          = PASS_CNV.attach_pheno_info(h1$CNV_s,          h2$Pheno, "Sample_Group",             "Participant.ID", c("Demise.Type", "Specimen.Type"))
h2$CNV_s          = PASS_CNV.attach_pheno_info(h2$CNV_s,          h2$Pheno, "Sample_Group",             "Participant.ID", c("Demise.Type", "Specimen.Type"))
h1$SSheet         = PASS_CNV.attach_pheno_info(h1$SSheet,         h2$Pheno, "Sample_Group",             "Participant.ID", c("Demise.Type", "Specimen.Type"))
h2$SSheet         = PASS_CNV.attach_pheno_info(h2$SSheet,         h2$Pheno, "Sample_Group",             "Participant.ID", c("Demise.Type", "Specimen.Type"))
h1$Summary_report = PASS_CNV.attach_pheno_info(h1$Summary_report, h2$Pheno, "Collaborator.Participant", "Participant.ID", c("Demise.Type", "Specimen.Type"))
h2$Summary_report = PASS_CNV.attach_pheno_info(h2$Summary_report, h2$Pheno, "Collaborator.Participant", "Participant.ID", c("Demise.Type", "Specimen.Type"))

# s = unique(c(which(is.na(h2$SSheet[,"Demise.Type"])), which(is.na(h2$SSheet[,"Specimen.Type"]))))
# s = unique(c(which(is.na(h2$Summary_report[,"Demise.Type"])), which(is.na(h2$Summary_report[,"Specimen.Type"]))))
# h2$Summary_report[s,]
# s = unique(c(which(is.na(h1$Summary_report[,"Demise.Type"])), which(is.na(h1$Summary_report[,"Specimen.Type"]))))
# h1$Summary_report[s,]
} ## ---- attach phenotype, rename one chip

# ---- ---- Create unified tables
u_Summary_report = PASS_CNV.create_u_Summary_report(h1, h2)
CNV_s            = rbind(h1$CNV_s, h2$CNV_s); # dim(CNV_s); dim(h1$CNV_s); dim(h2$CNV_s); 

# ---- ---- Attach Exposure
if(1) {
Exposure_fn    = "~/Documents/shared/proj/PASS_CNV/Data/2017_09_pheno/PASS Trajectory Groups.csv"  ## Exposure from Jacob, 1st & 2nd batch, 2017_09
Exposure       = read.csv(Exposure_fn,          stringsAsFactors=F, header=T, skip= 0)
x              = length(Exposure);  
names(Exposure)[(x-1):x] = c("E_Alc", "E_Smk")

u_Summary_report = PASS_CNV.attach_pheno_info(u_Summary_report, Exposure, "Collaborator.Participant", "patid", c("E_Alc", "E_Smk"))
CNV_s            = PASS_CNV.attach_pheno_info(CNV_s,            Exposure,             "Sample_Group", "patid", c("E_Alc", "E_Smk"))

if(0) {
s = is.na(CNV_s[,"E_Alc"])
CNV_s[s,]  ## two HapMap and 3 RETR (retrospective)
s = is.na(u_Summary_report[,"E_Alc"])
u_Summary_report[s,]  ## two HapMap and 3 RETR (retrospective)
} ## ---- show missing exposure: HapMap & retrospective (RETR)
} ## ---- Attach exposure infor to u_Summary_report, CNV_s


