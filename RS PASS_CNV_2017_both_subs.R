
## ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- functions
PASS_CNV.remove_AB <- function(x) {
    x = sub("-A$", "", x)
    x = sub("-B$", "", x)
    return(x)
}

## ---- ---- ---- ---- 
PASS_CNV.attach_pheno_info <- function(xxt, xxp, xxt_name, xxp_name, pheno_list) {
    i  = match(PASS_CNV.remove_AB(xxt[, xxt_name]),    PASS_CNV.remove_AB(xxp[, xxp_name]))
    for(iph in pheno_list) {
        xxt[, iph] = xxp[i, iph]
    }
    return(xxt)
}

## ---- ---- ---- ---- 
PASS_CNV.create_u_Summary_report <- function(h1, h2) {
    u_Summary_report           = h1$Summary_report[,c("Collaborator.Participant", "Demise.Type", "Specimen.Type", "Call.Rate")]
    u_Summary_report[,"Batch"] = "h1"
    xx                         = h2$Summary_report[,c("Collaborator.Participant", "Demise.Type", "Specimen.Type", "AutoCall")]
    names(xx)[which(names(xx) == "AutoCall")] = "Call.Rate"
    xx[,"Batch"]               = "h2"
    u_Summary_report           = rbind(u_Summary_report, xx)
    u_Summary_report[,"pass_QC"] = u_Summary_report[,"Call.Rate"] >= 98.0 
    u_Summary_report[,"Collaborator.Participant"] = PASS_CNV.remove_AB(u_Summary_report[, "Collaborator.Participant"])
    return(u_Summary_report)
}

## ---- ---- ---- ---- 
PASS_CNV.pheno_table_u_Summary_report <- function(u_Summary_report, s, fn="xx_pheno_sample_type_h1_2_2018_04_29.csv") {
#   I assume RESULTS_dir is present in the environment
    xx = table(u_Summary_report[s, c("Demise.Type", "Specimen.Type")], useNA="ifany")
    xx = cbind(xx, Total = rowSums(xx))
    xx = rbind(xx, Total = colSums(xx))
    if(fn != "") { fn = paste(RESULTS_dir, fn, sep="");   write.csv(xx, fn);  system(paste("open ", fn)) } ## 
    return(xx)
}

## ---- ---- ---- ---- 
PASS_CNV.select_CNV <- function(xx, cuts) {
    s = T
    s = s &   xx[,       "Chr"] %in%  cuts$Chr_range;  tNA(s)
    s = s & !(xx[,     "Value"] %in%  cuts$CN_exclude);  tNA(s)
    s = s &  findInterval(xx[,"Confidence"], cuts$Confidence) == 1;  tNA(s)
    s = s &  findInterval(xx[,"Size"],       cuts$Size      ) == 1;  tNA(s)
    return(s)
}

## ---- ---- ---- ---- 
PASS_CNV.define_CNV_cuts <- function() {
    CNV_cuts  = list()
    cuts      = list()

    cuts$Chr_range     = 1:22
    cuts$CN_exclude    = 2
    cuts$Confidence    = c(800,  1e6)
    cuts$Size          = c(.0, 1e6)
    CNV_cuts[[1]]      = cuts           # standard cuts

    cuts$Confidence    = c(1000,  1e6)
    CNV_cuts[[2]]      = cuts           # used to produce Klaus first list

    return(CNV_cuts)
}
CNV_cuts = PASS_CNV.define_CNV_cuts()

## ---- ---- ---- ---- 

## ---- ---- ---- ---- 

## ---- ---- ---- ---- 

## ---- ---- ---- ---- 
