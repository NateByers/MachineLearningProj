load("kirc_mirna.RData")
load("kirc_clin.RData")


# code tumor stage
clin_os$stage <- clin_os$ajcc_pathologic_tumor_stage
clin_os[clin_os$stage == "Stage I" | clin_os$stage == "Stage II", "stage"] <- "early"
clin_os[clin_os$stage == "Stage III" | clin_os$stage == "Stage IV", "stage"] <- "late"

# get p-values
ids <- names(mir_tumor)[-1]
ids <- strsplit(ids, "-", fixed = FALSE)
ids <- sapply(ids, function(x) paste(x[1:3], collapse = "-"))
ids <- data.frame(ids, stringsAsFactors = FALSE)
ids_stage <- merge(ids, clin_os[, c("bcr_patient_barcode", "stage")],
                   all.x = TRUE, by.x = "ids", by.y = "bcr_patient_barcode")
ids_stage$stage <- factor(ids_stage$stage)

pvalue <- apply(mir_tumor[, -1], 1, function(row, stage){
  t.test(as.numeric(row), as.numeric(stage))$p.value
}, stage = ids_stage$stage)

# rank miRNAs
miRNA_pvalue <- data.frame(miRNA_id = mir_tumor$miRNA_id, pvalue, stringsAsFactors = FALSE)
miRNA_pvalue <- miRNA_pvalue[order(miRNA_pvalue$pvalue),]
