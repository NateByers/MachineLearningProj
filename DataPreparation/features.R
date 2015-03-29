load("kirc_mirna.RData")
load("kirc_clin.RData")


# code tumor stage
clin_os$stage <- clin_os$ajcc_pathologic_tumor_stage
clin_os[clin_os$stage == "Stage I" | clin_os$stage == "Stage II", "stage"] <- "early"
clin_os[clin_os$stage == "Stage III" | clin_os$stage == "Stage IV", "stage"] <- "late"

#split into training and test sets: 80/20 from early stage and 80/20 from late stage
clin_early <- clin_os[clin_os$stage == "early", ]
set.seed(10); r_early_rows <- sample(1:nrow(clin_early), nrow(clin_early))
r_clin_early <- clin_early[r_early_rows, ]

clin_early_train <- r_clin_early[1:round(.8*nrow(r_clin_early)), ]
clin_early_test <- r_clin_early[(round(.8*nrow(r_clin_early)) + 1): nrow(r_clin_early), ]

clin_late <- clin_os[clin_os$stage == "late", ]
set.seed(23); r_late_rows <- sample(1:nrow(clin_late), nrow(clin_late))
r_clin_late <- clin_late[r_late_rows, ]
clin_late_train <- r_clin_late[1:round(.8*nrow(r_clin_late)), ]
clin_late_test <- r_clin_late[(round(.8*nrow(r_clin_late)) + 1):nrow(r_clin_late), ]

clin_train <- rbind(clin_early_train, clin_late_train)
clin_test <- rbind(clin_early_test, clin_late_test)


# get p-values
ids <- names(mir_tumor)[-1]
ids <- strsplit(ids, "-", fixed = FALSE)
ids <- sapply(ids, function(x) paste(x[1:3], collapse = "-"))
ids <- data.frame(ids, stringsAsFactors = FALSE)
ids_stage <- merge(ids, clin_os[, c("bcr_patient_barcode", "stage")],
                   all.x = TRUE, by.x = "ids", by.y = "bcr_patient_barcode")
ids_stage$stage <- factor(ids_stage$stage)

boxplot(as.numeric(mir_tumor[1, 2:ncol(mir_tumor)]) ~ ids_stage$stage)

pvalue <- apply(mir_tumor[, -1], 1, function(row, stage){
  early <- as.numeric(row)[stage == "early"]
  late <- as.numeric(row)[stage == "late"]
  t.test(early, late)$p.value
}, stage = ids_stage$stage)

# rank miRNAs
miRNA_pvalue <- data.frame(miRNA_id = mir_tumor$miRNA_id, pvalue, stringsAsFactors = FALSE)
miRNA_pvalue <- miRNA_pvalue[order(miRNA_pvalue$pvalue),]

