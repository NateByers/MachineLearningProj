load("Data/kirc_mirna.RData")
load("Data/kirc_clin.RData")


# code tumor stage
kirc_clin$stage <- kirc_clin$ajcc_pathologic_tumor_stage
kirc_clin[kirc_clin$stage == "Stage I" | kirc_clin$stage == "Stage II", "stage"] <- "early"
kirc_clin[kirc_clin$stage == "Stage III" | kirc_clin$stage == "Stage IV", "stage"] <- "late"


# merge clinical and tumor data
all_data <- merge(kirc_clin[, c("bcr_patient_barcode", "stage")], kirc_mirna)

# rename id column
names(all_data)[1] <- "id"



#split into training and test sets: 80/20 from early stage and 80/20 from late stage
early_data <- all_data[all_data$stage == "early", ]
set.seed(10); r_early_rows <- sample(1:nrow(early_data), nrow(early_data))
r_early_data <- early_data[r_early_rows, ]

early_train <- r_early_data[1:round(.8*nrow(r_early_data)), ]
early_test <- r_early_data[(round(.8*nrow(r_early_data)) + 1): nrow(r_early_data), ]

late_data <- all_data[all_data$stage == "late", ]
set.seed(23); r_late_rows <- sample(1:nrow(late_data), nrow(late_data))
r_late_data <- late_data[r_late_rows, ]
late_train <- r_late_data[1:round(.8*nrow(r_late_data)), ]
late_test <- r_late_data[(round(.8*nrow(r_late_data)) + 1):nrow(r_late_data), ]

train_data <- rbind(early_train, late_train)
test_data <- rbind(early_test, late_test)

# remove miRNAs that have 0 variance in the train data
no_variance <- sapply(train_data[, -c(1:2)], function(column){
  v <- var(column, na.rm = TRUE)
  if(v == 0){TRUE}else{FALSE}
  })
train_data <- train_data[, c(TRUE, TRUE, !no_variance)]


# get p-values for test set
# look at boxplot for first miRNA
boxplot(hsa.let.7a.1 ~ stage, data = train_data)

pvalue <- apply(train_data[, -c(1:2)], 2, function(column, stage){
  early <- column[stage == "early"]
  late <- column[stage == "late"]
  t.test(early, late)$p.value
}, stage = train_data$stage)

# rank train miRNAs
miRNA_pvalue <- data.frame(miRNA_id = names(train_data)[-c(1:2)], pvalue, stringsAsFactors = FALSE)
miRNA_pvalue <- miRNA_pvalue[order(miRNA_pvalue$pvalue),]


# look at variation
train_coef_var <- sapply(train_data[, -c(1:2)], function(column){
  sd(column, na.rm = TRUE)/mean(column, na.rm = TRUE)
  })
train_coef_var <- train_coef_var[order(train_coef_var)]
summary(train_coef_var)
boxplot(train_coef_var)

# remove correlated miRNAs

# function that takes a data frame with normalized columns and 
# a cutoff value for positive correlation
removeCorrFeatures <- function(data, cutoff){ 
  # data = train_data[, -c(1:2)]
  # cutoff = .75
  cor_m <- cor(data, use = "pairwise.complete.obs")
  cor_m[cor_m == 1] = 0
  col.max <- which.max(apply(cor_m, 2, max))
  row.max <- which.max(cor_m[, col.max])
  max <- cor_m[names(row.max), names(col.max)]
  
  while(max > cutoff){
    mirna1 <- cor_m[, names(col.max)]
    mirna1 <- mirna1[!(names(mirna1) %in% c(names(col.max), names(row.max)))]
    mean1 <- mean(mirna1)
    mirna2 <- cor_m[, names(row.max)]
    mirna2 <- mirna2[!(names(mirna2) %in% c(names(col.max), names(row.max)))]
    mean2 <- mean(mirna2)
    if(mean1 > mean2)
    {
      data <- data[names(data) != names(col.max)]
    }else{
      data <- data[names(data) != names(row.max)]
    }
    cor_m <- cor(data, use = "pairwise.complete.obs")
    cor_m[cor_m == 1] = 0
    col.max <- which.max(apply(cor_m, 2, max))
    row.max <- which.max(cor_m[, col.max])
    max <- cor_m[names(row.max), names(col.max)]
  }
  data
  
}

# then <- Sys.time()
# train_trimmed <- removeCorrFeatures(train_data[, -c(1:2)], cutoff = .75)
# Sys.time() - then

df <- nearZeroVar(train_data)
