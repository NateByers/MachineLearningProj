# Most R code is from Applied Predictive Modeling (2013) by Kuhn and Johnson

library(kernlab)
library(caret)
library(pROC)
# library the genefilter package (install from bioconductor first if it's not 
# installed already)
if(!require(genefilter)){
  source("http://bioconductor.org/biocLite.R")
  biocLite("genefilter")
  library(genefilter)
}

load("../Data/split_data.RData")
# make stage a factor
train_data$stage <- as.factor(train_data$stage)
test_data$stage <- as.factor(test_data$stage)

# filter using the genefilter package http://www.bioconductor.org/packages/release/bioc/vignettes/genefilter/inst/doc/howtogenefilter.pdf

# create a function that filters out miRNAs 
# (1) At least 20% of samples should have raw intensity greater than 100; 
# (2) The coefficient of variation (sd/mean) is between 0.7 and 10.
# https://www.biostars.org/p/86981/

# create the filter function
exprFilterFun <- filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))
# make a matrix from the train_data data frame
train_matrix <- train_data[, -c(1:2)]
# transpose the matrix so that rows are miRNAs and columns are samples
train_matrix_t <- t(train_matrix)
# create a logical vector that indicates TRUE for features that meet the 
# filter requirments
expr_filter <- genefilter(train_matrix_t, exprFilterFun)
# see how many features are left
sum(expr_filter)
# subset the training data down to the features that meet the requirements
train_data_filtered <- train_data[, c(TRUE, TRUE, expr_filter)]
# subset test data
test_data_filtered <- test_data[, names(train_data_filtered)]

## make a tuning grid  that covers a range of sigma and cost values
## (values from Kuhn and Johnson 2013)  
svmrGrid <- expand.grid(sigma = c(.00007, .00009, .0001, .0002),
                        C = 2^(-3:8))

## Evaluate the model using overall 10-fold cross-validation (default for method "cv")
ctrl <- trainControl(method = "cv",
                      summaryFunction = twoClassSummary,
                      classProbs = TRUE)

# fit model
set.seed(477)
fit <- train(stage ~ .,
             data = train_data_filtered[, -1],
             method = "svmRadial",
             tuneGrid = svmrGrid,
             #preProc = c("center", "scale"),
             metric = "ROC",
             trControl = ctrl)
fit

fit_pred <- predict(fit, newdata = test_data_filtered[, - c(1:2)])
fit_prob <- predict(fit, newdata = test_data_filtered[, - c(1:2)],
                                type = "prob")

fit.df <- data.frame(stage = test_data_filtered$stage, fit_pred, 
                     early_prob = fit_prob$early,
                     late_prob = fit_prob$late)

head(fit.df)

sensitivity(data = fit.df$fit_pred,
            reference = test_data_filtered$stage,
            positive = "late")

specificity(data = fit.df$fit_pred,
            reference = test_data_filtered$stage,
            negative = "early")

posPredValue(data = fit.df$fit_pred,
            reference = test_data_filtered$stage,
            positive = "late")

negPredValue(data = fit.df$fit_pred,
            reference = test_data_filtered$stage,
            negative = "early")

confusionMatrix(data = fit.df$fit_pred,
                reference = test_data_filtered$stage,
                positive = "late")

rocCurve <- roc(response = test_data_filtered$stage,
                predictor = fit.df$early_prob,
                ## This function assumes that the second
                ## class is the event of interest, so we
                ## reverse the labels
                levels = rev(levels(test_data_filtered$stage)))

auc(rocCurve)

ci.roc(rocCurve)

plot(rocCurve, legacy.axes = TRUE)

###############################################################
# Recursive feature elimination for SVM is computationally prohibitive

# # use recursive feature elimination (pg. 494-95, Kuhn and Johnson)
# 
# # first, remove columns that had a pvalue above .1
# sig_rna <- miRNA_pvalue[miRNA_pvalue$pvalue < .1, "miRNA_id"]
# train_data_sig <- data.frame(train_data[, 1:2], 
#                              train_data[, sig_rna])
# 
# ## This summary function is used to evaluate the models.
# fiveStats <- function(...) c(twoClassSummary(...), defaultSummary(...))
# 
# ## The candidate set of the number of predictors to evaluate
# varSeq <- seq(2, ncol(train_data_sig[, -c(1:2)]) - 1, by = 2)
# 
# ## The rfe() function in the caret package is used for recursive feature 
# ## elimiation. 
# 
# ctrl$functions <- caretFuncs
# ctrl$functions$summary <- fiveStats
# 
# 
# cvCtrl <- trainControl(method = "cv",
#                        verboseIter = FALSE,
#                        classProbs = TRUE,
#                        allowParallel = FALSE)
# then <- Sys.time()
# set.seed(721)
# svmRFE <- rfe(train_data_sig[, -c(1:2)],
#               train_data_sig$stage,
#               sizes = varSeq,
#               rfeControl = ctrl,
#               metric = "ROC",
#               ## Now arguments to train() are used.
#               method = "svmRadial",
#               tuneLength = 12,
#               #preProc = c("center", "scale"),
#               trControl = cvCtrl)
# Sys.time() - then
# svmRFE

######################################################################################

# Use p-values to trim and use PCA for preprocessing

# first, remove columns that had a pvalue above .1
sig_rna <- miRNA_pvalue[miRNA_pvalue$pvalue <= .1, "miRNA_id"]
train_data_sig <- data.frame(train_data[, 1:2], 
                             train_data[, sig_rna])
test_data_sig <- data.frame(test_data[, 1:2],
                            test_data[, sig_rna])

## make a tuning grid  that covers a range of sigma and cost values
## (values from Kuhn and Johnson 2013)  
svmrGrid <- expand.grid(sigma = c(.00007, .00009, .0001, .0002),
                        C = 2^(-3:8))

## Evaluate the model using overall 10-fold cross-validation (default for method "cv")
ctrl <- trainControl(method = "cv",
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE)
# fit model
then <- Sys.time()
set.seed(477)
fit_pca <- train(stage ~ .,
             data = train_data_sig[, -1],
             method = "svmRadial",
             tuneGrid = svmrGrid,
             preProc = c("center", "scale", "pca"),
             metric = "ROC",
             trControl = ctrl)
Sys.time() - then
fit_pca

fit_pca_pred <- predict(fit_pca, newdata = test_data_sig[, - c(1:2)])
fit_pca_prob <- predict(fit_pca, newdata = test_data_sig[, - c(1:2)],
                    type = "prob")

fit_pca.df <- data.frame(stage = test_data_sig$stage, fit_pca_pred, 
                     early_prob = fit_pca_prob$early,
                     late_prob = fit_pca_prob$late)

head(fit_pca.df)

confusionMatrix(data = fit_pca.df$fit_pca_pred,
                reference = test_data_sig$stage,
                positive = "late")

pca_rocCurve <- roc(response = test_data_sig$stage,
                predictor = fit_pca.df$early_prob,
                ## This function assumes that the second
                ## class is the event of interest, so we
                ## reverse the labels
                levels = rev(levels(test_data_sig$stage)))

auc(pca_rocCurve)

ci.roc(pca_rocCurve)

plot(pca_rocCurve, legacy.axes = TRUE)



