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

load("./Data/split_data.RData")
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

fit_pred <- predict(fit, newdata = test_data[, - c(1:2)])
fit_prob <- predict(fit, newdata = test_data[, - c(1:2)],
                                type = "prob")

fit.df <- data.frame(stage = test_data$stage, fit_pred, 
                     early_prob = fit_prob$early,
                     late_prob = fit_prob$late)

head(fit.df)

sensitivity(data = fit.df$fit_pred,
            reference = test_data$stage,
            positive = "late")

specificity(data = fit.df$fit_pred,
            reference = test_data$stage,
            negative = "early")

posPredValue(data = fit.df$fit_pred,
            reference = test_data$stage,
            positive = "late")

negPredValue(data = fit.df$fit_pred,
            reference = test_data$stage,
            negative = "early")

confusionMatrix(data = fit.df$fit_pred,
                reference = test_data$stage,
                positive = "late")

rocCurve <- roc(response = test_data$stage,
                predictor = fit.df$early,
                ## This function assumes that the second
                ## class is the event of interest, so we
                ## reverse the labels
                levels = rev(levels(test_data$stage)))

auc(rocCurve)

ci.roc(rocCurve)

plot(rocCurve, legacy.axes = TRUE)
