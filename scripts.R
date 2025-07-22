rm(list = ls())

# set work path

work.path <- getwd()
setwd(work.path) 

code.path <- file.path(work.path, "Codes")
data.path <- file.path(work.path, "InputData")
res.path <- file.path(work.path, "Results")
fig.path <- file.path(work.path, "Figures")

if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)

# install.packages("randomForestSRC")
# install.packages("snowfall")
# install.packages("xgboost")

library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)

# Loading scripts for model training and evaluation
source(file.path(code.path, "ML.R"))
# Select the final generated model type: panML represents generating models constructed by different algorithms; MultiLogistic represents extracting variables used in other models and establishing a multivariate logistic model
FinalModel <- c("panML", "multiLogistic")[2]

## Training Cohort ---------------------------------------------------------
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_class <- read.table(file.path(data.path, "Training_class.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
# Extract common samples from the training set
comsam <- intersect(rownames(Train_class), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_class <- Train_class[comsam,,drop = F]

## Validation Cohort -------------------------------------------------------
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Test_class <- read.table(file.path(data.path, "Testing_class.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
# Extract common samples from the test set
comsam <- intersect(rownames(Test_class), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_class <- Test_class[comsam,,drop = F]

# Extract common gene
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) 
Test_expr <- t(Test_expr[comgene,]) 

#Standardize data separately by queue (adjust centerFlags and scaleFlags according to the situation)
##Data: Requires expression profile data (behavioral samples, listed as genes)
##Cohort: The queue to which the sample belongs is a vector. When no value is input, the default full expression matrix comes from the same queue
##CenterFlag/scaleFlags: whether to standardize the gene mean/standard deviation to 1;
##The default parameter is NULL, indicating that no standardization will be performed;
##When T/F, it means to standardize/not standardize all queues
##When inputting a vector composed of T/F, process the queue in order, and the length of the vector should be the same as the number of queues
##If centerFlags=c (F, F, F, T, T), it means standardizing the 4th and 5th queues, and the order of the flags should be consistent with the queue order
##If centerFlags=c ("A"=F, "C"=T, "B"=F), it means that the queue C is standardized, and the order of the flags is not required to be consistent with the data
Train_set = scaleData(data = Train_expr, centerFlags = T, scaleFlags = T) 
names(x = split(as.data.frame(Test_expr), f = Test_class$Cohort))
Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = T, scaleFlags = T)
# summary(apply(Train_set, 2, var))
# summary(apply(Test_set, 2, var))
# lapply(split(as.data.frame(Test_set), Test_class$Cohort), function(x) summary(apply(x, 2, var))) # 测试scale结果

# Model training and validation -------------------------------------------

## method list --------------------------------------------------------
#Record the model that needs to be run here in the format of Algorithm 1 Name [Algorithm Parameters]+Algorithm 2 Name [Algorithm Parameters]
#Currently, only Stepglm and Enet support inputting algorithm parameters
methods <- read.xlsx(file.path(code.path, "methods.xlsx"), startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)

## Train the model --------------------------------------------------------

classVar = "outcome"
min.selected.var = 2 

##Pre training summarizes the variable screening process used by each method to reduce computational complexity
Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+") # Review all methods and analyze whether each method requires variable pre screening
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1]) # Delete the algorithms used to construct classification models and retain the algorithms used for variable filtering
preTrain.method = unique(unlist(preTrain.method)) # Summarize all variable filtering algorithms and remove duplicate calculations

preTrain.var <- list() #Used to store variables selected by various algorithms
set.seed(seed = 12345) # Set modeling seeds to make the results reproducible
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, 
                                 Train_set = Train_set,
                                 Train_label = Train_class, 
                                 mode = "Variable",     
                                 classVar = classVar) 
}
preTrain.var[["simple"]] <- colnames(Train_set)

## Model training
model <- list() 
set.seed(seed = 777) 
Train_set_bk = Train_set 
for (i in 1:length(methods)){
  method=methods[i]
  cat(match(method, methods), ":", method, "\n")
  method_name = method
  method <- strsplit(method, "\\+")[[1]] 
  
  if (length(method) == 1) method <- c("simple", method) #If this method does not pre screen variables, it is considered to have used the simple method for variable screening
  
  Variable = preTrain.var[[method[1]]] # Call the result filtered by the previous variable based on the first value of the method name
  Train_set = Train_set_bk[, Variable]   # Taking a subset of the training set, because there is an algorithm written by the original author that has some issues and cannot pass parameters properly
  Train_label = Train_class            # So it is necessary to modify the variable name here to avoid the function calling the object incorrectly
  model[[method_name]] <- RunML(method = method[2],       
                                Train_set = Train_set,     # Variables with potential predictive value in the training set
                                Train_label = Train_label,
                                mode = "Model",           
                                classVar = classVar)      

  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
Train_set = Train_set_bk; rm(Train_set_bk)
saveRDS(model, file.path(res.path, "model.rds")) 

if (FinalModel == "multiLogistic"){
  logisticmodel <- lapply(model, function(fit){ 
    tmp <- glm(formula = Train_class[[classVar]] ~ .,
               family = "binomial", 
               data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit) # Extract the predictor variables that will ultimately be used in the logistic model
    return(tmp)
  })
}
saveRDS(logisticmodel, file.path(res.path, "logisticmodel.rds")) 

## Evaluate the model -----------------------------------------------------

model <- readRDS(file.path(res.path, "model.rds"))
# model <- readRDS(file.path(res.path, "logisticmodel.rds"))
methodsValid <- names(model)

#Calculate the sample risk score based on the given expression level
#Prediction probability
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set)) # 2.0
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(res.path, "RS_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) 

Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set))
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
#Class_mat <- cbind.data.frame(Test_class, Class_mat[rownames(Class_mat),])
write.table(Class_mat, file.path(res.path, "Class_mat.txt"), 
            sep = "\t", row.names = T, col.names = NA, quote = F)

fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}

fea_df <- lapply(model, function(fit){
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file.path(res.path, "fea_df.txt"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

# Calculate C-index for each model
AUC_list <- list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],     
                                Test_set = Test_set,   
                                Test_label = Test_class,   
                                Train_set = Train_set,   
                                Train_label = Train_class, 
                                Train_name = "Merge",      
                                cohortVar = "Cohort",     #Important: Used to specify the variable for the queue, the column must exist and be specified [default is' Cohort '], otherwise an error will be reported
                                classVar = classVar)       
}
AUC_mat <- do.call(rbind, AUC_list)
write.table(AUC_mat, file.path(res.path, "AUC_mat.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------

AUC_mat <- read.table(file.path(res.path, "AUC_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_AUC <- apply(AUC_mat, 1, mean)       
avg_AUC <- sort(avg_AUC, decreasing = T)   
AUC_mat <- AUC_mat[names(avg_AUC), ]  
fea_sel <- fea_list[[rownames(AUC_mat)[1]]] #The features selected by the optimal model (with the highest AUC mean in the test set)
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))

if(ncol(AUC_mat) < 3) { 
  CohortCol <- c("red","blue") 
} else { 
  CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired") 
}
names(CohortCol) <- colnames(AUC_mat)

cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(AUC_mat, 
                    avg_AUC, 
                    CohortCol, "steelblue", 
                    cellwidth = cellwidth, cellheight = cellheight, 
                    cluster_columns = F, cluster_rows = F)

pdf(file.path(fig.path, "AUC.pdf"), width = cellwidth * ncol(AUC_mat) + 4.5, height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm)
invisible(dev.off())


fea_df1=dplyr::filter(fea_df,algorithm == "GBM")

fit <- gbm(formula = Train_label[[classVar]] ~ .,
           data = as.data.frame(Train_set),
           distribution = 'bernoulli',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
best <- which.min(fit$cv.error)
fit <- gbm(formula = Train_label[[classVar]] ~ .,
           data = as.data.frame(Train_set),
           distribution = 'bernoulli',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001, n.cores = 8)
fit$subFeature = colnames(Train_set)

# Extracting the relative influence (importance) of the model
importance <- summary(fit, plotit = FALSE)

importance_df <- data.frame(
  Feature = importance$var,
  Importance = importance$rel.inf
)

importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]

print(head(importance_df, 20))

pdf("gene_importance_plot.pdf", width = 10, height = 8)

top_n <- 10
top_features <- head(importance_df, top_n)

top_features <- top_features[order(top_features$Importance), ]
top_features$Feature <- factor(top_features$Feature, levels = top_features$Feature)

barplot(top_features$Importance, 
        names.arg = top_features$Feature, 
        horiz = TRUE, 
        las = 1, 
        cex.names = 0.7,
        col = "steelblue",
        main = "Top Gene Features by Importance",
        xlab = "Relative Importance (%)")

dev.off()
library(ggplot2)
top_n <- 10
top_features <- head(importance_df, top_n)

top_features$Feature <- factor(top_features$Feature, 
                               levels = top_features$Feature[order(top_features$Importance)])

p <- ggplot(top_features, aes(x = Importance, y = Feature)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Top Gene Features by Importance",
       x = "Relative Importance (%)",
       y = NULL) +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("gene_importance_ggplot.pdf", p, width = 10, height = 8)
hubgene=top_features$Feature

write.table(hubgene,file = "machine_learing_gene.txt",quote = F,row.names = F)
