# TCGA-BRCA ER Status Classification
# Author: Semiha Kartal
# Description:
# Supervised machine learning–based classification of estrogen receptor (ER)
# status in breast cancer using TCGA-BRCA RNA-seq data.


getwd()
setwd("/Users/semihakartal/Desktop/Mini Project")  
list.files()


#Load the matrix and metadata
expr <- read.delim("TCGA-BRCA.star_fpkm-uq.tsv.gz", row.names = 1, check.names = FALSE)
dim(expr)
expr[1:3, 1:3]

#Load metadata
untar("gdc_download_20251212_162227.550252.tar", exdir = "gdc_clinical")
list.files("gdc_clinical", recursive = TRUE)
metadata <- read.delim("gdc_clinical/nationwidechildrens.org_clinical_patient_brca.txt", stringsAsFactors = FALSE)

View(metadata)
head(metadata)
dim(metadata)

metadata <- metadata[-c(1,2), ]   #because first two rows repeats the column names




#Match the sample ID's in metadata and expression matrix
expr <- expr[, grepl("-01A$", colnames(expr))] #First we only keep the -01A samples which are primary tumors, there are also normal tissues o we only keep the tumors


expr_sample_ids <- colnames(expr)
expr_patient_ids <- substr(expr_sample_ids, 1, 12)  #extract ID's from expression matrix (expression ID's are sample level and metadata ID's are patient level)

metadata_patient_ids <- metadata$bcr_patient_barcode #extract ID's from metadata 

length(intersect(expr_patient_ids, metadata_patient_ids)) #check how many they overlap
table(expr_patient_ids %in% metadata_patient_ids) #basically same as above (to check) but different version 


colnames(expr) <- substr(colnames(expr), 1, 12) #we remove -01A info so that expr and metadata sample names match
head(colnames(expr)) #check

rownames(metadata) <- metadata$bcr_patient_barcode
common_ids <- intersect(colnames(expr), rownames(metadata))
length(common_ids)  

expr <- expr[, common_ids]
metadata <- metadata[common_ids, ]

all(colnames(expr) == rownames(metadata)) #Should Return TRUE -- now they are in the same order 





#Keep only ER +/- samples 

table(metadata$er_status_by_ihc) #Because we will build the classifier based on ER status, we check if it exist in the file

keep <- metadata$er_status_by_ihc %in% c("Positive", "Negative")

metadata_er <- metadata[keep, ]
expr_er <- expr[, keep]

table(metadata_er$er_status_by_ihc)  #Check
all(colnames(expr_er) == rownames(metadata_er)) #Check



#Arrange Gene names
gene_clean <- sub("\\..*$", "", rownames(expr_er))
sum(duplicated(gene_clean)) #if >0 then remove duplicates - keep one row per gene

# compute mean per row
row_means <- rowMeans(expr_er, na.rm = TRUE)

# for each gene_clean, keep the row with highest mean
keep_idx <- tapply(seq_along(gene_clean), gene_clean, function(ii) ii[which.max(row_means[ii])])
keep_idx <- unlist(keep_idx, use.names = FALSE)
expr_er <- expr_er[keep_idx, ]
rownames(expr_er) <- gene_clean[keep_idx]
any(grepl("\\.", rownames(expr_er)))     # should be FALSE
any(duplicated(rownames(expr_er)))       # should be FALSE
head(rownames(expr_er))




#Filter out ESR1 expression so that the model wouldn't cheat by directly using ESR1 itself to decide the label (label leakage)
expr_er <- expr_er[rownames(expr_er) != "ENSG00000091831", ]
any(grepl("ENSG00000091831", rownames(expr_er))) #verify ESR1 removal - should return FALSE

sum(duplicated(gene_clean))
dim(expr_er)



#There are duplicated genes so below we keep the row with the highest mean expression per gene
row_means <- rowMeans(expr_er) #here we compute average expression per row

keep_idx <- tapply(            #For each cleaned gene ID, select the row with highest mean (most expressed version is kept)
  seq_along(gene_clean),
  gene_clean,
  function(ii) ii[which.max(row_means[ii])]
)

keep_idx <- unlist(keep_idx, use.names = FALSE)   #convert to a vector because tapply() returns a list-like object.

expr_er <- expr_er[keep_idx, ]
rownames(expr_er) <- gene_clean[keep_idx]
row_means <- rowMeans(expr_er, na.rm = TRUE)
#sanity check
any(grepl("\\.", rownames(expr_er)))      # FALSE
any(duplicated(rownames(expr_er)))        # FALSE
dim(expr_er)
all(colnames(expr_er) == rownames(metadata_er))  #TRUE



#Make the lebel vector "y"
y <- factor(metadata_er$er_status_by_ihc, levels = c("Negative", "Positive"))  #we create the target labels for classification, "negative" is the baseline

#Make the future matrix "x" (samples x genes)
X <- t(expr_er)

#Select TOP 5000 VARIABLE GENES before scaling
gene_var <- apply(X, 2, var, na.rm = TRUE)
summary(gene_var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:3000]
X_top <- X[, top_genes]
dim(X_top)


# Now scale -> Z-score (standardize) each gene
X_top_scaled <- scale(X_top) #we make here mean=0 sd=1 Because many ML models behave better when features are on comparable scales.

# Combine into one data frame for modeling
data_for_ml <- data.frame(X_top_scaled, ER_status = y) #This creates a single object that has: all gene features + one column at the end with ER_status
#check
dim(X_top)        # should be 1030 x 5000
dim(data_for_ml)  # 1030 x 5001


#SPLIT INTO TRAIN/TEST
set.seed(1)
idx <- sample(seq_len(nrow(X_top)), size = floor(0.8*nrow(X_top)))

X_train <- X_top[idx, ]
X_test  <- X_top[-idx, ]

y_train <- y[idx]
y_test  <- y[-idx]

#Scale using training statistics only (avoid leakage)
mu <- colMeans(X_train)
sdv <- apply(X_train, 2, sd)
sdv[sdv == 0] <- 1

X_train_sc <- sweep(sweep(X_train, 2, mu, "-"), 2, sdv, "/")
X_test_sc  <- sweep(sweep(X_test,  2, mu, "-"), 2, sdv, "/")


#install.packages(c("glmnet", "pROC"))   # run once if not installed
library(glmnet)
library(pROC)

#PREPARE DATA FOR glmnet #glmnet wants a matrix for X and a 0/1 outcome (or factor, but 0/1 is simplest)
# Make sure X's are matrices
Xtr <- as.matrix(X_train_sc)
Xte <- as.matrix(X_test_sc)

# Convert labels to 0/1: Positive = 1, Negative = 0
ytr <- ifelse(y_train == "Positive", 1, 0)
yte <- ifelse(y_test  == "Positive", 1, 0)

#CROSS-VALIDATED ELASTIC NET REGRESSION
set.seed(1)
cvfit <- cv.glmnet(
  x = Xtr,
  y = ytr,
  family = "binomial",
  alpha = 0.5,
  nfolds = 5,
  type.measure = "auc"
)

plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se

#PREDICT ON THE TEST SET + COMPUTE AUC
prob_test <- predict(cvfit, newx = Xte, s = "lambda.1se", type = "response")[,1]

roc_obj <- roc(response = yte, predictor = prob_test)
auc(roc_obj)
plot(roc_obj)

#CONFUSION MATRIX (optional, uses a 0.5 threshold)
pred_class <- ifelse(prob_test >= 0.5, 1, 0)

table(
  Predicted = pred_class,
  True = yte
)

#How many genes did the model use?
coef_mat <- coef(cvfit, s = "lambda.1se")
nonzero <- sum(coef_mat != 0) - 1  # minus intercept
nonzero


coef_mat <- coef(cvfit, s = "lambda.1se")
selected_genes <- rownames(coef_mat)[coef_mat[,1] != 0]
selected_genes




