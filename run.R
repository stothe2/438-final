#----------------
# Setup libraries
#----------------
library(ade4) # For dudi.pca
library(adegenet) # For colors
library(dplyr) # For data frame manipulation
library(gplots) # For heatmap2

#-----------
# Setup data
#-----------
# GTEx, GWAS, and OMIM data accessed on April 27, 2016.

# Gene expression matrix with rows representing genes
Adipose.expr <- read.delim("data/GTEx-expr/Adipose_Subcutaneous_Analysis.expr.txt", row.names=1)
# Covariates matrix containing top 3 genotype PCs, gender, and PEER factors
Adipose.covariates <- read.delim("data/GTEx_Analysis_V6_eQTLInputFiles_covariates/Adipose_Subcutaneous_Analysis.covariates.txt", row.names=1)
# Age and Gender of all individuals who contributed data for all tissues
subjects <- read.delim("data/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt", row.names=1)
# Age and Gender of individuals who contributed to only adipose subcutaneous tissue data
Adipose.subjects <- subjects[gsub("[[:punct:]]", x=colnames(Adipose.covariates), replacement="-"),]
# Disease gene information
gwas <- read.delim("data/gwas_catalog_v1.0.1-associations_e84_r2016-04-24.tsv")

#---------------------------------------------------------------
# Obtain data statistics - %age of female/male, age distribution
#---------------------------------------------------------------
Adipose.dist <- table(Adipose.subjects[,1:2])
print("age and gender stats for adipose subcutaneous tissue samples")
print("1: male, 2: female")
print(Adipose.dist)

#------------------------------------------------
# Correcting for confounding factors based on PCA
#------------------------------------------------
Adipose <- as.data.frame(t(Adipose.expr)) # Rows now represent individuals
t.Adipose.covariates <- as.data.frame(t(Adipose.covariates))

pca <- dudi.pca(Adipose, center=TRUE, scale=TRUE, scannf=FALSE, nf=5)

eig.perc <- pca$eig/sum(pca$eig) * 100 # Percentage of variance
print("Plotting scree plot for original gene expression matrix...")
barplot(eig.perc[1:10], ylim=c(0,100), ylab="Percentage of Variance",
  xlab="Principal Component", main="Scree Plot (original GeneExpr matrix)")

loadings <- pca$l1 # PCs

# We correct for gender, top 3 genotype PCs, and top 5 gene expression PCs
M0 <- lm(as.matrix(Adipose) ~ as.matrix(t.Adipose.covariates$gender)
         + as.matrix(t.Adipose.covariates[,1:3]) + as.matrix(loadings))

# Following Pickrell et al. (2010) approach, we're going to use the residuals
# as our new gene expression matrix for age - gene expression analysis
E <- M0$residuals

pca1 <- dudi.pca(E, center=TRUE, scale=TRUE, scannf=FALSE, nf=5)
eig.perc1 <- pca1$eig/sum(pca1$eig) * 100 # Percentage of variance
print("Plotting scree plot for new gene expression matrix...")
barplot(eig.perc1[1:10], ylim=c(0,100), ylab="Percentage of Variance",
  xlab="Principal Component", main="Scree Plot (new GeneExpr matrix)")

#-------------------------------------------
# Linear regression for aging gene detection
#-------------------------------------------
all.M1 <- list()
age.genes <- vector()
age.pvals <- vector()
age.coefs <- vector()
print("Beginning linear regression for each gene...")
for (j in 1:ncol(E)) { # For each gene
  M1 <- lm(as.matrix(E[,j]) ~ as.matrix(as.numeric(Adipose.subjects$AGE))
           + as.matrix(t.Adipose.covariates$gender) + as.matrix(t.Adipose.covariates[,1:3])
           + as.matrix(loadings))
  all.M1[[j]] <- M1
  age.genes <- append(age.genes, colnames(E)[j])
  age.pvals <- append(age.pvals, anova(M1)$'Pr(>F)'[1])
  age.coefs <- append(age.coefs, as.numeric(M1$coefficients[2]))
}

cat("# of aging genes found using only p-values (<0.05) = ",
  length(age.pvals[age.pvals < 0.05]), "\n", sep=" ")

# False discovery rate adjustment using Benjamini Hochberg method
qvals <- p.adjust(age.pvals, method="BH")

cat("# of aging genes found after FDR adjustment (<0.05) = ",
  length(qvals[qvals < 0.05]), "\n", sep=" ")

df <- cbind.data.frame(age.genes, age.coefs, age.pvals, qvals)
df.q <- df[df$qvals < 0.05,]

t.E <- t(E) # Rows now represent genes
t.E.q <- t.E[df.q$age.genes,] # Expression levels for each aging gene
df.q <- mutate(df.q, mean=.rowMeans(t.E.q[,1:75], dim(t.E.q)[1], 75)) # Expression level means of top 25% samples (~75 samples)
df.q <- arrange(df.q, desc(mean)) # Sort in descending order of means
df.final <- df.q[1:(dim(df.q)[1]-48),] # Remove bottom 20% aging genes (~48 genes)

#--------------------------------------------------------------------
# Q-Q and other plots for a strongly expressed (low qvalue) aging gene
#--------------------------------------------------------------------
print("Plotting Q-Q and other plots for ENSG00000048540...")
index <- which(colnames(E)=="ENSG00000048540.10")
plot(all.M1[[index]])

#---------------------------------------------------------
# Plotting heatmap for top age-associated gene expressions
#---------------------------------------------------------
t.E.final <- t.E[df.final$age.genes,]
#distm <- dist(t.E.final, method="euclidean")
#hclustm <- hclust(distm, method="ward.D2")
#dendm <- as.dendrogram(hclustm)
#heatmap(t.E.final, Rowv=dendm, Colv=NA, scale="column")

#----------------------------------------------------
# Disease-aging gene link detection and KEGG pathways
#----------------------------------------------------
disease.genes <- read.delim("davidtools/disease_genes.txt" ) # From David tools
kegg.genes <- read.delim("davidtools/kegg_genes.txt") # From David tools

d <- select(disease.genes, Genes, Count, Term, PValue, FDR)
k <- select(kegg.genes, Genes, Count, Term, PValue, FDR)
