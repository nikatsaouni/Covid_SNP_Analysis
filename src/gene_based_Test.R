library(COMBAT)
library(dplyr)

# Function of VEGAS with different proportion tests and Fisher combination test
vegas_f = function(x, cor_G, vegas.pct=c(0.1,0.2,0.3,0.4,1), max.simulation=1e4){
  pval_vegas <- vegas.call(x=x, cor_G=cor_G, vegas.pct=vegas.pct, n_simul=max.simulation)
  pval_vegas
}

#this is not to be called directly
vegas.call = function(x, cor_G, vegas.pct, n_simul){
  stopifnot(length(x) == ncol(cor_G))
  vegas_vec <- ceiling(vegas.pct*ncol(cor_G))
  vegas_vec <- sort(vegas_vec)
  if(vegas_vec[1]>1){
    vegas.pct <- c(0,vegas.pct)
    vegas_vec <- c(1,vegas_vec)
  }
  chisq_vec <- qchisq(x,1,lower.tail=FALSE)
  chisq_vec[chisq_vec == Inf] <- 60
  n_snps <- length(x)
  n_tests <- length(vegas_vec)
  
  TS_obs <- rep(NA,n_tests)
  TS_obs[1] <- max(chisq_vec, na.rm=TRUE)
  chisq_vec <- sort(chisq_vec, decreasing = TRUE)
  for (j in 2:n_tests) TS_obs[j] <- sum(chisq_vec[1:vegas_vec[j]])
  
  rd  <- rmvnorm(n_simul, mean=rep(0,n_snps),sigma=cor_G)
  rd2 <- rd^2
  rd2 <- apply(rd2,1,sort,decreasing=TRUE)
  
  pPerm0 <- rep(NA,n_tests)
  T0s <- apply(rd2,2,max)
  pPerm0[1]<- (sum(T0s >= TS_obs[1])+1)/(length(T0s)+1)
  for(j in 2:n_tests){
    for (i in 1:n_simul) T0s[i] <- sum(rd2[1:vegas_vec[j],i])
    pPerm0[j] <- (sum(T0s >= TS_obs[j])+1)/(length(T0s)+1)
  }
  v1 <- paste0('VEGAS.p',vegas.pct)
  v1[vegas_vec==ncol(cor_G)] <- 'VEGAS.all'
  v1[vegas_vec==1] <- 'VEGAS.max'
  names(pPerm0) <- v1
  pPerm0
}

# Function of GATES
gates_f = function(x, cor_G){
  if(is.positive.definite(cor_G)==FALSE) stop('cor_G is not positive definite. Please re-calculate with ld.Rsquare function.\n')
  pval_sort  <- sort(x)
  pval_order <- order(x)
  n_snps     <- length(x)
  
  cor_P <- cor_G[pval_order, pval_order]
  cor_P <- 0.2982*cor_P^6 - 0.0127*cor_P^5 + 0.0588*cor_P^4 + 0.0099*cor_P^3 + 0.6281*cor_P^2 - 0.0009*cor_P
  
  p_gates <- ext_simes_f(pval_sort, cor_P)
  p_gates
}

# Function of extended Simes
ext_simes_f = function(x, cor_r){
  eff.snpcount.fun <- function(ldmat) {
    ldmat <- as.matrix(ldmat)
    snpcount.local <- dim(ldmat)[1]
    if (snpcount.local <= 1) return(1)
    ev <- eigen(ldmat, only.values = TRUE)$values
    if (sum(ev < 0) != 0) {
      ev <- ev[ev > 0]
      ev <- ev/sum(ev) * snpcount.local
    }
    ev <- ev[ev > 1]
    snpcount.local - sum(ev - 1)
  }
  eff.snpcount.global <- eff.snpcount.fun(cor_r)
  
  n_values <- length(x)
  candid <- sapply(1:n_values, function(i){
    (eff.snpcount.global * x[i])/eff.snpcount.fun(cor_r[1:i,1:i])
    
  })
  
  p_ext_simes <- min(candid)
  p_ext_simes
}


#### Read P-values and SNPs
data_snp_p <- read.csv2("/Users/nikatsaouni/Documents/Disease_Prediction/Covid_Proj/CovidSNPAnalysis/data/06_reporting/Case_1/log_Regression_Genes.csv", header=TRUE,sep=",")
snp.info   <- as.data.frame(data_snp_p)


#### Read reference genotype
snp.ref_   <- read.csv2("/Users/nikatsaouni/Documents/Disease_Prediction/Covid_Proj/CovidSNPAnalysis/data/06_reporting/Case_1/rs_correlat.csv", header = FALSE,sep=",")
### Keep genotype only for SNPs in a gene
snp.ref_   <- as.matrix(snp.ref_)
snp.ref_   <- as.data.frame(snp.ref_)

### make first row title (make SNP names as title)
colnames(snp.ref_) <- snp.ref_[1,]
snp.ref_ <- snp.ref_[-1, ] 

### Keep only 566 SNPs
snp.ref_ <- snp.ref_ %>%  select(snp.info$SNP)
snp.info   <- as.matrix(data_snp_p)




#### Read gene list
gene_list <- read.csv2("/Users/nikatsaouni/Documents/Disease_Prediction/Covid_Proj/Updated_Results/genes_of_interest_final.csv", header=TRUE,sep=",") ###Pablos genes of interest

gene.info   <- as.data.frame(gene_list[,1])

### For all genes
snp.pvals_ <- as.data.frame(snp.info)


#### Initialize dataframes 
gates_df <- matrix(ncol=1, nrow=nrow(gene.info))
vegas_df <- matrix(ncol=2, nrow=nrow(gene.info))
simpleM_df <- matrix(ncol=1, nrow=nrow(gene.info))
no_snps_per_gene <- matrix(ncol=1, nrow=nrow(gene.info))

for (i in 1:nrow(gene.info))
  # for (i in 1)
{
  # gene_to_test <- "KLRC3"
  gene_to_test <- gene.info[[1]][i]
  print(gene_to_test)
  snp.pvals <- snp.pvals_[grepl(gene_to_test, snp.pvals_$Gene),]
  
  ### check if gene has important SNPs
  if (dim(snp.pvals)[1] != 0)
  {
    ###ignore genes with less than 3 SNPs
    if (dim(snp.pvals)[1] > 0)
    {
      
      snplist <- list(snp.pvals$SNP)
      #snplist <- split(snplist, seq(nrow(snplist)))
      
      snp.ref1 <- snp.ref_ %>% select(snp.pvals$SNP)
      snp.ref1[] <- lapply(snp.ref1, as.numeric)
      
      
      snp.pvals <- as.matrix(snp.pvals)
      snp.pvals<- subset(snp.pvals, select = - c(Gene))
      snp.pvals <- subset(snp.pvals, select = - c(SNP))
      snp.pvals <- as.numeric(snp.pvals)
      
      if (length(snplist[[1]]) == 1)
      {
        pval_gates <- snp.pvals
        pval_vegas <- snp.pvals
        pval_simpleM <- snp.pvals
        
      }
      else {
        ### Calculate LD R SQUARED
        cor_G <- ld.Rsquare(snp.ref1)
        
        #### GATES TEST
        pval_gates <- gates_f(x=snp.pvals, cor_G=cor_G)
        pval_vegas <- vegas_f(x=snp.pvals, cor_G=cor_G, vegas.pct=c(0.1,0.2,0.3,0.4,1))
        pval_simpleM <- simpleM(x=snp.pvals, cor_G=cor_G)
        
      }
      ##gate df
      number_of_snps <- length(snplist[[1]])
      gates_df[i,] <- pval_gates
      ### number of snps
      no_snps_per_gene[i,] <- number_of_snps
      
      ##vegas df
      vegas_df[i,1] <- pval_vegas[[1]]
      a <- tail(pval_vegas, n=1)
      vegas_df[i,2] <- a[[1]]
      
      ##simpleM df
      simpleM_df[i,] <- pval_simpleM 
    }
  }  
}


#### Concatenate info
result_data <- data.frame(gene.info,gates_df,vegas_df,simpleM_df,no_snps_per_gene)
result_data <- result_data %>% mutate(across(where(is.numeric), ~ round(., digits = 6)))

result_data <- result_data[order(result_data$gates_df),]

## multiple testing corrected p-values

p.adjuste_values_G <- p.adjust(result_data$gates_df, method = p.adjust.methods, n = length(result_data$gates_df))
p.adjuste_values_vegas <- p.adjust(result_data$X1, method = p.adjust.methods, n = length(result_data$X1))
p.adjuste_values_vegas10 <- p.adjust(result_data$X2, method = p.adjust.methods, n = length(result_data$X2))
p.adjuste_values_simpleM <- p.adjust(result_data$simpleM_df, method = p.adjust.methods, n = length(result_data$simpleM_df))


result_data$GATES_adj <- p.adjuste_values_G
result_data$VEGAS_adj <- p.adjuste_values_vegas
result_data$VEGAS10_adj <- p.adjuste_values_vegas10
result_data$VEGAS_simpleM_adj <- p.adjuste_values_simpleM


write.csv(x=result_data, file="/Users/nikatsaouni/Documents/Disease_Prediction/Covid_Proj/CovidSNPAnalysis/data/06_reporting/Case_1/gene_based_filtered.csv")


