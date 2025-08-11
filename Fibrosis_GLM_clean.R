##This is the code to take our normalized head kidney reads and run our glm models looking at fibrosis##
##NOTE- these are two-way models (all treatments, time points condensed) considering the effects of fibrosis, population, and their interaction##
##at the end we assess overrepresentation of genes and then run post hoc testing##
##Specifically we focus on post-hoc tests of population level differences 
##We also consider population-specific associations between fibrosis and gene expression for genes with a sig interaction term 
##last updated Aug 11th, 2025 by LEF##

##load in relevant packages##
library(tibble)
library(data.table)

##next do initial data processing steps, first reading in normalized reads and metadata files##
  {
    # Read in data
  readcount <- read.csv("normalizedreads_hk.csv", check.names=FALSE)
  readcount=t(readcount)
  readcount=as.data.frame(readcount)
  readcount=tibble::rownames_to_column(readcount, "FishID")
  colnames(readcount) <- readcount[1,]
  readcount <- readcount[-1, ] 
  colnames(readcount)[1] = "Fish_ID"
  metadata <- read.csv("expdesign_fib_matchhund.csv")
  metadata <- na.omit(metadata)
  head(metadata)
  dim(readcount)
  dim(metadata)
  
  
  }
  ##check and make sure reads and meta are consistent with each other##
  metadata$Fish_ID=as.character(metadata$Fish_ID)
  identical(readcount[['Fish_ID']],metadata[['Fish_ID']])
  
  {
    # Read in gene info
  geneinfo <- read.csv("MapENSGACT_ENSGACG.csv")
  head(geneinfo)
  genenames <- read.csv("GeneNames.csv", header = F)
  head(genenames)
  names(genenames) <- c("ENSGACG", "genename")
  geneinfo2 <- merge(geneinfo, genenames, by = "ENSGACG")
  dim(geneinfo2)
  
  }
  
  {
    # Omit day 90
  readcount <- readcount[metadata$DPI < 90,]
  metadata <- metadata[metadata$DPI < 90,]
  }
  #Omit WA
  readcount <- readcount[metadata$Treatment %in% c("Alum", "PBS", "Worm Protein"),]
  metadata <- metadata[metadata$Treatment %in% c("Alum", "PBS", "Worm Protein"),]
  write.csv(metadata, "fib_model_meta.csv", row.names=FALSE)
  write.csv(readcount, "fib_model_reads.csv", row.names=FALSE)
  
  ##check and make sure reads and meta are consistent with each other one more time (cause I'm paranoid)
  identical(readcount[['Fish_ID']],metadata[['Fish_ID']])


######## STEP 1- model running


  {
  ####### SKIP THIS IF RERUNNING and instead reload output
  ###### GLM of alum effects to identify instances with a main effect of treatment, time, population*treatment, etc
  ##start by finding the actual quartile to use##
  #prop.matrix <- matrix(data = NA, nrow = ncol(readcount), ncol = 2)
  #prop.matrix <- as.data.frame(prop.matrix)
  #readcount <- readcount[,-1]
  #geneIDs <- names(readcount)
  #readcount=sapply(readcount,as.numeric)
  
  
  #for(gene_i in 1:ncol(readcount)){
    #readcount <- as.data.frame(readcount)
    #geneID <- geneIDs[gene_i]
    #y <- readcount[,gene_i]
    #Totalreads <- rowSums(readcount)
    #Y <- cbind(focalgene = y, otherreads = Totalreads - y)
    #Y <- as.matrix(Y)
    #Yprop <- Y[,1] / rowSums(Y)
    #mean_expr <- mean(  Yprop, na.rm = T)
    #prop.matrix[gene_i,]= c(geneID, mean_expr)}
  #prop.matrix$V2=lapply(prop.matrix$V2, as.numeric)
  #quantile(unlist(prop.matrix$V2), probs = 0.5, na.rm = TRUE)
  
  ##then actually run the model, using the detected quartile cutoff
    matrix <- matrix(data = NA, nrow = ncol(readcount), ncol = 12)
    matrix <- as.data.frame(matrix)
    readcount <- readcount[,-1]
    geneIDs <- names(readcount)
    readcount=sapply(readcount,as.numeric)
    metadata$Fib_L_count=sapply(metadata$Fib_L_count,as.numeric)
    metadata$Treatment<- relevel(factor(metadata$Treatment), ref = "PBS")
    metadata$Population<- relevel(factor(metadata$Population), ref = "SAY")
  
  for(gene_i in 1:ncol(readcount)){
    readcount <- as.data.frame(readcount)
    geneID <- geneIDs[gene_i]
    y <- readcount[,gene_i]
    Totalreads <- rowSums(readcount)
    Y <- cbind(focalgene = y, otherreads = Totalreads - y)
    Y <- as.matrix(Y)
    Yprop <- Y[,1] / rowSums(Y)
    mean_expr <- mean(  Yprop, na.rm = T)
    if(mean_expr > 2.493706e-05   ){ #50% quantile mean read count threshold
       model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      results <- anova(model, test = "Chisq")
      matrix[gene_i,] <- c( geneID, mean_expr, results$`Pr(>Chi)`,model$coefficients)
    }  else{
      matrix[gene_i,] <- c( geneID, mean_expr, NA, NA, NA, NA,NA, NA, NA,  NA, NA, NA)
    }
  }
  ##adjust matrix names
  names(matrix) <- c("gene", "meanreadprop", "NA", "Pop_P", "Fib_P", "Pop:Fib_P", "Int_C", "Pop_GvS_C","Pop_RvS_C","Fib_C", "GvS:Fib_C", "RvS:Fib_C")
  
  ##format and write out##
  head(matrix)
  matrix <- matrix[order(matrix$Fib_P, decreasing = F),]
  matrix[1:20,]
  matrix <- matrix[,-3]
  matrix <- matrix[,-6]
  dim(na.omit(matrix))
  matrix.final=na.omit(matrix)
  #dim(matrix)
  write.csv(matrix.final, "Fibrosis_GLM_matrix_results_fixed.csv", row.names=FALSE)
  write.csv(matrix, "Fibrosis_GLM_matrix_results_everything_fixed.csv", row.names=FALSE)
  }




######  RESTART HERE for sumarizing/testing overrepresentation
  matrix <- read.csv( "Fibrosis_GLM_matrix_results_fixed.csv")
  head(matrix)
  
  ### Are there more P < 0.01 than we expect?
  alpha <- 0.01
  # Number of significant treatment main effects
    Ntests <- nrow(matrix)
  pop.sig <- matrix$Pop_P < alpha
  fib.sig <- matrix$Fib_P < alpha
  popfib.sig <- matrix$Pop.Fib_P < alpha
  sum(pop.sig) # 10599
  sum(fib.sig) #2888
  sum(popfib.sig)  #389
  
  # number of main effects of treatment
  sum(fib.sig) # 2888 fibrosis effects
  prop.test(x = sum(fib.sig), n = Ntests, p = 0.01) # sig (.237> .1)
    #write out genes for later analysis##
    fib.goi=matrix$Fib_P < 0.01
    fib.goi.fin <- matrix$gene[fib.goi]
    write.csv(fib.goi.fin, "Fib_sig.csv", row.names=FALSE)
  ##number of interaction effects
  sum(popfib.sig) # 563 popxfib effects
  prop.test(x = sum(popfib.sig), n = Ntests, p = 0.01) # sig (0.003 > 1)
    #write out genes for later analysis##
    popfib.goi=matrix$Treat_P < 0.01
    popfib.goi.fin <- matrix$gene[popfib.goi]
    write.csv(popfib.goi.fin, "PopFib_sig.csv", row.names=FALSE)
    
    
##Post hoc testing##
    
    ##to start we will use emmeans package to do post-hoc tests on population level differences for fibrotic genes of interest##
    #initial process##
    readcount <- readcount[,-1]
    geneIDs <- names(readcount)
    readcount=sapply(readcount,as.numeric)
    metadata$Fib_L_count=sapply(metadata$Fib_L_count,as.numeric)
    metadata$Treatment<- relevel(factor(metadata$Treatment), ref = "PBS")
    metadata$Population<- relevel(factor(metadata$Population), ref = "SAY")
    
    ##models- run individually for each gene of interest with posthoc tests
      ##ada
      gene_i="ENSGACT00000015197"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##RSL sig different##
      #contrast  estimate      SE  df z.ratio p.value
      #SAY - GOS   0.0147 0.00952 Inf   1.548  0.2685
      #SAY - RSL  -0.0485 0.00966 Inf  -5.017  <.0001
      #GOS - RSL  -0.0632 0.00973 Inf  -6.493  <.0001
      
      ##aimp1
      gene_i="ENSGACT00000023165"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##All sig different##
      #contrast  estimate      SE  df z.ratio p.value
      #SAY - GOS   0.0417 0.00911 Inf   4.579  <.0001
      #SAY - RSL  -0.0589 0.00915 Inf  -6.440  <.0001
      #GOS - RSL  -0.1006 0.00928 Inf -10.843  <.0001
      
      
      ##anxa1a
      gene_i="ENSGACT00000015295"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##All sig different##
      #contrast  estimate      SE  df z.ratio p.value
      #SAY - GOS   0.0303 0.00827 Inf   3.662  0.0007
      #SAY - RSL  -0.0517 0.00835 Inf  -6.190  <.0001
      #GOS - RSL  -0.0820 0.00845 Inf  -9.701  <.0001
      
      ##anxa4
      gene_i="ENSGACT00000022600"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##all sig different##
      #contrast  estimate      SE  df z.ratio p.value
      #SAY - GOS   0.0264 0.00898 Inf   2.938  0.0093
      #SAY - RSL  -0.0617 0.00905 Inf  -6.816  <.0001
      #GOS - RSL  -0.0881 0.00915 Inf  -9.628  <.0001
      
      ##BCL11B
      gene_i="ENSGACT00000017387"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##RSL sig different##
      #contrast  estimate      SE  df z.ratio p.value
      #SAY - GOS   0.0401 0.0346 Inf   1.159  0.4776
      #SAY - RSL  -0.1350 0.0337 Inf  -4.003  0.0002
      #GOS - RSL  -0.1751 0.0344 Inf  -5.088  <.0001
      
      
      ##BTLA
      gene_i="ENSGACT00000019304"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##GOS sig different##
      #contrast  estimate      SE  df z.ratio p.value
      #SAY - GOS   0.1406 0.0261 Inf   5.395  <.0001
      #SAY - RSL  -0.0465 0.0254 Inf  -1.835  0.1582
      #GOS - RSL  -0.1871 0.0266 Inf  -7.046  <.0001
      
      
      ##ccdc22  
      gene_i="ENSGACT00000017674"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##all sig different##
      # contrast  estimate      SE  df z.ratio p.value
      #SAY - GOS   0.0356 0.00862 Inf   4.129  0.0001
      #SAY - RSL  -0.0859 0.00862 Inf  -9.963  <.0001
      #GOS - RSL  -0.1215 0.00873 Inf -13.913  <.0001
      
      
      ##CCR6  
      gene_i="ENSGACT00000018265"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##all sig different##
      # contrast  estimate      SE  df z.ratio p.value
      # SAY - GOS  0.00965 0.0489 Inf   0.198  0.9787
      #SAY - RSL  0.31421 0.0543 Inf   5.786  <.0001
      #GOS - RSL  0.30456 0.0548 Inf   5.556  <.0001
      
      
      ##cd40lg  
      gene_i="ENSGACT00000022829"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##all sig different##
      # contrast  estimate      SE  df z.ratio p.value
      #SAY - GOS   -0.118 0.0396 Inf  -2.985  0.0080
      #SAY - RSL    0.147 0.0431 Inf   3.409  0.0019
      #GOS - RSL    0.265 0.0425 Inf   6.233  <.0001
      
      
      ##cxcl12a  
      gene_i="ENSGACT00000009199"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##RSL sig different##
      # contrast  estimate      SE  df z.ratio p.value
      #SAY - GOS  -0.0246 0.0244 Inf  -1.011  0.5701
      #SAY - RSL  -0.2963 0.0234 Inf -12.688  <.0001
      #GOS - RSL  -0.2717 0.0233 Inf -11.643  <.0001
      
      #cxcl19  
      gene_i="ENSGACT00000025278"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##all sig different##
      # contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS   0.0413 0.0124 Inf   3.324  0.0026
      #SAY - RSL  -0.0558 0.0124 Inf  -4.492  <.0001
      #GOS - RSL  -0.0971 0.0126 Inf  -7.684  <.0001
      
      #cyldl
      gene_i="ENSGACT00000007096"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##GOS sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS   0.0530 0.0151 Inf   3.499  0.0014
      #SAY - RSL  -0.0198 0.0153 Inf  -1.291  0.4000
      #GOS - RSL  -0.0728 0.0156 Inf  -4.670  <.0001
      
      #gpnmb
      gene_i="ENSGACT00000009231"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##GOS sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS   0.0615 0.0277 Inf   2.216  0.0684
      #SAY - RSL  -0.1574 0.0268 Inf  -5.873  <.0001
      #GOS - RSL  -0.2189 0.0275 Inf  -7.958  <.0001
      
      
      #grb2b
      gene_i="ENSGACT00000009170"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##GOS sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS  0.00364 0.00766 Inf   0.475  0.8830
      #SAY - RSL -0.06746 0.00775 Inf  -8.706  <.0001
      #GOS - RSL -0.07109 0.00779 Inf  -9.126  <.0001
      
      
      #hmgb2a
      gene_i="ENSGACT00000021812"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##RSL sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS   0.0328 0.00996 Inf   3.297  0.0028
      #SAY - RSL  -0.0414 0.01010 Inf  -4.116  0.0001
      #GOS - RSL  -0.0743 0.01020 Inf  -7.287  <.0001
      
      
      #HMGB1
      gene_i="ENSGACT00000027215"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##RSL sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS  0.00454 0.0192 Inf   0.236  0.9698
      #SAY - RSL -0.25147 0.0185 Inf -13.577  <.0001
      #GOS - RSL -0.25600 0.0187 Inf -13.708  <.0001
      
      
      #iqgap1
      gene_i="ENSGACT00000023086"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##RSL sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS   0.0105 0.0119 Inf   0.882  0.6515
      #SAY - RSL   0.1923 0.0129 Inf  14.915  <.0001
      #GOS - RSL   0.1819 0.0130 Inf  14.025  <.0001
      
      
      #irf4a
      gene_i="ENSGACT00000021776"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##RSL sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS  -0.0476 0.0359 Inf  -1.326  0.3809
      #SAY - RSL  -0.3359 0.0343 Inf  -9.801  <.0001
      #GOS - RSL  -0.2883 0.0339 Inf  -8.494  <.0001
      
      
      ##klf2b
      gene_i="ENSGACT00000022354"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##RSL sig different##
      # contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS  0.00239 0.00989 Inf   0.242  0.9683
      #SAY - RSL -0.07833 0.00998 Inf  -7.846  <.0001
      #GOS - RSL -0.08072 0.01000 Inf  -8.048  <.0001
      
      #LAMP1
      gene_i="ENSGACT00000020286"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##RSL sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS   0.0053 0.0102 Inf   0.520  0.8617
      #SAY - RSL  -0.1101 0.0102 Inf -10.797  <.0001
      #GOS - RSL  -0.1154 0.0103 Inf -11.255  <.0001
      
      
      #lck
      gene_i="ENSGACT00000009597"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##RSL sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS   0.0126 0.0194 Inf   0.651  0.7920
      #SAY - RSL  -0.0694 0.0195 Inf  -3.558  0.0011
      #GOS - RSL  -0.0821 0.0197 Inf  -4.169  0.0001
      
      
      #MAF
      gene_i="ENSGACT00000021127"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##RSL sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS   0.0360 0.0282 Inf   1.278  0.4075
      #SAY - RSL  -0.0994 0.0278 Inf  -3.579  0.0010
      #GOS - RSL  -0.1354 0.0283 Inf  -4.784  <.0001
      
      
      #mrc1b
      gene_i="ENSGACT00000022556"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##All sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS   0.0671 0.0157 Inf   4.275  0.0001
      #SAY - RSL  -0.0463 0.0157 Inf  -2.945  0.0091
      #GOS - RSL  -0.1134 0.0160 Inf  -7.071  <.0001
      
      
      #pikfyve
      gene_i="ENSGACT00000011654"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##GOS sig different##
      #contrast  estimate     SE  df z.ratio p.value
      # SAY - GOS  0.03236 0.00862 Inf   3.754  0.0005
      #SAY - RSL  0.00878 0.00885 Inf   0.993  0.5815
      #GOS - RSL -0.02358 0.00895 Inf  -2.635  0.0229
      
      
      #prkd2
      gene_i="ENSGACT00000012461"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##SAY sig different##
      ##contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS  -0.2211 0.0429 Inf  -5.159  <.0001
      #SAY - RSL  -0.1778 0.0444 Inf  -4.007  0.0002
      #GOS - RSL   0.0433 0.0421 Inf   1.028  0.5592
      
      
      #rbx1
      gene_i="ENSGACT00000010177"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##all sig different##
      #contrast  estimate      SE  df z.ratio p.value
      #SAY - GOS   0.0409 0.00836 Inf   4.901  <.0001
      #SAY - RSL  -0.0728 0.00836 Inf  -8.704  <.0001
      #GOS - RSL  -0.1137 0.00848 Inf -13.407  <.0001
      
      
      #spi1b
      gene_i="ENSGACT00000020522"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##All sig different##
      # contrast  estimate      SE  df z.ratio p.value
      #SAY - GOS   0.0211 0.00890 Inf   2.365  0.0473
      #SAY - RSL  -0.0762 0.00896 Inf  -8.504  <.0001
      #GOS - RSL  -0.0972 0.00904 Inf -10.757  <.0001
      
      #tnfaip8l2b
      gene_i="ENSGACT00000017961"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##RSL sig different##
      # contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS   0.0398 0.0203 Inf   1.961  0.1220
      #SAY - RSL  -0.0898 0.0203 Inf  -4.420  <.0001
      #GOS - RSL  -0.1297 0.0206 Inf  -6.304  <.0001
      
      
      
      #ufd1l
      gene_i="ENSGACT00000005239"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##All sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS   0.0299 0.00831 Inf   3.595  0.0009
      #SAY - RSL  -0.0633 0.00836 Inf  -7.570  <.0001
      #GOS - RSL  -0.0932 0.00846 Inf -11.014  <.0001
      
      
      #zap70
      gene_i="ENSGACT00000010195"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##SAY sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS  -0.1359 0.0355 Inf  -3.826  0.0004
      #SAY - RSL  -0.1707 0.0359 Inf  -4.749  <.0001
      #GOS - RSL  -0.0348 0.0351 Inf  -0.992  0.5820
      
      
      #zbtb1
      gene_i="ENSGACT00000018236"
      readcount <- as.data.frame(readcount)
      geneID <- geneIDs[gene_i]
      y <- readcount[,gene_i]
      Totalreads <- rowSums(readcount)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      model <- glm(Y ~ metadata$Population * metadata$Fib_L_count, family = "quasibinomial")
      em <- emmeans(model, "Population")
      contrast(em, "pairwise", adjust = "Tukey") ##SAY sig different##
      #contrast  estimate     SE  df z.ratio p.value
      #SAY - GOS   0.0801 0.0614 Inf   1.304  0.3926
      #SAY - RSL  -0.1778 0.0588 Inf  -3.024  0.0071
      #GOS - RSL  -0.2580 0.0602 Inf  -4.285  0.0001
    
    
    
    
    
    
    ######Next we want to dig deeper into those with significant interaction terms##
    ##we will run population specific GLMs to test for fibrosis dependence of each gene within each population##
      ### pull relevant data and relevel
      readcount.SAY <- readcount[metadata$Population %in% c("SAY"),]
      metadata.SAY <- metadata[metadata$Population %in% c("SAY"),]
      readcount.SAY<- as.data.frame(readcount.SAY)
      readcount.SAY <- readcount.SAY[,-1]
      geneIDs <- names(readcount.SAY)
      readcount.SAY=sapply(readcount.SAY,as.numeric)
      
      readcount.GOS <- readcount[metadata$Population %in% c("GOS"),]
      metadata.GOS <- metadata[metadata$Population %in% c("GOS"),]
      readcount.GOS<- as.data.frame(readcount.GOS)
      readcount.GOS <- readcount.GOS[,-1]
      geneIDs <- names(readcount.GOS)
      readcount.GOS=sapply(readcount.GOS,as.numeric)
      
      readcount.RSL <- readcount[metadata$Population %in% c("RSL"),]
      metadata.RSL <- metadata[metadata$Population %in% c("RSL"),]
      readcount.RSL<- as.data.frame(readcount.RSL)
      readcount.RSL <- readcount.RSL[,-1]
      geneIDs <- names(readcount.RSL)
      readcount.RSL=sapply(readcount.RSL,as.numeric)
      
      
      ##now run your models, one gene and population at a time##
        ##cyldl: + in GOS/SAY; ns in RSL
        gene_i="ENSGACT00000007096"
        
        readcount.SAY <- as.data.frame(readcount.SAY)
        geneID <- geneIDs[gene_i]
        y <- readcount.SAY[,gene_i]
        Totalreads <- rowSums(readcount.SAY)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.SAY$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p<0.001
        model$coefficients #0.0357
        
        readcount.GOS <- as.data.frame(readcount.GOS)
        geneID <- geneIDs[gene_i]
        y <- readcount.GOS[,gene_i]
        Totalreads <- rowSums(readcount.GOS)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.GOS$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p<0.001
        model$coefficients #0.0347
        
        readcount.RSL <- as.data.frame(readcount.RSL)
        geneID <- geneIDs[gene_i]
        y <- readcount.RSL[,gene_i]
        Totalreads <- rowSums(readcount.RSL)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.RSL$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p=0.901
        model$coefficients #-0.0009
        
        
        ##irf4a: - in GOS/SAY; ns in RSL
        gene_i="ENSGACT00000021776"
        
        readcount.SAY <- as.data.frame(readcount.SAY)
        geneID <- geneIDs[gene_i]
        y <- readcount.SAY[,gene_i]
        Totalreads <- rowSums(readcount.SAY)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.SAY$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p<0.001
        model$coefficients #-0.101
        
        readcount.GOS <- as.data.frame(readcount.GOS)
        geneID <- geneIDs[gene_i]
        y <- readcount.GOS[,gene_i]
        Totalreads <- rowSums(readcount.GOS)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.GOS$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p=0.0384
        model$coefficients #-0.0461
        
        readcount.RSL <- as.data.frame(readcount.RSL)
        geneID <- geneIDs[gene_i]
        y <- readcount.RSL[,gene_i]
        Totalreads <- rowSums(readcount.RSL)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.RSL$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p=0.601
        model$coefficients #-0.00785
        
        
        ##klf2b: + in GOS/SAY; ns in RSL
        gene_i="ENSGACT00000022354"
        
        readcount.SAY <- as.data.frame(readcount.SAY)
        geneID <- geneIDs[gene_i]
        y <- readcount.SAY[,gene_i]
        Totalreads <- rowSums(readcount.SAY)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.SAY$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p=0.0011
        model$coefficients #0.017
        
        readcount.GOS <- as.data.frame(readcount.GOS)
        geneID <- geneIDs[gene_i]
        y <- readcount.GOS[,gene_i]
        Totalreads <- rowSums(readcount.GOS)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.GOS$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p<0.001
        model$coefficients #0.0258
        
        readcount.RSL <- as.data.frame(readcount.RSL)
        geneID <- geneIDs[gene_i]
        y <- readcount.RSL[,gene_i]
        Totalreads <- rowSums(readcount.RSL)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.RSL$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p=0.702
        model$coefficients #0.00188
        
        
        ##mrc1b: + in GOS only
        gene_i="ENSGACT00000022556"
        
        readcount.SAY <- as.data.frame(readcount.SAY)
        geneID <- geneIDs[gene_i]
        y <- readcount.SAY[,gene_i]
        Totalreads <- rowSums(readcount.SAY)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.SAY$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p=0.0694
        model$coefficients #0.0131
        
        readcount.GOS <- as.data.frame(readcount.GOS)
        geneID <- geneIDs[gene_i]
        y <- readcount.GOS[,gene_i]
        Totalreads <- rowSums(readcount.GOS)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.GOS$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p<0.001
        model$coefficients #0.0445
        
        readcount.RSL <- as.data.frame(readcount.RSL)
        geneID <- geneIDs[gene_i]
        y <- readcount.RSL[,gene_i]
        Totalreads <- rowSums(readcount.RSL)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.RSL$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p=0.7301
        model$coefficients #0.00282
        
        
        ##prkd2: - in SAY only
        gene_i="ENSGACT00000012461"
        
        readcount.SAY <- as.data.frame(readcount.SAY)
        geneID <- geneIDs[gene_i]
        y <- readcount.SAY[,gene_i]
        Totalreads <- rowSums(readcount.SAY)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.SAY$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p<0.001
        model$coefficients #-0.0916
        
        readcount.GOS <- as.data.frame(readcount.GOS)
        geneID <- geneIDs[gene_i]
        y <- readcount.GOS[,gene_i]
        Totalreads <- rowSums(readcount.GOS)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.GOS$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p=0.526
        model$coefficients #-0.0150
        
        readcount.RSL <- as.data.frame(readcount.RSL)
        geneID <- geneIDs[gene_i]
        y <- readcount.RSL[,gene_i]
        Totalreads <- rowSums(readcount.RSL)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.RSL$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p=0.664
        model$coefficients #0.00923
        
        
        
        ##rbx1: + in GOS/SAY only
        gene_i="ENSGACT00000010177"
        
        readcount.SAY <- as.data.frame(readcount.SAY)
        geneID <- geneIDs[gene_i]
        y <- readcount.SAY[,gene_i]
        Totalreads <- rowSums(readcount.SAY)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.SAY$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p=0.0163
        model$coefficients #0.01034
        
        readcount.GOS <- as.data.frame(readcount.GOS)
        geneID <- geneIDs[gene_i]
        y <- readcount.GOS[,gene_i]
        Totalreads <- rowSums(readcount.GOS)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.GOS$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p<0.001
        model$coefficients #0.0120
        
        readcount.RSL <- as.data.frame(readcount.RSL)
        geneID <- geneIDs[gene_i]
        y <- readcount.RSL[,gene_i]
        Totalreads <- rowSums(readcount.RSL)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.RSL$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p=0.9439
        model$coefficients #-0.0003
        
        
        
        ##spi1b: + in GOS/SAY only
        gene_i="ENSGACT00000020522"
        
        readcount.SAY <- as.data.frame(readcount.SAY)
        geneID <- geneIDs[gene_i]
        y <- readcount.SAY[,gene_i]
        Totalreads <- rowSums(readcount.SAY)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.SAY$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p<0.001
        model$coefficients #0.0188
        
        readcount.GOS <- as.data.frame(readcount.GOS)
        geneID <- geneIDs[gene_i]
        y <- readcount.GOS[,gene_i]
        Totalreads <- rowSums(readcount.GOS)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.GOS$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p<0.001
        model$coefficients #0.0273
        
        readcount.RSL <- as.data.frame(readcount.RSL)
        geneID <- geneIDs[gene_i]
        y <- readcount.RSL[,gene_i]
        Totalreads <- rowSums(readcount.RSL)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.RSL$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p=0.113
        model$coefficients #0.00540
        
        
        ##tnfaip8l2b: + in GOS/SAY; ns in RSL
        gene_i="ENSGACT00000017961"
        
        readcount.SAY <- as.data.frame(readcount.SAY)
        geneID <- geneIDs[gene_i]
        y <- readcount.SAY[,gene_i]
        Totalreads <- rowSums(readcount.SAY)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.SAY$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p-0.00926
        model$coefficients #0.0284
        
        readcount.GOS <- as.data.frame(readcount.GOS)
        geneID <- geneIDs[gene_i]
        y <- readcount.GOS[,gene_i]
        Totalreads <- rowSums(readcount.GOS)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.GOS$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p<0.001
        model$coefficients #-0.0649
        
        readcount.RSL <- as.data.frame(readcount.RSL)
        geneID <- geneIDs[gene_i]
        y <- readcount.RSL[,gene_i]
        Totalreads <- rowSums(readcount.RSL)
        Y <- cbind(focalgene = y, otherreads = Totalreads - y)
        Y <- as.matrix(Y)
        Yprop <- Y[,1] / rowSums(Y)
        mean_expr <- mean(  Yprop, na.rm = T)
        model <- glm(Y ~ metadata.RSL$Fib_L_count, family = "quasibinomial")
        anova(model, test = "Chisq") #p=0.0507
        model$coefficients #0.0171      
        