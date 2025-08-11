##This is the code to take our normalized head kidney reads and run our glm models on Alum and Cestode Protein Individually##
##It also includes steps to statistically evaluate overrepresentation of DEGs and compare DEGs across treatments##
##NOTE- these are three-way models (all pops, time points) split by treatment##
##last updated Aug 5th, 2025 by LEF##

##start by loading in necessary pacakges##
library(tibble)
library(data.table)
library(tidyverse)

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
  metadata <- read.csv("ExperimentalDesign.csv")
  metadata <- na.omit(metadata)
  head(metadata)
  dim(readcount)
  dim(metadata)
  
  }

##read in gene annotation data in case you want it##
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

##remove day 90 from these analyses##
  {
    # Omit day 90
  readcount <- readcount[metadata$DPI < 90,]
  metadata <- metadata[metadata$DPI < 90,]
  }

######## STEP 1- run models

##start with alum, the following code isolaltes just alum and control (PBS), and sets appropriate reference levels##
  {
    ### Focusing on a particular treatment, alum
  readcount.al <- readcount[metadata$Treatment %in% c("Alum", "PBS"),]
  metadata.al <- metadata[metadata$Treatment %in% c("Alum", "PBS"),]
  metadata.al$treatmentColor <- as.numeric(factor(metadata.al$Treatment))*-1 + 3
  metadata.al$Treatment<- relevel(factor(metadata.al$Treatment), ref = "PBS")
  metadata.al$Population<- relevel(factor(metadata.al$Population), ref = "SAY")
  readcount.al<- as.data.frame(readcount.al)
  }

###### GLM of alum effects to identify instances with a main effect of treatment, time, population*treatment, etc
####### SKIP THIS IF RERUNNING and instead reload output
  {
  ##start by finding the actual quartile to use##
  #prop.matrix <- matrix(data = NA, nrow = ncol(readcount.al), ncol = 2)
  #prop.matrix <- as.data.frame(prop.matrix)
  #readcount.al <- readcount.al[,-1]
  #geneIDs <- names(readcount.al)
  #readcount.al=sapply(readcount.al,as.numeric)
  
  #for(gene_i in 1:ncol(readcount.al)){
    #readcount.al <- as.data.frame(readcount.al)
    #geneID <- geneIDs[gene_i]
    #y <- readcount.al[,gene_i]
    #Totalreads <- rowSums(readcount.al)
    #Y <- cbind(focalgene = y, otherreads = Totalreads - y)
    #Y <- as.matrix(Y)
    #Yprop <- Y[,1] / rowSums(Y)
    #mean_expr <- mean(  Yprop, na.rm = T)
    #prop.matrix[gene_i,]= c(geneID, mean_expr)}
  #prop.matrix$V2=lapply(prop.matrix$V2, as.numeric)
  #quantile(unlist(prop.matrix$V2), probs = 0.5, na.rm = TRUE)
  
  ##then do the rest, actually running the model##
    alum.matrix <- matrix(data = NA, nrow = ncol(readcount.al), ncol = 22)
    alum.matrix <- as.data.frame(alum.matrix)
    readcount.al <- readcount.al[,-1]
    geneIDs <- names(readcount.al)
    readcount.al=sapply(readcount.al,as.numeric)
  
  for(gene_i in 1:ncol(readcount.al)){
    readcount.al <- as.data.frame(readcount.al)
    geneID <- geneIDs[gene_i]
    y <- readcount.al[,gene_i]
    Totalreads <- rowSums(readcount.al)
    Y <- cbind(focalgene = y, otherreads = Totalreads - y)
    Y <- as.matrix(Y)
    Yprop <- Y[,1] / rowSums(Y)
    mean_expr <- mean(  Yprop, na.rm = T)
    if(mean_expr > 2.487893e-05  ){ #50% quantile mean read count threshold
       model <- glm(Y ~ metadata.al$Population * metadata.al$Treatment * metadata.al$DPI, family = "quasibinomial")
      results <- anova(model, test = "Chisq")
      alum.matrix[gene_i,] <- c( geneID, mean_expr, results$`Pr(>Chi)`,model$coefficients)
    }  else{
      alum.matrix[gene_i,] <- c( geneID, mean_expr, NA, NA, NA, NA,NA, NA, NA,  NA, NA, NA, NA, NA,NA, NA, NA,NA, NA, NA, NA,NA)
    }
  }
  
  ## rename columns of matrix
  names(alum.matrix) <- c("gene", "meanreadprop", "NA", "Pop_P", "Treat_P", "DPI_P", "Pop_Trt_P", "Pop_DPI_P", "Treat_DPI_P", "Pop_Trt_DPI_P","Int_C",
                          "Pop_GvS_C","Pop_RvS_C","Treat_C", "DPI_C", "GvS:Treat_C", "RvS:Treat_C", "GvS:DPI_C", "RvS:DPI_C", "Treat:DPI_C","GvS:Treat:DPI_C", "RvS:Treat:DPI_C")
  
  ##check and format output matrix##
  head(alum.matrix)
  alum.matrix <- alum.matrix[order(alum.matrix$Pop_Trt_DPI_P, decreasing = F),]
  alum.matrix[1:20,]
  alum.matrix <- alum.matrix[,-3]
  alum.matrix <- alum.matrix[,-10]
  dim(na.omit(alum.matrix))
  alum.matrix.final=na.omit(alum.matrix)
  #dim(alum.matrix)
  
  ##some outputs##
  write.csv(alum.matrix.final, "Alum_GLM_matrix_results.csv", row.names=FALSE)
  write.csv(alum.matrix, "Alum_GLM_matrix_results_everything.csv", row.names=FALSE)
  }




##next do cestode protein, the following code isolaltes just cp and control (PBS), and sets appropriate reference levels##
  ### Focusing on a particular treatment, CP
  readcount.cp <- readcount[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  metadata.cp <- metadata[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  metadata.cp$treatmentColor <- as.numeric(factor(metadata.cp$Treatment))
  metadata.cp$Treatment<- relevel(factor(metadata.cp$Treatment), ref = "PBS")
  metadata.cp$Population<- relevel(factor(metadata.cp$Population), ref = "SAY")
  readcount.cp <- readcount.cp[,-1]
  geneIDs <- names(readcount.cp)
  readcount.cp=sapply(readcount.cp,as.numeric)

  ###### GLM of alum effects to identify instances with a main effect of treatment, time, population*treatment, etc
  ####### SKIP THIS IF RERUNNING and instead reload output
  {##start by finding the actual quartile to use##
    #prop.matrix <- matrix(data = NA, nrow = ncol(readcount.cp), ncol = 2)
    #prop.matrix <- as.data.frame(prop.matrix)
    #readcount.cp <- readcount.cp[,-1]
    #geneIDs <- names(readcount.cp)
    #readcount.cp=sapply(readcount.cp,as.numeric)
    
    #for(gene_i in 1:ncol(readcount.cp)){
      #readcount.cp <- as.data.frame(readcount.cp)
      #geneID <- geneIDs[gene_i]
      #y <- readcount.cp[,gene_i]
      #Totalreads <- rowSums(readcount.cp)
      #Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      #Y <- as.matrix(Y)
      #Yprop <- Y[,1] / rowSums(Y)
      #mean_expr <- mean(  Yprop, na.rm = T)
      #prop.matrix[gene_i,]= c(geneID, mean_expr)}
    #prop.matrix$V2=lapply(prop.matrix$V2, as.numeric)
    #quantile(unlist(prop.matrix$V2), probs = 0.5, na.rm=TRUE)
    
  ###### then run the actual GLM
  cp.matrix <- matrix(data = NA, nrow = ncol(readcount.cp), ncol = 22)
  cp.matrix <- as.data.frame(cp.matrix)
  readcount.cp <- readcount.cp[,-1]
  geneIDs <- names(readcount.cp)
  readcount.cp=sapply(readcount.cp,as.numeric)
  
  for(gene_i in 1:ncol(readcount.cp)){
    readcount.cp <- as.data.frame(readcount.cp)
    geneID <- "ENSGACT00000020986"
    y <- readcount.cp[,gene_i]
    Totalreads <- rowSums(readcount.cp)
    Y <- cbind(focalgene = y, otherreads = Totalreads - y)
    Y <- as.matrix(Y)
    Yprop <- Y[,1] / rowSums(Y)
    mean_expr <- mean(  Yprop, na.rm = T)
    if(mean_expr > 2.487893e-05  ){ #50% quantile mean read count threshold
      model <- glm(Y ~ metadata.cp$Population * metadata.cp$Treatment * metadata.cp$DPI, family = "quasibinomial")
      results <- anova(model, test = "Chisq")
      cp.matrix[gene_i,] <- c( geneID, mean_expr, results$`Pr(>Chi)`,model$coefficients)
    }  else{
      cp.matrix[gene_i,] <- c( geneID, mean_expr, NA, NA, NA, NA,NA, NA, NA,  NA, NA, NA, NA, NA,NA, NA, NA,NA, NA, NA, NA,NA)
    }
  }
  ##rename matrix columns
  names(cp.matrix) <- c("gene", "meanreadprop", "NA", "Pop_P", "Treat_P", "DPI_P", "Pop_Trt_P", "Pop_DPI_P", "Treat_DPI_P", "Pop_Trt_DPI_P","Int_C",
                          "Pop_GvS_C","Pop_RvS_C","Treat_C", "DPI_C", "GvS:Treat_C", "RvS:Treat_C", "GvS:DPI_C", "RvS:DPI_C", "Treat:DPI_C","GvS:Treat:DPI_C", "RvS:Treat:DPI_C")
  
  ##format and output matrices##
  cp.matrix <- cp.matrix[order(cp.matrix$Pop_Trt_DPI_P, decreasing = F),]
  cp.matrix[1:20,]
  cp.matrix <- cp.matrix[,-3]
  cp.matrix <- cp.matrix[,-10]
  dim(na.omit(cp.matrix))
  dim(cp.matrix)
  cp.matrix.final=na.omit(cp.matrix)
  
  write.csv(cp.matrix, "CestodeProt_GLM_matrix_results_all.csv", row.names=FALSE)
  write.csv(cp.matrix.final, "CestodeProt_GLM_matrix_results.csv", row.names=FALSE)
  }


##restart HERE
###STEP 2- ANALYSIS AND COMPARISON OF DEGS##
  #####reload results##
  cp.matrix <- read.csv("CestodeProt_GLM_matrix_results_fixed.csv")
  alum.matrix <- read.csv( "Alum_GLM_matrix_results_fixed.csv")

##First question##
  ### Are there more P < 0.01 than we expect?
  ### Building a venn diagram
  cp.matrix2<- na.omit(cp.matrix)
  alum.matrix2 <- na.omit( alum.matrix)
  alpha <- 0.01 ##can change if interested; we like conservative##
  merged.matrix2 <- merge(alum.matrix2, cp.matrix2, by = "gene", all = T)
  merged.matrix2 <- na.omit(merged.matrix2)
  
  # Number of significant for each model term for each model##
  {
    Ntests.cp <- nrow(cp.matrix2)
    Ntests.al <- nrow(alum.matrix2)
  cp.sig <- cp.matrix2$Treat_P < alpha
  a.sig <- alum.matrix2$Treat_P < alpha
  cp_geno.sig <- cp.matrix2$Pop_Trt_P < alpha
  a_geno.sig <- alum.matrix2$Pop_Trt_P < alpha
  cp_day.sig <- cp.matrix2$Treat_DPI_P < alpha
  a_day.sig <- alum.matrix2$Treat_DPI_P < alpha
  cp_3w.sig <- cp.matrix2$Pop_Trt_DPI_P < alpha
  a_3w.sig <- alum.matrix2$Pop_Trt_DPI_P < alpha
   
  pop.day.cpsig <- cp.matrix2$Pop_P < alpha
  pop.cpsig <- cp.matrix2$Pop_DPI_P < alpha
  day.cpsig <- cp.matrix2$DPI_P < alpha
  pop.day.alsig <- alum.matrix2$Pop_P < alpha
  pop.alsig <- alum.matrix2$Pop_DPI_P < alpha
  day.alsig <- alum.matrix2$DPI_P < alpha
  
  }
  
  # For venn diagram but not key otherwise
  sum(pop.cpsig) # 1806
  sum(day.cpsig) #3063
  sum(pop.day.cpsig)  #10087
  sum(pop.alsig) # 1590
  sum(day.alsig)  #3395
  sum(pop.day.alsig) # 9900

  # number of main effects of treatment
    ##alum##
    sum(a.sig) # 1365 alum effects
    ##more than expected?##
    prop.test(x = sum(a.sig), n = Ntests.al, p = 0.01) # sig (.105 > .01)
    #write out genes for later analysis##
      alum.goi=alum.matrix2$Treat_P < 0.01
      alum.goi.fin <- alum.matrix2$gene[alum.goi]
      write.csv(alum.goi.fin, "Alum_treat_sig.csv", row.names=FALSE)  
    ##cp##
    sum(cp.sig) # 76 cestode protein effects
    ##more than expected??##
    prop.test(x = sum(cp.sig), n = Ntests.cp, p = 0.01) # sig (.00578 < .01)
      #write out genes for later analysis##
      cp.goi=cp.matrix2$Treat_P < 0.01
      cp.goi.fin <- cp.matrix2$gene[cp.goi]
      write.csv(cp.goi.fin, "CP_treat_sig.csv", row.names=FALSE)



# number of effects of treatment*time
  ##alum##
  sum(a_day.sig) # 195 alum effects
  ##more than expected?
  prop.test(x = sum(a_day.sig), n = Ntests.al, p = 0.01) # sig (0.015 > 0.01)
    #write out genes for later analysis##
    alum.day.goi=alum.matrix2$Treat_DPI_P < 0.01
    alum.day.goi.fin <- alum.matrix2$gene[alum.day.goi]
    write.csv(alum.day.goi.fin, "alum_treatxDay_sig.csv", row.names=FALSE)
  ##cp##
  sum(cp_day.sig) # 378 cestode protein effects
  ##more than expected?##
  prop.test(x = sum(cp_day.sig), n = Ntests.cp, p = 0.01) # sig (0.0287 > 0.01)
    #write out genes for later analysis##
    cp.day.goi=cp.matrix2$Treat_DPI_P < 0.01
    cp.day.goi.fin <- cp.matrix2$gene[cp.day.goi]
    write.csv(cp.day.goi.fin, "CP_treatxDay_sig.csv", row.names=FALSE)


# number of genes showing treatment OR treatment*time interaction effects
  sum(a.sig | a_day.sig)  #1544 genes responded to alum, with or without time effect
  sum(a.sig & a_day.sig) # 16 genes had both main and time-dependent alum effects
  sum(cp.sig | cp_day.sig)  #453 genes responded to protein (with or without time effect)
  sum(cp.sig & cp_day.sig)  # 1 genes had both main and time-dependent protein effects

# number of effects of genotype*treatment
  sum(cp_geno.sig) # 82 cestode protein effects
  prop.test(x = sum(cp_geno.sig), n = Ntests.cp, p = 0.01) # sig (0.00623 < 0.01)
    #write out genes for later analysis##
    cp_geno.goi=cp.matrix2$Pop_Trt_P < 0.01
    cp_geno.goi.fin <- cp.matrix2$gene[cp_geno.goi]
    write.csv(cp_geno.goi.fin, "cp_treatxpop_sig.csv", row.names=FALSE)
  sum(a_geno.sig) # 57 alum effects
  prop.test(x = sum(a_geno.sig), n = Ntests.al, p = 0.01) # Significant (0.004 < 0.01)
    #write out genes for later analysis##
    alum_geno.goi=alum.matrix2$Pop_Trt_P < 0.01
    alum_geno.goi.fin <- alum.matrix2$gene[alum_geno.goi]
    write.csv(alum_geno.goi.fin, "alum_treatxpop_sig.csv", row.names=FALSE)

# number of effects of genotype*treatment*time
  sum(cp_3w.sig) # 35 cestode protein effects
  prop.test(x = sum(cp_3w.sig), n = Ntests.cp, p = 0.01) # Significant (.002 < .01)
    #write out genes for later analysis##
    cp_3w.goi=cp.matrix2$Pop_Trt_DPI_P < 0.01
    cp_3w.goi.fin <- cp.matrix2$gene[cp_3w.goi]
    write.csv(cp_3w.goi.fin, "cp_threeway_sig.csv", row.names=FALSE)
  sum(a_3w.sig) # 47 alum effects
  prop.test(x = sum(a_3w.sig), n = Ntests.al, p = 0.01) # Significant (.003 > .01)
    #write out genes for later analysis##
    alum_3w.goi=alum.matrix2$Pop_Trt_DPI_P < 0.01
    alum_3w.goi.fin <- alum.matrix2$gene[alum_3w.goi]
    write.csv(alum_3w.goi.fin, "alum_threeway_sig.csv", row.names=FALSE)
  
# main effects without population interactions
  sum(a.sig &  a_geno.sig == F &  a_3w.sig == F) #1363
  sum(cp.sig &  cp_geno.sig == F &  cp_3w.sig == F) #76
  
##let's consider shared overlap##
# Number shared main effects
  shared.trt <- merged.matrix2$Treat_P.x < 0.01 & merged.matrix2$Treat_P.y < 0.01
  sum(shared.trt)# 13
  #write out genes for later analysis##
  shared.goi <- merged.matrix2$gene[shared.trt]
  write.csv(shared.goi, "Shared_treat_sig.csv", row.names=FALSE)
  ##check proportions##
  al.prop <- 1365/13033
  cp.prop <- 76/13150
  random_expect <- al.prop * cp.prop
  random_expect*13033  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) # ns
 
# Number of shared time*treatment interactions
  shared.trt.time <- merged.matrix2$Treat_DPI_P.x < 0.01 & merged.matrix2$Treat_DPI_P.y < 0.01
  sum(shared.trt.time) #21
  #write out genes for later analysis##
  shared.time.goi <- merged.matrix2$gene[shared.trt.time]
  write.csv(shared.time.goi, "Shared_treatxtime_sig.csv", row.names=FALSE)
  ##check proportions##
  al.prop <- 195/13033
  cp.prop <- 378/13150
  random_expect <- al.prop * cp.prop
  random_expect*13033  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt.time), n = length(shared.trt.time), p = random_expect) #   sig (0.0016 > 0.000430)
  
  # merging alum and cp- further exploration
  sum(shared.trt & shared.trt.time)  # 0 genes had injection and injection*time effects for both alum and cp
  sum(shared.trt | shared.trt.time)  # 34 genes had injection or injection*time effects for both alum and cp
  sum(shared.trt == T & shared.trt.time == F ) # of 13 shared main effects, 13 had no time interaction
  sum(shared.trt == F & shared.trt.time == T ) # of 21 shared time*treatment interactions, 21 had no main effect
  
  
  
# Number shared treatment*genotype effects
  shared.trt <- merged.matrix2$Pop_Trt_P.x < 0.01 & merged.matrix2$Pop_Trt_P.y < 0.01
  sum(shared.trt)# 5
  #write out genes for later analysis##
  shared.pop.goi <- merged.matrix2$gene[shared.trt]
  write.csv(shared.pop.goi, "Shared_treatxpop_sig.csv", row.names=FALSE)
  ##do we deviate from null?
  alum.prop <- 57/13033
  cestode.prop <- 82/13034
  random_expect <- alum.prop * cestode.prop
  random_expect*13033  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) #   sig (0.00387 > 0.0000275)

# Number shared treatment*genotype effects
  shared.3wtrt <- merged.matrix2$Pop_Trt_DPI_P.x < 0.01 & merged.matrix2$Pop_Trt_DPI_P.y < 0.01
  sum(shared.3wtrt)# 1
    #write out genes for later analysis##
    shared.3w.goi <- merged.matrix2$gene[shared.3wtrt]
    write.csv(shared.3w.goi, "Shared_threeway_sig.csv", row.names=FALSE)
  ##do we deviate from null?##
  alum.prop <- 47/13033
  cestode.prop <- 35/13034
  random_expect <- alum.prop * cestode.prop
  random_expect*13033  # number of genes expected to be shared under null
  #do we deviate from expectation?
  prop.test(x = sum(shared.3wtrt), n = length(shared.3wtrt), p = random_expect) #   ns 
  
  
# Number shared alum main vs cp*time
  shared.3wtrt <- merged.matrix2$Treat_P.x < 0.01 & merged.matrix2$Treat_DPI_P.y < 0.01
  sum(shared.3wtrt)# 54
  alum.prop <- 1365/13033
  cestode.prop <- 378/13034
  random_expect <- alum.prop * cestode.prop
  random_expect*13033  # number of genes expected to be shared under null
  #do we deviate from expectation?
  prop.test(x = sum(shared.3wtrt), n = length(shared.3wtrt), p = random_expect) #   sig (0.00418 > 0.00304)

  
  
  
