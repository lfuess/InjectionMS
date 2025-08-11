##This is the code to take our normalized head kidney reads and run our glm models on Alum and Cestode Protein Individually for each population##
##It also includes steps to statistically evaluate overrepresentation of DEGs and compare DEGs across treatments##
##NOTE- these are two-way models (all pops, time points) split by population and treatment##
##last updated Aug 5th, 2025 by LEF##

##load in relevant packages##
library(tibble)
library(data.table)
library(dplyr)

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

##gene info if you want it##
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

##remove day 90
  {
    # Omit day 90
    readcount <- readcount[metadata$DPI < 90,]
    metadata <- metadata[metadata$DPI < 90,]
  }

######## STEP 1- each model##

##Start with SAY
##and Alum specifically

### pull relevant data and relevel
  readcount.al <- readcount[metadata$Treatment %in% c("Alum", "PBS"),]
  metadata.al <- metadata[metadata$Treatment %in% c("Alum", "PBS"),]
  readcount.al.SAY <- readcount.al[metadata.al$Population %in% c("SAY"),]
  metadata.al.SAY <- metadata.al[metadata.al$Population %in% c("SAY"),]
  metadata.al.SAY$Treatment<- relevel(factor(metadata.al.SAY$Treatment), ref = "PBS")
  readcount.al.SAY<- as.data.frame(readcount.al.SAY)

  ####### SKIP THIS IF RERUNNING and instead reload output
  ###### GLM of alum effects to identify instances with a main effect of treatment, time, treatment*time, etc
  {##start by finding the actual quartile to use##
    #prop.matrix <- matrix(data = NA, nrow = ncol(readcount.al), ncol = 2)
    #prop.matrix <- as.data.frame(prop.matrix)
    #readcount.al.SAY <- readcount.al.SAY[,-1]
    #geneIDs <- names(readcount.al.SAY)
    #readcount.al.SAY=sapply(readcount.al.SAY,as.numeric)
    
    #for(gene_i in 1:ncol(readcount.al.SAY)){
    #readcount.al.SAY <- as.data.frame(readcount.al.SAY)
    #geneID <- geneIDs[gene_i]
    #y <- readcount.al.SAY[,gene_i]
    #Totalreads <- rowSums(readcount.al.SAY)
    #Y <- cbind(focalgene = y, otherreads = Totalreads - y)
    #Y <- as.matrix(Y)
    #Yprop <- Y[,1] / rowSums(Y)
    #mean_expr <- mean(  Yprop, na.rm = T)
    #prop.matrix[gene_i,]= c(geneID, mean_expr)}
    #prop.matrix$V2=lapply(prop.matrix$V2, as.numeric)
    #quantile(unlist(prop.matrix$V2), probs = 0.5, na.rm = TRUE)
    
    ##then do the rest, actually running the model##
    alum.matrix <- matrix(data = NA, nrow = ncol(readcount.al.SAY), ncol = 10)
    alum.matrix <- as.data.frame(alum.matrix)
    readcount.al.SAY <- readcount.al.SAY[,-1]
    geneIDs <- names(readcount.al.SAY)
    readcount.al.SAY=sapply(readcount.al.SAY,as.numeric)
    
    for(gene_i in 1:ncol(readcount.al.SAY)){
      readcount.al <- as.data.frame(readcount.al.SAY)
      geneID <- geneIDs[gene_i]
      y <- readcount.al.SAY[,gene_i]
      Totalreads <- rowSums(readcount.al.SAY)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      if(mean_expr > 2.514422e-05    ){ #50% quantile mean read count threshold
        model <- glm(Y ~ metadata.al.SAY$Treatment * metadata.al.SAY$DPI, family = "quasibinomial")
        results <- anova(model, test = "Chisq")
        alum.matrix[gene_i,] <- c( geneID, mean_expr, results$`Pr(>Chi)`,model$coefficients)
      }  else{
        alum.matrix[gene_i,] <- c( geneID, mean_expr, NA, NA, NA, NA,NA, NA, NA,  NA, NA, NA, NA, NA)
      }
    }
    
    ##adjust matrix names
    names(alum.matrix) <- c("gene", "meanreadprop", "NA", "Treat_P", "DPI_P", "Treat_DPI_P", "Int_C",
                            "Treat_C","DPI_C","Treat_DPI_C")
    ##format and output matrix##
    head(alum.matrix)
    alum.matrix <- alum.matrix[order(alum.matrix$Treat_P, decreasing = F),]
    alum.matrix[1:20,]
    alum.matrix <- alum.matrix[,-3]
    alum.matrix <- alum.matrix[,-6]
    dim(na.omit(alum.matrix))
    alum.matrix.final=na.omit(alum.matrix)
    #dim(alum.matrix)
    write.csv(alum.matrix.final, "SAY_Alum_GLM_matrix_results.csv", row.names=FALSE)
    write.csv(alum.matrix, "SAY_Alum_GLM_matrix_results_everything.csv", row.names=FALSE)
  }

##Next we do SAY CP

  ### pull relevant data and relevel
  readcount.cp <- readcount[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  metadata.cp <- metadata[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  readcount.cp.SAY <- readcount.cp[metadata.cp$Population %in% c("SAY"),]
  metadata.cp.SAY <- metadata.cp[metadata.cp$Population %in% c("SAY"),]
  metadata.cp.SAY$Treatment<- relevel(factor(metadata.cp.SAY$Treatment), ref = "PBS")
  readcount.cp.SAY<- as.data.frame(readcount.cp.SAY)
  readcount.cp.SAY <- readcount.cp.SAY[,-1]
  geneIDs <- names(readcount.cp.SAY)
  readcount.cp.SAY=sapply(readcount.cp.SAY,as.numeric)
  
  ####### SKIP THIS IF RERUNNING and instead reload output
  ###### GLM of alum effects to identify instances with a main effect of treatment, time, treatment*time, etc
  {##start by finding the actual quartile to use##
    prop.matrix <- matrix(data = NA, nrow = ncol(readcount.cp.SAY), ncol = 2)
    prop.matrix <- as.data.frame(prop.matrix)
    
    for(gene_i in 1:ncol(readcount.cp.SAY)){
      readcount.cp.SAY <- as.data.frame(readcount.cp.SAY)
      geneID <- geneIDs[gene_i]
      y <- readcount.cp.SAY[,gene_i]
      Totalreads <- rowSums(readcount.cp.SAY)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      prop.matrix[gene_i,]= c(geneID, mean_expr)}
    prop.matrix$V2=lapply(prop.matrix$V2, as.numeric)
    quantile(unlist(prop.matrix$V2), probs = 0.5, na.rm=TRUE)
    
    ##now run actual GLM##
    cp.matrix <- matrix(data = NA, nrow = ncol(readcount.cp.SAY), ncol = 10)
    cp.matrix <- as.data.frame(cp.matrix)
    readcount.cp.SAY <- readcount.cp.SAY[,-1]
    geneIDs <- names(readcount.cp.SAY)
    readcount.cp.SAY=sapply(readcount.cp.SAY,as.numeric)
    
    for(gene_i in 1:ncol(readcount.cp.SAY)){
      readcount.cp <- as.data.frame(readcount.cp.SAY)
      geneID <- geneIDs[gene_i]
      y <- readcount.cp.SAY[,gene_i]
      Totalreads <- rowSums(readcount.cp.SAY)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      if(mean_expr > 2.542707e-05  ){ #50% quantile mean read count threshold
        model <- glm(Y ~ metadata.cp.SAY$Treatment * metadata.cp.SAY$DPI, family = "quasibinomial")
        results <- anova(model, test = "Chisq")
        cp.matrix[gene_i,] <- c( geneID, mean_expr, results$`Pr(>Chi)`,model$coefficients)
      }  else{
        cp.matrix[gene_i,] <- c( geneID, mean_expr, NA, NA, NA, NA,NA, NA, NA,  NA, NA, NA, NA, NA)
      }
    }
    
    ##adjust matrix names
    names(cp.matrix) <- c("gene", "meanreadprop", "NA", "Treat_P", "DPI_P", "Treat_DPI_P", "Int_C",
                          "Treat_C","DPI_C","Treat_DPI_C")
    ##format and output matrices##
    head(cp.matrix)
    cp.matrix <- cp.matrix[order(cp.matrix$Treat_P, decreasing = F),]
    cp.matrix[1:20,]
    cp.matrix <- cp.matrix[,-3]
    cp.matrix <- cp.matrix[,-6]
    dim(na.omit(cp.matrix))
    cp.matrix.final=na.omit(cp.matrix)
    dim(cp.matrix.final)
    write.csv(cp.matrix.final, "SAY_CP_GLM_matrix_results.csv", row.names=FALSE)
    write.csv(cp.matrix, "SAY_CP_GLM_matrix_results_everything.csv", row.names=FALSE)
  }


##Next up is RSL
##First Alum
  ### pull relevant data and relevel
  readcount.al <- readcount[metadata$Treatment %in% c("Alum", "PBS"),]
  metadata.al <- metadata[metadata$Treatment %in% c("Alum", "PBS"),]
  readcount.al.RSL <- readcount.al[metadata.al$Population %in% c("RSL"),]
  metadata.al.RSL <- metadata.al[metadata.al$Population %in% c("RSL"),]
  metadata.al.RSL$Treatment<- relevel(factor(metadata.al.RSL$Treatment), ref = "PBS")
  readcount.al.RSL<- as.data.frame(readcount.al.RSL)

  ####### SKIP THIS IF RERUNNING and instead reload output
  ###### GLM of alum effects to identify instances with a main effect of treatment, time, treatment*time, etc
  {##start by finding the actual quartile to use##
    prop.matrix <- matrix(data = NA, nrow = ncol(readcount.al), ncol = 2)
    prop.matrix <- as.data.frame(prop.matrix)
    readcount.al.RSL <- readcount.al.RSL[,-1]
    geneIDs <- names(readcount.al.RSL)
    readcount.al.RSL=sapply(readcount.al.RSL,as.numeric)
    
    for(gene_i in 1:ncol(readcount.al.RSL)){
    readcount.al.RSL <- as.data.frame(readcount.al.RSL)
    geneID <- geneIDs[gene_i]
    y <- readcount.al.RSL[,gene_i]
    Totalreads <- rowSums(readcount.al.RSL)
    Y <- cbind(focalgene = y, otherreads = Totalreads - y)
    Y <- as.matrix(Y)
    Yprop <- Y[,1] / rowSums(Y)
    mean_expr <- mean(  Yprop, na.rm = T)
    prop.matrix[gene_i,]= c(geneID, mean_expr)}
    prop.matrix$V2=lapply(prop.matrix$V2, as.numeric)
    quantile(unlist(prop.matrix$V2), probs = 0.5, na.rm = TRUE)
    
    ##then actually run the GLM##
    alum.matrix <- matrix(data = NA, nrow = ncol(readcount.al.RSL), ncol = 10)
    alum.matrix <- as.data.frame(alum.matrix)
    readcount.al.RSL <- readcount.al.RSL[,-1]
    geneIDs <- names(readcount.al.RSL)
    readcount.al.RSL=sapply(readcount.al.RSL,as.numeric)
    
    for(gene_i in 1:ncol(readcount.al.RSL)){
      readcount.al <- as.data.frame(readcount.al.RSL)
      geneID <- geneIDs[gene_i]
      y <- readcount.al.RSL[,gene_i]
      Totalreads <- rowSums(readcount.al.RSL)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      if(mean_expr > 2.248659e-05   ){ #50% quantile mean read count threshold
        model <- glm(Y ~ metadata.al.RSL$Treatment * metadata.al.RSL$DPI, family = "quasibinomial")
        results <- anova(model, test = "Chisq")
        alum.matrix[gene_i,] <- c( geneID, mean_expr, results$`Pr(>Chi)`,model$coefficients)
      }  else{
        alum.matrix[gene_i,] <- c( geneID, mean_expr, NA, NA, NA, NA,NA, NA, NA,  NA, NA, NA, NA, NA)
      }
    }
    ##set up column names##
    names(alum.matrix) <- c("gene", "meanreadprop", "NA", "Treat_P", "DPI_P", "Treat_DPI_P", "Int_C",
                            "Treat_C","DPI_C","Treat_DPI_C")
    ##format and output matrix##
    head(alum.matrix)
    alum.matrix <- alum.matrix[order(alum.matrix$Treat_P, decreasing = F),]
    alum.matrix[1:20,]
    alum.matrix <- alum.matrix[,-3]
    alum.matrix <- alum.matrix[,-6]
    dim(na.omit(alum.matrix))
    alum.matrix.final=na.omit(alum.matrix)
    #dim(alum.matrix)
    write.csv(alum.matrix.final, "RSL_Alum_GLM_matrix_results.csv", row.names=FALSE)
    write.csv(alum.matrix, "RSL_Alum_GLM_matrix_results_everything.csv", row.names=FALSE)
  }
  
##Next is RSL, CP

  ### pull relevant data and relevel
  readcount.cp <- readcount[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  metadata.cp <- metadata[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  readcount.cp.RSL <- readcount.cp[metadata.cp$Population %in% c("RSL"),]
  metadata.cp.RSL <- metadata.cp[metadata.cp$Population %in% c("RSL"),]
  metadata.cp.RSL$Treatment<- relevel(factor(metadata.cp.RSL$Treatment), ref = "PBS")
  readcount.cp.RSL<- as.data.frame(readcount.cp.RSL)
  readcount.cp.RSL <- readcount.cp.RSL[,-1]
  geneIDs <- names(readcount.cp.RSL)
  readcount.cp.RSL=sapply(readcount.cp.RSL,as.numeric)
  
  ###### SKIP THIS IF RERUNNING< skip ahead and reload output
  ###### GLM of alum effects to identify instances with a main effect of treatment, time, time*treatment, etc
  { ##start by finding the actual quartile to use##
    #prop.matrix <- matrix(data = NA, nrow = ncol(readcount.cp.RSL), ncol = 2)
    #prop.matrix <- as.data.frame(prop.matrix)
    
    #for(gene_i in 1:ncol(readcount.cp.RSL)){
      #readcount.cp.RSL <- as.data.frame(readcount.cp.RSL)
      #geneID <- geneIDs[gene_i]
      #y <- readcount.cp.RSL[,gene_i]
      #Totalreads <- rowSums(readcount.cp.RSL)
      #Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      #Y <- as.matrix(Y)
      #Yprop <- Y[,1] / rowSums(Y)
      #mean_expr <- mean(  Yprop, na.rm = T)
      #prop.matrix[gene_i,]= c(geneID, mean_expr)}
    #prop.matrix$V2=lapply(prop.matrix$V2, as.numeric)
    #quantile(unlist(prop.matrix$V2), probs = 0.5, na.rm=TRUE)
    
    ##then actually run the GLM##
    cp.matrix <- matrix(data = NA, nrow = ncol(readcount.cp.RSL), ncol = 10)
    cp.matrix <- as.data.frame(cp.matrix)
    readcount.cp.RSL <- readcount.cp.RSL[,-1]
    geneIDs <- names(readcount.cp.RSL)
    readcount.cp.RSL=sapply(readcount.cp.RSL,as.numeric)
    
    for(gene_i in 1:ncol(readcount.cp.RSL)){
      readcount.cp <- as.data.frame(readcount.cp.RSL)
      geneID <- geneIDs[gene_i]
      y <- readcount.cp.RSL[,gene_i]
      Totalreads <- rowSums(readcount.cp.RSL)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      if(mean_expr > 2.293504e-05  ){ #50% quantile mean read count threshold
        model <- glm(Y ~ metadata.cp.RSL$Treatment * metadata.cp.RSL$DPI, family = "quasibinomial")
        results <- anova(model, test = "Chisq")
        cp.matrix[gene_i,] <- c( geneID, mean_expr, results$`Pr(>Chi)`,model$coefficients)
      }  else{
        cp.matrix[gene_i,] <- c( geneID, mean_expr, NA, NA, NA, NA,NA, NA, NA,  NA, NA, NA, NA, NA)
      }
    }
    
    ##rename matrix columns
    names(cp.matrix) <- c("gene", "meanreadprop", "NA", "Treat_P", "DPI_P", "Treat_DPI_P", "Int_C",
                          "Treat_C","DPI_C","Treat_DPI_C")
    ##format and output matrices
    head(cp.matrix)
    cp.matrix <- cp.matrix[order(cp.matrix$Treat_P, decreasing = F),]
    cp.matrix[1:20,]
    cp.matrix <- cp.matrix[,-3]
    cp.matrix <- cp.matrix[,-6]
    dim(na.omit(cp.matrix))
    cp.matrix.final=na.omit(cp.matrix)
    dim(cp.matrix.final)
    
    write.csv(cp.matrix.final, "RSL_CP_GLM_matrix_results.csv", row.names=FALSE)
    write.csv(cp.matrix, "RSL_CP_GLM_matrix_results_everything.csv", row.names=FALSE)}


##Finally we have GOS
##starting with Alum

  ### Pull relevant data and relevel
  readcount.al <- readcount[metadata$Treatment %in% c("Alum", "PBS"),]
  metadata.al <- metadata[metadata$Treatment %in% c("Alum", "PBS"),]
  readcount.al.GOS <- readcount.al[metadata.al$Population %in% c("GOS"),]
  metadata.al.GOS <- metadata.al[metadata.al$Population %in% c("GOS"),]
  metadata.al.GOS$Treatment<- relevel(factor(metadata.al.GOS$Treatment), ref = "PBS")
  readcount.al.GOS<- as.data.frame(readcount.al.GOS)
  
  ####### SKIP THIS IF RERUNNING and instead reload output
  ###### GLM of alum effects to identify instances with a main effect of treatment, time, treatment*time, etc
  {##start by finding the actual quartile to use##
    prop.matrix <- matrix(data = NA, nrow = ncol(readcount.al), ncol = 2)
    prop.matrix <- as.data.frame(prop.matrix)
    readcount.al.GOS <- readcount.al.GOS[,-1]
    geneIDs <- names(readcount.al.GOS)
    readcount.al.GOS=sapply(readcount.al.GOS,as.numeric)
    
    for(gene_i in 1:ncol(readcount.al.GOS)){
    readcount.al.GOS <- as.data.frame(readcount.al.GOS)
    geneID <- geneIDs[gene_i]
    y <- readcount.al.GOS[,gene_i]
    Totalreads <- rowSums(readcount.al.GOS)
    Y <- cbind(focalgene = y, otherreads = Totalreads - y)
    Y <- as.matrix(Y)
    Yprop <- Y[,1] / rowSums(Y)
    mean_expr <- mean(  Yprop, na.rm = T)
    prop.matrix[gene_i,]= c(geneID, mean_expr)}
    prop.matrix$V2=lapply(prop.matrix$V2, as.numeric)
    quantile(unlist(prop.matrix$V2), probs = 0.5, na.rm = TRUE)
    
    ##GLM actually starts here##
    alum.matrix <- matrix(data = NA, nrow = ncol(readcount.al.GOS), ncol = 10)
    alum.matrix <- as.data.frame(alum.matrix)
    readcount.al.GOS <- readcount.al.GOS[,-1]
    geneIDs <- names(readcount.al.GOS)
    readcount.al.GOS=sapply(readcount.al.GOS,as.numeric)
    
    for(gene_i in 1:ncol(readcount.al.GOS)){
      readcount.al <- as.data.frame(readcount.al.GOS)
      geneID <- geneIDs[gene_i]
      y <- readcount.al.GOS[,gene_i]
      Totalreads <- rowSums(readcount.al.GOS)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      if(mean_expr > 2.573403e-05    ){ #50% quantile mean read count threshold
        model <- glm(Y ~ metadata.al.GOS$Treatment * metadata.al.GOS$DPI, family = "quasibinomial")
        results <- anova(model, test = "Chisq")
        alum.matrix[gene_i,] <- c( geneID, mean_expr, results$`Pr(>Chi)`,model$coefficients)
      }  else{
        alum.matrix[gene_i,] <- c( geneID, mean_expr, NA, NA, NA, NA,NA, NA, NA,  NA, NA, NA, NA, NA)
      }
    }
    
    ##fix matrix names
    names(alum.matrix) <- c("gene", "meanreadprop", "NA", "Treat_P", "DPI_P", "Treat_DPI_P", "Int_C",
                            "Treat_C","DPI_C","Treat_DPI_C")
    ##format and output matrix
    head(alum.matrix)
    alum.matrix <- alum.matrix[order(alum.matrix$Treat_P, decreasing = F),]
    alum.matrix[1:20,]
    alum.matrix <- alum.matrix[,-3]
    alum.matrix <- alum.matrix[,-6]
    dim(na.omit(alum.matrix))
    alum.matrix.final=na.omit(alum.matrix)
    #dim(alum.matrix)
    
    write.csv(alum.matrix.final, "GOS_Alum_GLM_matrix_results.csv", row.names=FALSE)
    write.csv(alum.matrix, "GOS_Alum_GLM_matrix_results_everything.csv", row.names=FALSE)
  }

####### Then GOS CP
  ### Select relevant data and relevel
  readcount.cp <- readcount[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  metadata.cp <- metadata[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  readcount.cp.GOS <- readcount.cp[metadata.cp$Population %in% c("GOS"),]
  metadata.cp.GOS <- metadata.cp[metadata.cp$Population %in% c("GOS"),]
  metadata.cp.GOS$Treatment<- relevel(factor(metadata.cp.GOS$Treatment), ref = "PBS")
  readcount.cp.GOS<- as.data.frame(readcount.cp.GOS)
  readcount.cp.GOS <- readcount.cp.GOS[,-1]
  geneIDs <- names(readcount.cp.GOS)
  readcount.cp.GOS=sapply(readcount.cp.GOS,as.numeric)
  
  ####### SKIP THIS IF RERUNNING and instead reload output
  ###### GLM of alum effects to identify instances with a main effect of treatment, time, treatment*time, etc
  
  {##start by finding the actual quartile to use##
    prop.matrix <- matrix(data = NA, nrow = ncol(readcount.cp.GOS), ncol = 2)
    prop.matrix <- as.data.frame(prop.matrix)
    
    for(gene_i in 1:ncol(readcount.cp.GOS)){
      readcount.cp.GOS <- as.data.frame(readcount.cp.GOS)
      geneID <- geneIDs[gene_i]
      y <- readcount.cp.GOS[,gene_i]
      Totalreads <- rowSums(readcount.cp.GOS)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      prop.matrix[gene_i,]= c(geneID, mean_expr)}
    prop.matrix$V2=lapply(prop.matrix$V2, as.numeric)
    quantile(unlist(prop.matrix$V2), probs = 0.5, na.rm=TRUE)
    
    ##run the actual GLM##
    cp.matrix <- matrix(data = NA, nrow = ncol(readcount.cp.GOS), ncol = 10)
    cp.matrix <- as.data.frame(cp.matrix)
    readcount.cp.GOS <- readcount.cp.GOS[,-1]
    geneIDs <- names(readcount.cp.GOS)
    readcount.cp.GOS=sapply(readcount.cp.GOS,as.numeric)
    
    for(gene_i in 1:ncol(readcount.cp.GOS)){
      readcount.cp <- as.data.frame(readcount.cp.GOS)
      geneID <- geneIDs[gene_i]
      y <- readcount.cp.GOS[,gene_i]
      Totalreads <- rowSums(readcount.cp.GOS)
      Y <- cbind(focalgene = y, otherreads = Totalreads - y)
      Y <- as.matrix(Y)
      Yprop <- Y[,1] / rowSums(Y)
      mean_expr <- mean(  Yprop, na.rm = T)
      if(mean_expr > 2.603938e-05  ){ #50% quantile mean read count threshold
        model <- glm(Y ~ metadata.cp.GOS$Treatment * metadata.cp.GOS$DPI, family = "quasibinomial")
        results <- anova(model, test = "Chisq")
        cp.matrix[gene_i,] <- c( geneID, mean_expr, results$`Pr(>Chi)`,model$coefficients)
      }  else{
        cp.matrix[gene_i,] <- c( geneID, mean_expr, NA, NA, NA, NA,NA, NA, NA,  NA, NA, NA, NA, NA)
      }
    }
    
    ##rename columns
    names(cp.matrix) <- c("gene", "meanreadprop", "NA", "Treat_P", "DPI_P", "Treat_DPI_P", "Int_C",
                          "Treat_C","DPI_C","Treat_DPI_C")
    ##format and output matrix
    head(cp.matrix)
    cp.matrix <- cp.matrix[order(cp.matrix$Treat_P, decreasing = F),]
    cp.matrix[1:20,]
    cp.matrix <- cp.matrix[,-3]
    cp.matrix <- cp.matrix[,-6]
    dim(na.omit(cp.matrix))
    cp.matrix.final=na.omit(cp.matrix)
    dim(cp.matrix.final)
    write.csv(cp.matrix.final, "GOS_CP_GLM_matrix_results.csv", row.names=FALSE)
    write.csv(cp.matrix, "GOS_CP_GLM_matrix_results_everything.csv", row.names=FALSE)
    }


##So now that we've done all that we can use the results to check for correlations between groups of interest in LFC of significant genes##
##start by importing and filtering our data from each GLM##
  ##GOS Alum##
  GOS_al=read.csv("GOS_Alum_GLM_matrix_results_everything.csv")
  GOS_al_treat=GOS_al[,c(1,3,6)]    
  GOS_al_treattime=GOS_al[,c(1,5,8)]
  GOS_al_treat_sig=GOS_al_treat[GOS_al_treat$Treat_P < .01, ] ##selecting significant treat genes
  GOS_al_treattime_sig=GOS_al_treattime[GOS_al_treattime$Treat_DPI_P < .01, ] ##selecting significant treat*time genes##
  ##format each one
  GOS_al_treat_sig=as.data.frame(GOS_al_treat_sig[,c(1)])
  GOS_al_treat_sig=as.data.frame(GOS_al_treat_sig[complete.cases(GOS_al_treat_sig),])
  names(GOS_al_treat_sig)[1] <- "gene"
  GOS_al_treattime_sig=as.data.frame(GOS_al_treattime_sig[,c(1)])
  GOS_al_treattime_sig=as.data.frame(GOS_al_treattime_sig[complete.cases(GOS_al_treattime_sig),])
  names(GOS_al_treattime_sig)[1] <- "gene"
  
  ##GOS CP##
  GOS_cp=read.csv("GOS_CP_GLM_matrix_results_everything.csv")
  GOS_cp_treat=GOS_cp[,c(1,3,6)]    
  GOS_cp_treattime=GOS_cp[,c(1,5,8)]
  GOS_cp_treat_sig=GOS_cp_treat[GOS_cp_treat$Treat_P < .01, ] ##selecting significant treat genes
  GOS_cp_treattime_sig=GOS_cp_treattime[GOS_cp_treattime$Treat_DPI_P < .01, ] ##selecting significant treat*time genes
  ##format each one
  GOS_cp_treat_sig=as.data.frame(GOS_cp_treat_sig[,c(1)])
  GOS_cp_treat_sig=as.data.frame(GOS_cp_treat_sig[complete.cases(GOS_cp_treat_sig),])
  names(GOS_cp_treat_sig)[1] <- "gene"
  GOS_cp_treattime_sig=as.data.frame(GOS_cp_treattime_sig[,c(1)])
  GOS_cp_treattime_sig=as.data.frame(GOS_cp_treattime_sig[complete.cases(GOS_cp_treattime_sig),])
  names(GOS_cp_treattime_sig)[1] <- "gene"
  
  ##RSL Alum##
  RSL_al=read.csv("RSL_Alum_GLM_matrix_results_everything.csv")
  RSL_al_treat=RSL_al[,c(1,3,6)]    
  RSL_al_treattime=RSL_al[,c(1,5,8)]
  RSL_al_treat_sig=RSL_al_treat[RSL_al_treat$Treat_P < .01, ] ##selecting significant treat genes
  RSL_al_treattime_sig=RSL_al_treattime[RSL_al_treattime$Treat_DPI_P < .01, ] ##selecting significant treat*time genes genes
  ##format each one##
  RSL_al_treat_sig=as.data.frame(RSL_al_treat_sig[,c(1)])
  RSL_al_treat_sig=as.data.frame(RSL_al_treat_sig[complete.cases(RSL_al_treat_sig),])
  names(RSL_al_treat_sig)[1] <- "gene"
  RSL_al_treattime_sig=as.data.frame(RSL_al_treattime_sig[,c(1)])
  RSL_al_treattime_sig=as.data.frame(RSL_al_treattime_sig[complete.cases(RSL_al_treattime_sig),])
  names(RSL_al_treattime_sig)[1] <- "gene"
  
  ##RSL CP##
  RSL_cp=read.csv("RSL_CP_GLM_matrix_results_everything.csv")
  RSL_cp_treat=RSL_cp[,c(1,3,6)]    
  RSL_cp_treattime=RSL_cp[,c(1,5,8)]
  RSL_cp_treat_sig=RSL_cp_treat[RSL_cp_treat$Treat_P < .01, ] ##selecting significant treat genes
  RSL_cp_treattime_sig=RSL_cp_treattime[RSL_cp_treattime$Treat_DPI_P < .01, ] ##selecting significant treat*time genes
  ##format each one##
  RSL_cp_treat_sig=as.data.frame(RSL_cp_treat_sig[,c(1)])
  RSL_cp_treat_sig=as.data.frame(RSL_cp_treat_sig[complete.cases(RSL_cp_treat_sig),])
  names(RSL_cp_treat_sig)[1] <- "gene"
  RSL_cp_treattime_sig=as.data.frame(RSL_cp_treattime_sig[,c(1)])
  RSL_cp_treattime_sig=as.data.frame(RSL_cp_treattime_sig[complete.cases(RSL_cp_treattime_sig),])
  names(RSL_cp_treattime_sig)[1] <- "gene"
  
  ##SAY ALUM##
  SAY_al=read.csv("SAY_Alum_GLM_matrix_results_everything.csv")
  SAY_al_treat=SAY_al[,c(1,3,6)]    
  SAY_al_treattime=SAY_al[,c(1,5,8)]
  SAY_al_treat_sig=SAY_al_treat[SAY_al_treat$Treat_P < .01, ] ##selecting significant treat genes
  SAY_al_treattime_sig=SAY_al_treattime[SAY_al_treattime$Treat_DPI_P < .01, ] ##selecting significant treat*time genes
  ##format each one##
  SAY_al_treat_sig=as.data.frame(SAY_al_treat_sig[,c(1)])
  SAY_al_treat_sig=as.data.frame(SAY_al_treat_sig[complete.cases(SAY_al_treat_sig),])
  names(SAY_al_treat_sig)[1] <- "gene"
  SAY_al_treattime_sig=as.data.frame(SAY_al_treattime_sig[,c(1)])
  SAY_al_treattime_sig=as.data.frame(SAY_al_treattime_sig[complete.cases(SAY_al_treattime_sig),])
  names(SAY_al_treattime_sig)[1] <- "gene"
  
  ##SAY CP##
  SAY_cp=read.csv("SAY_CP_GLM_matrix_results_everything.csv")
  SAY_cp_treat=SAY_cp[,c(1,3,6)]    
  SAY_cp_treattime=SAY_cp[,c(1,5,8)]
  SAY_cp_treat_sig=SAY_cp_treat[SAY_cp_treat$Treat_P < .01, ] ##selecting significant treat genes
  SAY_cp_treattime_sig=SAY_cp_treattime[SAY_cp_treattime$Treat_DPI_P < .01, ] ##selecting significant treat*time genes
  ##format each one##
  SAY_cp_treat_sig=as.data.frame(SAY_cp_treat_sig[,c(1)])
  SAY_cp_treat_sig=as.data.frame(SAY_cp_treat_sig[complete.cases(SAY_cp_treat_sig),])
  names(SAY_cp_treat_sig)[1] <- "gene"
  SAY_cp_treattime_sig=as.data.frame(SAY_cp_treattime_sig[,c(1)])
  SAY_cp_treattime_sig=as.data.frame(SAY_cp_treattime_sig[complete.cases(SAY_cp_treattime_sig),])
  names(SAY_cp_treattime_sig)[1] <- "gene"

##merge our lists of sig genes together for each pairwise combo##
  ##across pop alum##
  GvR_al_treat=rbind(GOS_al_treat_sig,RSL_al_treat_sig)
  GvS_al_treat=rbind(GOS_al_treat_sig,SAY_al_treat_sig)
  RvS_al_treat=rbind(SAY_al_treat_sig,RSL_al_treat_sig)
  GvR_al_treattime=rbind(GOS_al_treattime_sig,RSL_al_treattime_sig)
  GvS_al_treattime=rbind(GOS_al_treattime_sig,SAY_al_treattime_sig)
  RvS_al_treattime=rbind(SAY_al_treattime_sig,RSL_al_treattime_sig)

  ##across pop cp##
  GvR_cp_treat=rbind(GOS_cp_treat_sig,RSL_cp_treat_sig)
  GvS_cp_treat=rbind(GOS_cp_treat_sig,SAY_cp_treat_sig)
  RvS_cp_treat=rbind(SAY_cp_treat_sig,RSL_cp_treat_sig)
  GvR_cp_treattime=rbind(GOS_cp_treattime_sig,RSL_cp_treattime_sig)
  GvS_cp_treattime=rbind(GOS_cp_treattime_sig,SAY_cp_treattime_sig)
  RvS_cp_treattime=rbind(SAY_cp_treattime_sig,RSL_cp_treattime_sig)


  ##within pop, across treatments##
  G_treat=rbind(GOS_cp_treat_sig,GOS_al_treat_sig)
  R_treat=rbind(RSL_cp_treat_sig,RSL_al_treat_sig)
  S_treat=rbind(SAY_cp_treat_sig,SAY_al_treat_sig)  
  G_treattime=rbind(GOS_cp_treattime_sig,GOS_al_treattime_sig)
  R_treattime=rbind(RSL_cp_treattime_sig,RSL_al_treattime_sig)
  S_treattime=rbind(SAY_cp_treattime_sig,SAY_al_treattime_sig) 

##now merge it with the LFC data, then run cors##
##done individually for each comparison##
##start with main effects of treatments##
  
  ##GvR alum treat main##
    ##process##
    GvR_al_treat_cor=merge(GvR_al_treat, GOS_al_treat)
    GvR_al_treat_cor=GvR_al_treat_cor[,c(1,3)]
    names(GvR_al_treat_cor)[2] <- "GOS.lfc"
    GvR_al_treat_cor=merge(GvR_al_treat_cor, RSL_al_treat)
    GvR_al_treat_cor=GvR_al_treat_cor[,c(1:2,4)]
    names(GvR_al_treat_cor)[3] <- "RSL.lfc"
    GvR_al_treat_cor=GvR_al_treat_cor[!duplicated(GvR_al_treat_cor), ]
    GvR_al_treat_cor<- na.omit(GvR_al_treat_cor) ##812
    ##stats and graphs##
    cor.test(GvR_al_treat_cor$GOS.lfc,GvR_al_treat_cor$RSL.lfc, method="pearson")##sig, p<0.001, r=0.446
    
      ##graph if you want it##
      GvR_al_treat_cor$effect=GvR_al_treat_cor$GOS.lfc*GvR_al_treat_cor$RSL.lfc
      GvR_al_treat_cor=GvR_al_treat_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
      GvR_al_treat_cor=GvR_al_treat_cor[,c(1:3,5)]
      ggplot(GvR_al_treat_cor, aes(y=GOS.lfc, x=RSL.lfc,color=color)) + 
        geom_point() + 
        geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
        theme_classic() +
        scale_color_manual(values=c("#004f7a", "#cc3d24")) +
        theme(legend.background = element_rect(fill="white",
                                               size=0.5, linetype="solid", 
                                               colour ="black")) +
        xlab("RSL Coefficient") + ylab("GOS Coefficient") +
        theme(legend.position="none") 
    
    ##GvR cp main
    GvR_cp_treat_cor=merge(GvR_cp_treat, GOS_cp_treat)
    GvR_cp_treat_cor=GvR_cp_treat_cor[,c(1,3)]
    names(GvR_cp_treat_cor)[2] <- "GOS.lfc"
    GvR_cp_treat_cor=merge(GvR_cp_treat_cor, RSL_cp_treat)
    GvR_cp_treat_cor=GvR_cp_treat_cor[,c(1:2,4)]
    names(GvR_cp_treat_cor)[3] <- "RSL.lfc"
    GvR_cp_treat_cor=GvR_cp_treat_cor[!duplicated(GvR_cp_treat_cor), ]
    GvR_cp_treat_cor<- na.omit(GvR_cp_treat_cor) ##134
    ##stats and graphs##
    cor.test(GvR_cp_treat_cor$GOS.lfc,GvR_cp_treat_cor$RSL.lfc, method="pearson")##sig, p=0.01183, r=-0.217
    
      ##graph if you want it##
      GvR_cp_treat_cor$effect=GvR_cp_treat_cor$GOS.lfc*GvR_cp_treat_cor$RSL.lfc
      GvR_cp_treat_cor=GvR_cp_treat_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
      GvR_cp_treat_cor=GvR_cp_treat_cor[,c(1:3,5)]
      ggplot(GvR_cp_treat_cor, aes(y=GOS.lfc, x=RSL.lfc,color=color)) + 
        geom_point() + 
        geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
        theme_classic() +
        scale_color_manual(values=c("#004f7a", "#cc3d24")) +
        theme(legend.background = element_rect(fill="white",
                                               size=0.5, linetype="solid", 
                                               colour ="black")) +
        xlab("RSL Coefficient") + ylab("GOS Coefficient") +
        theme(legend.position="none") 
  
  
  ##GvS alum treat main##
    ##process##
    GvS_al_treat_cor=merge(GvS_al_treat, GOS_al_treat)
    GvS_al_treat_cor=GvS_al_treat_cor[,c(1,3)]
    names(GvS_al_treat_cor)[2] <- "GOS.lfc"
    GvS_al_treat_cor=merge(GvS_al_treat_cor, SAY_al_treat)
    GvS_al_treat_cor=GvS_al_treat_cor[,c(1:2,4)]
    names(GvS_al_treat_cor)[3] <- "SAY.lfc"
    GvS_al_treat_cor=GvS_al_treat_cor[!duplicated(GvS_al_treat_cor), ]
    GvS_al_treat_cor<- na.omit(GvS_al_treat_cor) ##882
    ##stats and graphs##
    cor.test(GvS_al_treat_cor$GOS.lfc,GvS_al_treat_cor$SAY.lfc, method="pearson")##sig, p<0.001, r=0.418
    
    ##graph if you want##
      GvS_al_treat_cor$effect=GvS_al_treat_cor$GOS.lfc*GvS_al_treat_cor$SAY.lfc
      GvS_al_treat_cor=GvS_al_treat_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
      GvS_al_treat_cor=GvS_al_treat_cor[,c(1:3,5)]
      ggplot(GvS_al_treat_cor, aes(y=GOS.lfc, x=SAY.lfc,color=color)) + 
        geom_point() + 
        geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
        theme_classic() +
        scale_color_manual(values=c("#004f7a", "#cc3d24")) +
        theme(legend.background = element_rect(fill="white",
                                               size=0.5, linetype="solid", 
                                               colour ="black")) +
        xlab("SAY Coefficient") + ylab("GOS Coefficient") +
        theme(legend.position="none") 
    
    ##GvS cp main
    GvS_cp_treat_cor=merge(GvS_cp_treat, GOS_cp_treat)
    GvS_cp_treat_cor=GvS_cp_treat_cor[,c(1,3)]
    names(GvS_cp_treat_cor)[2] <- "GOS.lfc"
    GvS_cp_treat_cor=merge(GvS_cp_treat_cor, SAY_cp_treat)
    GvS_cp_treat_cor=GvS_cp_treat_cor[,c(1:2,4)]
    names(GvS_cp_treat_cor)[3] <- "SAY.lfc"
    GvS_cp_treat_cor=GvS_cp_treat_cor[!duplicated(GvS_cp_treat_cor), ]
    GvS_cp_treat_cor<- na.omit(GvS_cp_treat_cor) ##153
    ##stats and graphs##
    cor.test(GvS_cp_treat_cor$GOS.lfc,GvS_cp_treat_cor$SAY.lfc, method="pearson")##ns (p=0.277, r=-0.089)

    ##graphing code, but it's not sig##
      #GvS_cp_treat_cor$effect=GvS_cp_treat_cor$GOS.lfc*GvS_cp_treat_cor$SAY.lfc
      #GvS_cp_treat_cor=GvS_cp_treat_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
      #GvS_cp_treat_cor=GvS_cp_treat_cor[,c(1:3,5)]
      #ggplot(GvS_cp_treat_cor, aes(y=GOS.lfc, x=SAY.lfc,color=color)) + 
      #geom_point() + 
      #geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
      #theme_classic() +
      #scale_color_manual(values=c("#004f7a", "#cc3d24")) +
      #theme(legend.background = element_rect(fill="white",
      #size=0.5, linetype="solid", 
      #colour ="black")) +
      #xlab("SAY Coefficient") + ylab("GOS Coefficient") +
      #theme(legend.position="none") 
  
  ##RvS alum treat main##
    ##process##
    RvS_al_treat_cor=merge(RvS_al_treat, RSL_al_treat)
    RvS_al_treat_cor=RvS_al_treat_cor[,c(1,3)]
    names(RvS_al_treat_cor)[2] <- "RSL.lfc"
    RvS_al_treat_cor=merge(RvS_al_treat_cor, SAY_al_treat)
    RvS_al_treat_cor=RvS_al_treat_cor[,c(1:2,4)]
    names(RvS_al_treat_cor)[3] <- "SAY.lfc"
    RvS_al_treat_cor=RvS_al_treat_cor[!duplicated(RvS_al_treat_cor), ]
    RvS_al_treat_cor<- na.omit(RvS_al_treat_cor) ##393
    ##stats and graphs##
    cor.test(RvS_al_treat_cor$RSL.lfc,RvS_al_treat_cor$SAY.lfc, method="pearson")##sig, p<0.001, r=0.241
    
    ##graph if you want it##
      RvS_al_treat_cor$effect=RvS_al_treat_cor$RSL.lfc*RvS_al_treat_cor$SAY.lfc
      RvS_al_treat_cor=RvS_al_treat_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
      RvS_al_treat_cor=RvS_al_treat_cor[,c(1:3,5)]
      ggplot(RvS_al_treat_cor, aes(y=RSL.lfc, x=SAY.lfc,color=color)) + 
        geom_point() + 
        geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
        theme_classic() +
        scale_color_manual(values=c("#004f7a", "#cc3d24")) +
        theme(legend.background = element_rect(fill="white",
                                               size=0.5, linetype="solid", 
                                               colour ="black")) +
        xlab("SAY Coefficient") + ylab("RSL Coefficient") +
        theme(legend.position="none") 
    
    ##RvS cp main
    RvS_cp_treat_cor=merge(RvS_cp_treat, RSL_cp_treat)
    RvS_cp_treat_cor=RvS_cp_treat_cor[,c(1,3)]
    names(RvS_cp_treat_cor)[2] <- "RSL.lfc"
    RvS_cp_treat_cor=merge(RvS_cp_treat_cor, SAY_cp_treat)
    RvS_cp_treat_cor=RvS_cp_treat_cor[,c(1:2,4)]
    names(RvS_cp_treat_cor)[3] <- "SAY.lfc"
    RvS_cp_treat_cor=RvS_cp_treat_cor[!duplicated(RvS_cp_treat_cor), ]
    RvS_cp_treat_cor<- na.omit(RvS_cp_treat_cor) ##158
    ##stats and graphs##
    cor.test(RvS_cp_treat_cor$RSL.lfc,RvS_cp_treat_cor$SAY.lfc, method="pearson")##sig, p=0.0317, r=-0.171
    
    ##graph if you want it
      RvS_cp_treat_cor$effect=RvS_cp_treat_cor$RSL.lfc*RvS_cp_treat_cor$SAY.lfc
      RvS_cp_treat_cor=RvS_cp_treat_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
      RvS_cp_treat_cor=RvS_cp_treat_cor[,c(1:3,5)]
      ggplot(RvS_cp_treat_cor, aes(y=RSL.lfc, x=SAY.lfc,color=color)) + 
        geom_point() + 
        geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
        theme_classic() +
        scale_color_manual(values=c("#004f7a", "#cc3d24")) +
        theme(legend.background = element_rect(fill="white",
                                               size=0.5, linetype="solid", 
                                               colour ="black")) +
        xlab("SAY Coefficient") + ylab("RSL Coefficient") +
        theme(legend.position="none") 
    

##Now do the same thing for treat*time##

  ##GvR alum treat int##
    ##process##
    GvR_al_treattime_cor=merge(GvR_al_treattime, GOS_al_treattime)
    GvR_al_treattime_cor=GvR_al_treattime_cor[,c(1,3)]
    names(GvR_al_treattime_cor)[2] <- "GOS.lfc"
    GvR_al_treattime_cor=merge(GvR_al_treattime_cor, RSL_al_treattime)
    GvR_al_treattime_cor=GvR_al_treattime_cor[,c(1:2,4)]
    names(GvR_al_treattime_cor)[3] <- "RSL.lfc"
    GvR_al_treattime_cor=GvR_al_treattime_cor[!duplicated(GvR_al_treattime_cor), ]
    GvR_al_treattime_cor<- na.omit(GvR_al_treattime_cor) ##93
    ##stats and graphs##
    cor.test(GvR_al_treattime_cor$GOS.lfc,GvR_al_treattime_cor$RSL.lfc, method="pearson")##ns (p=0.53, r=0.066)
    
    ##graph code- but it's not sig##
      #GvR_al_treattime_cor$effect=GvR_al_treattime_cor$GOS.lfc*GvR_al_treattime_cor$RSL.lfc
      #GvR_al_treattime_cor=GvR_al_treattime_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
      #GvR_al_treattime_cor=GvR_al_treattime_cor[,c(1:3,5)]
      #ggplot(GvR_al_treattime_cor, aes(y=GOS.lfc, x=RSL.lfc,color=color)) + 
      #geom_point() + 
      #geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
      #theme_classic() +
      #scale_color_manual(values=c("#004f7a", "#cc3d24")) +
      #theme(legend.background = element_rect(fill="white",
      #size=0.5, linetype="solid", 
      #colour ="black")) +
      #xlab("RSL Coefficient") + ylab("GOS Coefficient") +
      #theme(legend.position="none") 
    
    ##GvR cp int
    GvR_cp_treattime_cor=merge(GvR_cp_treattime, GOS_cp_treattime)
    GvR_cp_treattime_cor=GvR_cp_treattime_cor[,c(1,3)]
    names(GvR_cp_treattime_cor)[2] <- "GOS.lfc"
    GvR_cp_treattime_cor=merge(GvR_cp_treattime_cor, RSL_cp_treattime)
    GvR_cp_treattime_cor=GvR_cp_treattime_cor[,c(1:2,4)]
    names(GvR_cp_treattime_cor)[3] <- "RSL.lfc"
    GvR_cp_treattime_cor=GvR_cp_treattime_cor[!duplicated(GvR_cp_treattime_cor), ]
    GvR_cp_treattime_cor<- na.omit(GvR_cp_treattime_cor) ##93
    ##stats and graphs##
    cor.test(GvR_cp_treattime_cor$GOS.lfc,GvR_cp_treattime_cor$RSL.lfc, method="pearson")##ns (0.0592, 0.120)
    
    ##graph code- but it's not sig##
      #GvR_cp_treattime_cor$effect=GvR_cp_treattime_cor$GOS.lfc*GvR_cp_treattime_cor$RSL.lfc
      #GvR_cp_treattime_cor=GvR_cp_treattime_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
      #GvR_cp_treattime_cor=GvR_cp_treattime_cor[,c(1:3,5)]
      #ggplot(GvR_cp_treattime_cor, aes(y=GOS.lfc, x=RSL.lfc,color=color)) + 
      #geom_point() + 
      #geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
      #theme_classic() +
      #scale_color_manual(values=c("#004f7a", "#cc3d24")) +
      #theme(legend.background = element_rect(fill="white",
      #size=0.5, linetype="solid", 
      #colour ="black")) +
      #xlab("RSL Coefficient") + ylab("GOS Coefficient") +
      #theme(legend.position="none") 
  
  
  ##GvS alum treat int##
    ##process##
    GvS_al_treattime_cor=merge(GvS_al_treattime, GOS_al_treattime)
    GvS_al_treattime_cor=GvS_al_treattime_cor[,c(1,3)]
    names(GvS_al_treattime_cor)[2] <- "GOS.lfc"
    GvS_al_treattime_cor=merge(GvS_al_treattime_cor, SAY_al_treattime)
    GvS_al_treattime_cor=GvS_al_treattime_cor[,c(1:2,4)]
    names(GvS_al_treattime_cor)[3] <- "SAY.lfc"
    GvS_al_treattime_cor=GvS_al_treattime_cor[!duplicated(GvS_al_treattime_cor), ]
    GvS_al_treattime_cor<- na.omit(GvS_al_treattime_cor) ##93
    ##stats and graphs##
    cor.test(GvS_al_treattime_cor$GOS.lfc,GvS_al_treattime_cor$SAY.lfc, method="pearson")##sig, p<0.001, r=0.158
    
    ##graph if you want##
      GvS_al_treattime_cor$effect=GvS_al_treattime_cor$GOS.lfc*GvS_al_treattime_cor$SAY.lfc
      GvS_al_treattime_cor=GvS_al_treattime_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
      GvS_al_treattime_cor=GvS_al_treattime_cor[,c(1:3,5)]
      ggplot(GvS_al_treattime_cor, aes(y=GOS.lfc, x=SAY.lfc,color=color)) + 
        geom_point() + 
        geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
        theme_classic() +
        scale_color_manual(values=c("#004f7a", "#cc3d24")) +
        theme(legend.background = element_rect(fill="white",
                                               size=0.5, linetype="solid", 
                                               colour ="black")) +
        xlab("SAY Coefficient") + ylab("GOS Coefficient") +
        theme(legend.position="none") 
    
    ##GvS cp int
    GvS_cp_treattime_cor=merge(GvS_cp_treattime, GOS_cp_treattime)
    GvS_cp_treattime_cor=GvS_cp_treattime_cor[,c(1,3)]
    names(GvS_cp_treattime_cor)[2] <- "GOS.lfc"
    GvS_cp_treattime_cor=merge(GvS_cp_treattime_cor, SAY_cp_treattime)
    GvS_cp_treattime_cor=GvS_cp_treattime_cor[,c(1:2,4)]
    names(GvS_cp_treattime_cor)[3] <- "SAY.lfc"
    GvS_cp_treattime_cor=GvS_cp_treattime_cor[!duplicated(GvS_cp_treattime_cor), ]
    GvS_cp_treattime_cor<- na.omit(GvS_cp_treattime_cor) ##93
    ##stats and graphs##
    cor.test(GvS_cp_treattime_cor$GOS.lfc,GvS_cp_treattime_cor$SAY.lfc, method="pearson")##sig, p<0.001, r=0.275
    
    ##graph if you want##
      GvS_cp_treattime_cor$effect=GvS_cp_treattime_cor$GOS.lfc*GvS_cp_treattime_cor$SAY.lfc
      GvS_cp_treattime_cor=GvS_cp_treattime_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
      GvS_cp_treattime_cor=GvS_cp_treattime_cor[,c(1:3,5)]
      ggplot(GvS_cp_treattime_cor, aes(y=GOS.lfc, x=SAY.lfc,color=color)) + 
        geom_point() + 
        geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
        theme_classic() +
        scale_color_manual(values=c("#004f7a", "#cc3d24")) +
        theme(legend.background = element_rect(fill="white",
                                               size=0.5, linetype="solid", 
                                               colour ="black")) +
        xlab("SAY Coefficient") + ylab("GOS Coefficient") +
        theme(legend.position="none") 
  
  ##RvS alum treat int##
    ##process##
    RvS_al_treattime_cor=merge(RvS_al_treattime, RSL_al_treattime)
    RvS_al_treattime_cor=RvS_al_treattime_cor[,c(1,3)]
    names(RvS_al_treattime_cor)[2] <- "RSL.lfc"
    RvS_al_treattime_cor=merge(RvS_al_treattime_cor, SAY_al_treattime)
    RvS_al_treattime_cor=RvS_al_treattime_cor[,c(1:2,4)]
    names(RvS_al_treattime_cor)[3] <- "SAY.lfc"
    RvS_al_treattime_cor=RvS_al_treattime_cor[!duplicated(RvS_al_treattime_cor), ]
    RvS_al_treattime_cor<- na.omit(RvS_al_treattime_cor) ##251
    ##stats and graphs##
    cor.test(RvS_al_treattime_cor$RSL.lfc,RvS_al_treattime_cor$SAY.lfc, method="pearson")##ns (p=0.052, r=0.098)
    
    ##graph code if you want
      #RvS_al_treattime_cor$effect=RvS_al_treattime_cor$RSL.lfc*RvS_al_treattime_cor$SAY.lfc
      #RvS_al_treattime_cor=RvS_al_treattime_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
      #RvS_al_treattime_cor=RvS_al_treattime_cor[,c(1:3,5)]
      #ggplot(RvS_al_treattime_cor, aes(y=RSL.lfc, x=SAY.lfc,color=color)) + 
      #geom_point() + 
      #geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
      #theme_classic() +
      #scale_color_manual(values=c("#004f7a", "#cc3d24")) +
      #theme(legend.background = element_rect(fill="white",
      #size=0.5, linetype="solid", 
      #colour ="black")) +
      #xlab("SAY Coefficient") + ylab("RSL Coefficient") +
      #theme(legend.position="none") 
    
    ##RvS cp int
    RvS_cp_treattime_cor=merge(RvS_cp_treattime, RSL_cp_treattime)
    RvS_cp_treattime_cor=RvS_cp_treattime_cor[,c(1,3)]
    names(RvS_cp_treattime_cor)[2] <- "RSL.lfc"
    RvS_cp_treattime_cor=merge(RvS_cp_treattime_cor, SAY_cp_treattime)
    RvS_cp_treattime_cor=RvS_cp_treattime_cor[,c(1:2,4)]
    names(RvS_cp_treattime_cor)[3] <- "SAY.lfc"
    RvS_cp_treattime_cor=RvS_cp_treattime_cor[!duplicated(RvS_cp_treattime_cor), ]
    RvS_cp_treattime_cor<- na.omit(RvS_cp_treattime_cor) ##251
    ##stats and graphs##
    cor.test(RvS_cp_treattime_cor$RSL.lfc,RvS_cp_treattime_cor$SAY.lfc, method="pearson")##ns (p=0.072, r=0.152)
    
    ##graph code- but it's not sig##
      #RvS_cp_treattime_cor$effect=RvS_cp_treattime_cor$RSL.lfc*RvS_cp_treattime_cor$SAY.lfc
      #RvS_cp_treattime_cor=RvS_cp_treattime_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
      #RvS_cp_treattime_cor=RvS_cp_treattime_cor[,c(1:3,5)]
      #ggplot(RvS_cp_treattime_cor, aes(y=RSL.lfc, x=SAY.lfc,color=color)) + 
      #geom_point() + 
      #geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
      #theme_classic() +
      #scale_color_manual(values=c("#004f7a", "#cc3d24")) +
      #theme(legend.background = element_rect(fill="white",
      #size=0.5, linetype="solid", 
      #colour ="black")) +
      #xlab("SAY Coefficient") + ylab("RSL Coefficient") +
      #theme(legend.position="none") 
    
    ##we need to pull some numbers for tables, looking specifically at total number of genes tested for corrs##
    dim(GvR_al_treat_cor[complete.cases(GvR_al_treat_cor), ])#812
    dim(GvR_cp_treat_cor[complete.cases(GvR_cp_treat_cor), ])#134
    dim(GvR_al_treattime_cor[complete.cases(GvR_al_treattime_cor), ])#93
    dim(GvR_cp_treattime_cor[complete.cases(GvR_cp_treattime_cor), ])#247
    
    dim(GvS_al_treat_cor[complete.cases(GvS_al_treat_cor), ])#883
    dim(GvS_cp_treat_cor[complete.cases(GvS_cp_treat_cor), ])#153
    dim(GvS_al_treattime_cor[complete.cases(GvS_al_treattime_cor), ])#247
    dim(GvS_cp_treattime_cor[complete.cases(GvS_cp_treattime_cor), ])#292
    
    dim(RvS_al_treat_cor[complete.cases(RvS_al_treat_cor), ])#394
    dim(RvS_cp_treat_cor[complete.cases(RvS_cp_treat_cor), ])#158
    dim(RvS_al_treattime_cor[complete.cases(RvS_al_treattime_cor), ])#251
    dim(RvS_cp_treattime_cor[complete.cases(RvS_cp_treattime_cor), ])#141

##we'll also quickly consider corrs for all genes, not just sig##
    ##initial process##
    GOS_al_treat=GOS_al_treat[,c(1,3)]
    GOS_al_treattime=GOS_al_treattime[,c(1,3)]
    names(GOS_al_treat)[2] <- "GOS.lfc"
    names(GOS_al_treattime)[2] <- "GOS.lfc"
    
    RSL_al_treat=RSL_al_treat[,c(1,3)]
    RSL_al_treattime=RSL_al_treattime[,c(1,3)]
    names(RSL_al_treat)[2] <- "RSL.lfc"
    names(RSL_al_treattime)[2] <- "RSL.lfc"
    
    SAY_al_treat=SAY_al_treat[,c(1,3)]
    SAY_al_treattime=SAY_al_treattime[,c(1,3)]
    names(SAY_al_treat)[2] <- "SAY.lfc"
    names(SAY_al_treattime)[2] <- "SAY.lfc"
    
    GOS_cp_treat=GOS_cp_treat[,c(1,3)]
    GOS_cp_treattime=GOS_cp_treattime[,c(1,3)]
    names(GOS_cp_treat)[2] <- "GOS.lfc"
    names(GOS_cp_treattime)[2] <- "GOS.lfc"
    
    RSL_cp_treat=RSL_cp_treat[,c(1,3)]
    RSL_cp_treattime=RSL_cp_treattime[,c(1,3)]
    names(RSL_cp_treat)[2] <- "RSL.lfc"
    names(RSL_cp_treattime)[2] <- "RSL.lfc"
    
    SAY_cp_treat=SAY_cp_treat[,c(1,3)]
    SAY_cp_treattime=SAY_cp_treattime[,c(1,3)]
    names(SAY_cp_treat)[2] <- "SAY.lfc"
    names(SAY_cp_treattime)[2] <- "SAY.lfc"
    
    ##running actual stats##
    
    ##GvR alum treat main##
    GvR_al_treat_cor=merge(GOS_al_treat,RSL_al_treat)
    GvR_al_treat_cor<- na.omit(GvR_al_treat_cor) ##11770
    cor.test(GvR_al_treat_cor$GOS.lfc,GvR_al_treat_cor$RSL.lfc, method="pearson")##sig, p<0.001, r=0.215
    
    ##GvS alum main
    GvS_al_treat_cor=merge(GOS_al_treat,SAY_al_treat)
    GvS_al_treat_cor<- na.omit(GvS_al_treat_cor) ##12459
    cor.test(GvS_al_treat_cor$GOS.lfc,GvS_al_treat_cor$SAY.lfc, method="pearson")##sig, p<0.001, r=0.282
    
    ##RvS alum main
    RvS_al_treat_cor=merge(RSL_al_treat,SAY_al_treat)
    RvS_al_treat_cor<- na.omit(RvS_al_treat_cor) ##11919
    cor.test(RvS_al_treat_cor$RSL.lfc,RvS_al_treat_cor$SAY.lfc, method="pearson")##sig, p<0.001, r=0.158
    
    ##GvR cp  main##
    GvR_cp_treat_cor=merge(GOS_cp_treat,RSL_cp_treat)
    GvR_cp_treat_cor<- na.omit(GvR_cp_treat_cor) ##11817
    cor.test(GvR_cp_treat_cor$GOS.lfc,GvR_cp_treat_cor$RSL.lfc, method="pearson")##sig, p<0.001, r=-0.112
    
    ##GvS cp main
    GvS_cp_treat_cor=merge(GOS_cp_treat,SAY_cp_treat)
    GvS_cp_treat_cor<- na.omit(GvS_cp_treat_cor) ##12480
    cor.test(GvS_cp_treat_cor$GOS.lfc,GvS_cp_treat_cor$SAY.lfc, method="pearson")##sig, p<0.001, r=-0.044
    
    ##RvS cp main
    RvS_cp_treat_cor=merge(RSL_cp_treat,SAY_cp_treat)
    RvS_cp_treat_cor<- na.omit(RvS_cp_treat_cor) ##11918
    cor.test(RvS_cp_treat_cor$RSL.lfc,RvS_cp_treat_cor$SAY.lfc, method="pearson")##sig, p<0.001, r=0.046
    
    ##GvR alum*time ##
    GvR_al_treattime_cor=merge(GOS_al_treattime,RSL_al_treattime)
    GvR_al_treattime_cor<- na.omit(GvR_al_treattime_cor) ##11770
    cor.test(GvR_al_treattime_cor$GOS.lfc,GvR_al_treattime_cor$RSL.lfc, method="pearson")##sig, p<0.001, r=0.059
    
    ##GvS alum*time
    GvS_al_treattime_cor=merge(GOS_al_treattime,SAY_al_treattime)
    GvS_al_treattime_cor<- na.omit(GvS_al_treattime_cor) ##12459
    cor.test(GvS_al_treattime_cor$GOS.lfc,GvS_al_treattime_cor$SAY.lfc, method="pearson")##sig, p<0.001, r=0.120
    
    ##RvS alum*time
    RvS_al_treattime_cor=merge(RSL_al_treattime,SAY_al_treattime)
    RvS_al_treattime_cor<- na.omit(RvS_al_treattime_cor) ##11919
    cor.test(RvS_al_treattime_cor$RSL.lfc,RvS_al_treattime_cor$SAY.lfc, method="pearson")##sig, p<0.001, r=0.185
    
    ##GvR cp*time
    GvR_cp_treattime_cor=merge(GOS_cp_treattime,RSL_cp_treattime)
    GvR_cp_treattime_cor<- na.omit(GvR_cp_treattime_cor) ##11817
    cor.test(GvR_cp_treattime_cor$GOS.lfc,GvR_cp_treattime_cor$RSL.lfc, method="pearson")##sig, p<0.001, r=0.172
    
    ##GvS cp*time
    GvS_cp_treattime_cor=merge(GOS_cp_treattime,SAY_cp_treattime)
    GvS_cp_treattime_cor<- na.omit(GvS_cp_treattime_cor) ##12480
    cor.test(GvS_cp_treattime_cor$GOS.lfc,GvS_cp_treattime_cor$SAY.lfc, method="pearson")##sig, p<0.001, r=0.274
    
    
    ##RvS cp*time
    RvS_cp_treattime_cor=merge(RSL_cp_treattime,SAY_cp_treattime)
    RvS_cp_treattime_cor<- na.omit(RvS_cp_treattime_cor) ##11918
    cor.test(RvS_cp_treattime_cor$RSL.lfc,RvS_cp_treattime_cor$SAY.lfc, method="pearson")##sig, p<0.001, r=0.199
    

##finally we're going to get numbers of shared sig genes and stats on statistical overrepresentation##
##we'll reprocess data fiels for this##
  GOS.al.matrix=read.csv("GOS_Alum_GLM_matrix_results.csv")
  GOS.cp.matrix=read.csv("GOS_CP_GLM_matrix_results.csv")
  RSL.al.matrix=read.csv("RSL_Alum_GLM_matrix_results.csv")
  RSL.cp.matrix=read.csv("RSL_CP_GLM_matrix_results.csv")
  SAY.al.matrix=read.csv("SAY_Alum_GLM_matrix_results.csv")
  SAY.cp.matrix=read.csv("SAY_CP_GLM_matrix_results.csv")
  GOS.al.matrix<- na.omit(GOS.al.matrix)
  GOS.cp.matrix<- na.omit(GOS.cp.matrix)
  RSL.al.matrix<- na.omit(RSL.al.matrix)
  RSL.cp.matrix<- na.omit(RSL.cp.matrix)
  SAY.al.matrix<- na.omit(SAY.al.matrix)
  SAY.cp.matrix<- na.omit(SAY.cp.matrix)
  alpha <- 0.01
  
  ##look at how many genes are shared in each pairwise for our stats##
  GvR_alum=merge(GOS.al.matrix, RSL.al.matrix, by="gene")
  GvR_alum=na.omit(GvR_alum)
  dim(GvR_alum) ##11770
  GvS_alum=merge(GOS.al.matrix, SAY.al.matrix, by="gene")
  GvS_alum=na.omit(GvS_alum)
  dim(GvS_alum) ##12459
  RvS_alum=merge(RSL.al.matrix, SAY.al.matrix, by="gene")
  RvS_alum=na.omit(RvS_alum)
  dim(RvS_alum) ##11919
  GvR_cp=merge(GOS.cp.matrix, RSL.cp.matrix, by="gene")
  GvR_cp=na.omit(GvR_cp)
  dim(GvR_cp) ##11817
  GvS_cp=merge(GOS.cp.matrix, SAY.cp.matrix, by="gene")
  GvS_cp=na.omit(GvS_cp)
  dim(GvS_cp) ##12480
  RvS_cp=merge(RSL.cp.matrix, SAY.cp.matrix, by="gene")
  RvS_cp=na.omit(RvS_cp)
  dim(RvS_cp) ##11918

  ##merge and format data frames for downstream work##
  merged.matrix.al=merge(GOS.al.matrix,RSL.al.matrix, by="gene",all=T)
  merged.matrix.al=merged.matrix.al[,c(1,3,5,10,12)]
  names(merged.matrix.al)[2] <- "GOS.al.P"
  names(merged.matrix.al)[3] <- "GOS.aldpi.P"
  names(merged.matrix.al)[4] <- "RSL.al.P"
  names(merged.matrix.al)[5] <- "RSL.aldpi.P"
  merged.matrix.al2=merge(merged.matrix.al,SAY.al.matrix, by="gene",all=T)
  merged.matrix.al2=merged.matrix.al2[,c(1:5,7,9)]
  names(merged.matrix.al2)[6] <- "SAY.al.P"
  names(merged.matrix.al2)[7] <- "SAY.aldpi.P"
  merged.matrix.al2 <- na.omit(merged.matrix.al2)
  
  merged.matrix.cp=merge(GOS.cp.matrix,RSL.cp.matrix, by="gene",all=T)
  merged.matrix.cp=merged.matrix.cp[,c(1,3,5,10,12)]
  names(merged.matrix.cp)[2] <- "GOS.cp.P"
  names(merged.matrix.cp)[3] <- "GOS.cpdpi.P"
  names(merged.matrix.cp)[4] <- "RSL.cp.P"
  names(merged.matrix.cp)[5] <- "RSL.cpdpi.P"
  merged.matrix.cp2=merge(merged.matrix.cp,SAY.cp.matrix, by="gene",all=T)
  merged.matrix.cp2=merged.matrix.cp2[,c(1:5,7,9)]
  names(merged.matrix.cp2)[6] <- "SAY.cp.P"
  names(merged.matrix.cp2)[7] <- "SAY.cpdpi.P"
  merged.matrix.cp2 <- na.omit(merged.matrix.cp2)
  
  # Look at overlap between all pairwise comparisons for alum main effects
  ##GOS vs RSL
  shared.trt <- merged.matrix.al2$GOS.al.P < 0.01 & merged.matrix.al2$RSL.al.P < 0.01
  sum(shared.trt)# 35
  #write out genes for later analysis##
  shared.goi <- merged.matrix.al2$gene[shared.trt]
  write.csv(shared.goi, "Shared_pop_GvR_al_sig.csv", row.names=FALSE)
  ##check proportions##
  GOS.prop <- 710/13035
  RSL.prop <- 207/13034
  random_expect <- GOS.prop * RSL.prop
  random_expect*11770  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) # p<0.001
  
  ##GOS vs. SAY
  shared.trt <- merged.matrix.al2$GOS.al.P < 0.01 & merged.matrix.al2$SAY.al.P < 0.01
  sum(shared.trt)# 44
  #write out genes for later analysis##
  shared.goi <- merged.matrix.al2$gene[shared.trt]
  write.csv(shared.goi, "Shared_pop_GvS_al_sig.csv", row.names=FALSE)
  ##check proportions##
  GOS.prop <- 710/13035
  SAY.prop <- 235/13034
  random_expect <- GOS.prop * SAY.prop
  random_expect*12459  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) # p<0.001
  
  ##RSL vs. SAY
  shared.trt <- merged.matrix.al2$RSL.al.P < 0.01 & merged.matrix.al2$SAY.al.P < 0.01
  sum(shared.trt)# 21
  #write out genes for later analysis##
  shared.goi <- merged.matrix.al2$gene[shared.trt]
  write.csv(shared.goi, "Shared_pop_RvS_al_sig.csv", row.names=FALSE)
  ##check proportions##
  RSL.prop <- 207/13034
  SAY.prop <- 235/13034
  random_expect <- RSL.prop * SAY.prop
  random_expect*11919  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) # p<0.001
  
  ##Look at overlap between all pairwise comparisons for cp main effects
  ##GOS vs RSL
  shared.trt <- merged.matrix.cp2$GOS.cp.P < 0.01 & merged.matrix.cp2$RSL.cp.P < 0.01
  sum(shared.trt)# 2
  #write out genes for later ancpysis##
  shared.goi <- merged.matrix.cp2$gene[shared.trt]
  write.csv(shared.goi, "Shared_pop_GvR_cp_sig.csv", row.names=FALSE)
  ##check proportions##
  GOS.prop <- 66/13035
  RSL.prop <- 92/13036
  random_expect <- GOS.prop * RSL.prop
  random_expect*11817  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) #ns (p=0.0915)
  
  ##GOS vs. SAY
  shared.trt <- merged.matrix.cp2$GOS.cp.P < 0.01 & merged.matrix.cp2$SAY.cp.P < 0.01
  sum(shared.trt)# 0
  #write out genes for later ancpysis##
  shared.goi <- merged.matrix.cp2$gene[shared.trt]
  write.csv(shared.goi, "Shared_pop_GvS_cp_sig.csv", row.names=FALSE)
  ##check proportions##
  GOS.prop <- 66/13035
  SAY.prop <- 94/13035
  random_expect <- GOS.prop * SAY.prop
  random_expect*12480  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) # ns (p=1)
  
  ##RSL vs. SAY
  shared.trt <- merged.matrix.cp2$RSL.cp.P < 0.01 & merged.matrix.cp2$SAY.cp.P < 0.01
  sum(shared.trt)# 0
  #write out genes for later ancpysis##
  shared.goi <- merged.matrix.cp2$gene[shared.trt]
  write.csv(shared.goi, "Shared_pop_RvS_cp_sig.csv", row.names=FALSE)
  ##check proportions##
  RSL.prop <- 92/13036
  SAY.prop <- 94/13035
  random_expect <- RSL.prop * SAY.prop
  random_expect*11918  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) # ns (p=0.907)
  
  ##Look at overlap between all pairwise comparisons for alumxtime effects
  ##GOS vs RSL
  shared.trt <- merged.matrix.al2$GOS.aldpi.P < 0.01 & merged.matrix.al2$RSL.aldpi.P < 0.01
  sum(shared.trt)# 0
  #write out genes for later analysis##
  shared.goi <- merged.matrix.al2$gene[shared.trt]
  write.csv(shared.goi, "Shared_pop_GvR_alxdpi_sig.csv", row.names=FALSE)
  ##check proportions##
  GOS.prop <- 41/13035
  RSL.prop <- 66/13034
  random_expect <- GOS.prop * RSL.prop
  random_expect*11770  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) # ns (p=1)
  
  ##GOS vs. SAY
  shared.trt <- merged.matrix.al2$GOS.aldpi.P < 0.01 & merged.matrix.al2$SAY.aldpi.P < 0.01
  sum(shared.trt)# 0
  #write out genes for later analysis##
  shared.goi <- merged.matrix.al2$gene[shared.trt]
  write.csv(shared.goi, "Shared_pop_GvS_alxdpi_sig.csv", row.names=FALSE)
  ##check proportions##
  GOS.prop <- 41/13035
  SAY.prop <- 210/13034
  random_expect <- GOS.prop * SAY.prop
  random_expect*12459  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) # ns (p=.91)
  
  ##RSL vs SAY
  shared.trt <- merged.matrix.al2$RSL.aldpi.P < 0.01 & merged.matrix.al2$SAY.aldpi.P < 0.01
  sum(shared.trt)# 0
  #write out genes for later analysis##
  shared.goi <- merged.matrix.al2$gene[shared.trt]
  write.csv(shared.goi, "Shared_pop_RvS_alxdpi_sig.csv", row.names=FALSE)
  ##check proportions##
  RSL.prop <- 66/13034
  SAY.prop <- 210/13034
  random_expect <- RSL.prop * SAY.prop
  random_expect*11919  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) # ns(p=.648)
  
  #Look at overlap between all pairwise comparisons for cp*time effects
  ##GOS vs RSL
  shared.trt <- merged.matrix.cp2$GOS.cpdpi.P < 0.01 & merged.matrix.cp2$RSL.cpdpi.P < 0.01
  sum(shared.trt)# 0
  #write out genes for later analysis##
  shared.goi <- merged.matrix.cp2$gene[shared.trt]
  write.csv(shared.goi, "Shared_pop_GvR_cpxdpi_sig.csv", row.names=FALSE)
  ##check proportions##
  GOS.prop <- 221/13035
  RSL.prop <- 70/13036
  random_expect <- GOS.prop * RSL.prop
  random_expect*11817  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) #ns (p=0.589)
  
  ##GOS vs SAY
  shared.trt <- merged.matrix.cp2$GOS.cpdpi.P < 0.01 & merged.matrix.cp2$SAY.cpdpi.P < 0.01
  sum(shared.trt)# 1
  #write out genes for later analysis##
  shared.goi <- merged.matrix.cp2$gene[shared.trt]
  write.csv(shared.goi, "Shared_pop_GvS_cpxdpi_sig.csv", row.names=FALSE)
  ##check proportions##
  GOS.prop <- 221/13035
  SAY.prop <- 88/13035
  random_expect <- GOS.prop * SAY.prop
  random_expect*12480  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) # ns (p=1)
  
  ##RSL vs SAY
  shared.trt <- merged.matrix.cp2$RSL.cpdpi.P < 0.01 & merged.matrix.cp2$SAY.cpdpi.P < 0.01
  sum(shared.trt)# 1
  #write out genes for later analysis##
  shared.goi <- merged.matrix.cp2$gene[shared.trt]
  write.csv(shared.goi, "Shared_pop_RvS_cpxdpi_sig.csv", row.names=FALSE)
  ##check proportions##
  RSL.prop <- 71/13036
  SAY.prop <- 88/13035
  random_expect <- RSL.prop * SAY.prop
  random_expect*11918  # number of genes expected to be shared under null
  ##do we deviate from expectation?
  prop.test(x = sum(shared.trt), n = length(shared.trt), p = random_expect) # ns (p=0.910)


##last, out of curiousity, how do responses correlate across treatments within pops?
  ##GOS treat##
    ##process##
    G_treat_cor=merge(G_treat, GOS_al_treat)
    G_treat_cor=G_treat_cor[,c(1,3)]
    names(G_treat_cor)[2] <- "GOS.al"
    G_treat_cor=merge(G_treat_cor, GOS_cp_treat)
    G_treat_cor=G_treat_cor[,c(1:2,4)]
    names(G_treat_cor)[3] <- "GOS.cp"
    G_treat_cor=G_treat_cor[!duplicated(G_treat_cor), ]
    ##stats and graphs##
    cor.test(G_treat_cor$GOS.al,G_treat_cor$GOS.cp, method="pearson")##sig, p<0.001, r=0.679
    
    ##graph if you want##
    G_treat_cor$effect=G_treat_cor$GOS.al*G_treat_cor$GOS.cp
    G_treat_cor=G_treat_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
    G_treat_cor=G_treat_cor[,c(1:3,5)]
    ggplot(G_treat_cor, aes(y=GOS.al, x=GOS.cp,color=color)) + 
      geom_point() + 
      geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
      theme_classic() +
      scale_color_manual(values=c("#004f7a", "#cc3d24")) +
      theme(legend.background = element_rect(fill="white",
                                             size=0.5, linetype="solid", 
                                             colour ="black")) +
      xlab("CP Coefficient") + ylab("Alum Coefficient") +
      theme(legend.position="none") 
  
  ##GOS treatxtime##
    ##process##
    G_treattime_cor=merge(G_treattime, GOS_al_treattime)
    G_treattime_cor=G_treattime_cor[,c(1,3)]
    names(G_treattime_cor)[2] <- "GOS.al"
    G_treattime_cor=merge(G_treattime_cor, GOS_cp_treattime)
    G_treattime_cor=G_treattime_cor[,c(1:2,4)]
    names(G_treattime_cor)[3] <- "GOS.cp"
    G_treattime_cor=G_treattime_cor[!duplicated(G_treattime_cor), ]
    ##stats and graphs##
    cor.test(G_treattime_cor$GOS.al,G_treattime_cor$GOS.cp, method="pearson")##sig, p<0.001, r=0.587
    
    ##graph if you want##
    G_treattime_cor$effect=G_treattime_cor$GOS.al*G_treattime_cor$GOS.cp
    G_treattime_cor=G_treattime_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
    G_treattime_cor=G_treattime_cor[,c(1:3,5)]
    ggplot(G_treattime_cor, aes(y=GOS.al, x=GOS.cp,color=color)) + 
      geom_point() + 
      geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
      theme_classic() +
      scale_color_manual(values=c("#004f7a", "#cc3d24")) +
      theme(legend.background = element_rect(fill="white",
                                             size=0.5, linetype="solid", 
                                             colour ="black")) +
      xlab("CP Coefficient") + ylab("Alum Coefficient") +
      theme(legend.position="none") 
  
  
  ##RSL treat##
  ##process##
    R_treat_cor=merge(R_treat, RSL_al_treat)
    R_treat_cor=R_treat_cor[,c(1,3)]
    names(R_treat_cor)[2] <- "RSL.al"
    R_treat_cor=merge(R_treat_cor, RSL_cp_treat)
    R_treat_cor=R_treat_cor[,c(1:2,4)]
    names(R_treat_cor)[3] <- "RSL.cp"
    R_treat_cor=R_treat_cor[!duplicated(R_treat_cor), ]
    ##stats and graphs##
    cor.test(R_treat_cor$RSL.al,R_treat_cor$RSL.cp, method="pearson")##sig, p<0.001, r=0.600
    
    ##graph if you want##
    R_treat_cor$effect=R_treat_cor$RSL.al*R_treat_cor$RSL.cp
    R_treat_cor=R_treat_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
    R_treat_cor=R_treat_cor[,c(1:3,5)]
    ggplot(R_treat_cor, aes(y=RSL.al, x=RSL.cp,color=color)) + 
      geom_point() + 
      geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
      theme_classic() +
      scale_color_manual(values=c("#004f7a", "#cc3d24")) +
      theme(legend.background = element_rect(fill="white",
                                             size=0.5, linetype="solid", 
                                             colour ="black")) +
      xlab("CP Coefficient") + ylab("Alum Coefficient") +
      theme(legend.position="none") 
  
  ##RSL treatxtime##
    ##process##
    R_treattime_cor=merge(R_treattime, RSL_al_treattime)
    R_treattime_cor=R_treattime_cor[,c(1,3)]
    names(R_treattime_cor)[2] <- "RSL.al"
    R_treattime_cor=merge(R_treattime_cor, RSL_cp_treattime)
    R_treattime_cor=R_treattime_cor[,c(1:2,4)]
    names(R_treattime_cor)[3] <- "RSL.cp"
    R_treattime_cor=R_treattime_cor[!duplicated(R_treattime_cor), ]
    ##stats and graphs##
    cor.test(R_treattime_cor$RSL.al,R_treattime_cor$RSL.cp, method="pearson")##sig, p<0.001, r=0.803
    
    ##graph if you want##
    R_treattime_cor$effect=R_treattime_cor$RSL.al*R_treattime_cor$RSL.cp
    R_treattime_cor=R_treattime_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
    R_treattime_cor=R_treattime_cor[,c(1:3,5)]
    ggplot(R_treattime_cor, aes(y=RSL.al, x=RSL.cp,color=color)) + 
      geom_point() + 
      geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
      theme_classic() +
      scale_color_manual(values=c("#004f7a", "#cc3d24")) +
      theme(legend.background = element_rect(fill="white",
                                             size=0.5, linetype="solid", 
                                             colour ="black")) +
      xlab("CP Coefficient") + ylab("Alum Coefficient") +
      theme(legend.position="none") 
  
  ##SAY treat##
    ##process##
    S_treat_cor=merge(S_treat, SAY_al_treat)
    S_treat_cor=S_treat_cor[,c(1,3)]
    names(S_treat_cor)[2] <- "SAY.al"
    S_treat_cor=merge(S_treat_cor, SAY_cp_treat)
    S_treat_cor=S_treat_cor[,c(1:2,4)]
    names(S_treat_cor)[3] <- "SAY.cp"
    S_treat_cor=S_treat_cor[!duplicated(S_treat_cor), ]
    ##stats and graphs##
    cor.test(S_treat_cor$SAY.al,S_treat_cor$SAY.cp, method="pearson")##sig, p<0.001, r=0.524
    
    ##graph if you want##
    S_treat_cor$effect=S_treat_cor$SAY.al*S_treat_cor$SAY.cp
    S_treat_cor=S_treat_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
    S_treat_cor=S_treat_cor[,c(1:3,5)]
    ggplot(S_treat_cor, aes(y=SAY.al, x=SAY.cp,color=color)) + 
      geom_point() + 
      geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
      theme_classic() +
      scale_color_manual(values=c("#004f7a", "#cc3d24")) +
      theme(legend.background = element_rect(fill="white",
                                             size=0.5, linetype="solid", 
                                             colour ="black")) +
      xlab("CP Coefficient") + ylab("Alum Coefficient") +
      theme(legend.position="none") 
  
  ##SAY treatxtime##
    ##process##
    S_treattime_cor=merge(S_treattime, SAY_al_treattime)
    S_treattime_cor=S_treattime_cor[,c(1,3)]
    names(S_treattime_cor)[2] <- "SAY.al"
    S_treattime_cor=merge(S_treattime_cor, SAY_cp_treattime)
    S_treattime_cor=S_treattime_cor[,c(1:2,4)]
    names(S_treattime_cor)[3] <- "SAY.cp"
    S_treattime_cor=S_treattime_cor[!duplicated(S_treattime_cor), ]
    ##stats and graphs##
    cor.test(S_treattime_cor$SAY.al,S_treattime_cor$SAY.cp, method="pearson")##sig, p<0.001, r=0.750
    
    ##graph if you want
    S_treattime_cor$effect=S_treattime_cor$SAY.al*S_treattime_cor$SAY.cp
    S_treattime_cor=S_treattime_cor %>% mutate(color = ifelse(effect > 0, "blue", "red"))
    S_treattime_cor=S_treattime_cor[,c(1:3,5)]
    ggplot(S_treattime_cor, aes(y=SAY.al, x=SAY.cp,color=color)) + 
      geom_point() + 
      geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
      theme_classic() +
      scale_color_manual(values=c("#004f7a", "#cc3d24")) +
      theme(legend.background = element_rect(fill="white",
                                             size=0.5, linetype="solid", 
                                             colour ="black")) +
      xlab("CP Coefficient") + ylab("Alum Coefficient") +
      theme(legend.position="none") 

