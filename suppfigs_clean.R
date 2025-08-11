##This is the code to create our supplementary figures for the injection experiment pub (minus the experimental design figure, which was made in BioRender)##
##last updated Aug 7th, 2025 by LEF##

##start by importing necessary packages##
library(tibble)
library(data.table)
library(dplyr)
library(plotrix)
library(ggpubr)
library(VennDiagram)
library(pheatmap)
library(ellipse)
library(RColorBrewer)
library(MASS)
library(car)

##first, data import/processing##
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



##next we need to parse and format data##
  ##alum subset
  readcount.alum <- readcount[metadata$Treatment %in% c("Alum", "PBS"),]
  metadata.alum <- metadata[metadata$Treatment %in% c("Alum", "PBS"),]
  metadata.alum$Treatment<- relevel(factor(metadata.alum$Treatment), ref = "PBS")
  metadata.alum$Population<- relevel(factor(metadata.alum$Population), ref = "SAY")
  names.al=readcount.alum$Fish_ID
  readcount.alum <- readcount.alum[,-1]
  geneIDs <- names(readcount.alum)
  readcount.alum=sapply(readcount.alum,as.numeric)
  
  ##cp subset
  readcount.cp <- readcount[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  metadata.cp <- metadata[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  metadata.cp$Treatment<- relevel(factor(metadata.cp$Treatment), ref = "PBS")
  metadata.cp$Population<- relevel(factor(metadata.cp$Population), ref = "SAY")
  names.cp=readcount.cp$Fish_ID
  readcount.cp <- readcount.cp[,-1]
  geneIDs <- names(readcount.cp)
  readcount.cp=sapply(readcount.cp,as.numeric)
  
  ##fibrosis subset
  names_fib=readcount$Fish_ID
  readcount_fib <- readcount[,-1]
  geneIDs <- names(readcount_fib)
  readcount_fib=sapply(readcount_fib,as.numeric)
  metadata_fib <- na.omit(metadata)
  metadata_fib  <- metadata_fib[metadata_fib$DPI < 90,]
  metadata_fib <- metadata_fib[metadata_fib$Treatment %in% c("Alum", "PBS", "Worm Protein"),]
  metadata_fib$Treatment<- relevel(factor(metadata_fib$Treatment), ref = "PBS")
  metadata_fib$Population<- relevel(factor(metadata_fib$Population), ref = "SAY")
  metadata_fib$Fib_Avg=as.numeric(metadata_fib$Fib_Avg)


## Next set some plotting functions we'll use repetitively##
##efinie a plotting function for proportional reads over time, split by population
timecourseplot <- function(genename, name,Metadata, readcounts, names, linescol = c(2,1)){
  ### Plot time course for a focal population
  # Pick a gene
  y <- readcounts[,genename]  
  y
  Totalreads <- rowSums(readcounts)
  Y <- cbind(focalgene = y, Totalreads = Totalreads)
  Y <- as.data.frame(Y)
  Y <- Y %>%
    mutate(pop = focalgene/Totalreads)
  names=as.matrix(names)
  colnames(names)[1] <- "Fish_ID"
  Y=cbind(names,Y)
  final=Y[,c(1,4)]
  final$Yvals <- log(final$pop + 0.0000001)
  final=merge(final,Metadata,by="Fish_ID")
  final_raw=final
  means=final[,c(3:6)]
  lines=means %>%
    group_by(Treatment,DPI,Population) %>%
    summarize(SE = sd(Yvals)/n(), Yvals=mean(Yvals))
  
  final$tempvar <- name
  
  ggplot(final, aes(x=DPI, y=Yvals, colour=Population, shape=Treatment)) +
    scale_color_manual(values=c('#C3A016FF','#58A787FF', '#0C1F4BFF')) +
    theme_classic()+
    geom_line(data=lines, aes(x=DPI, y=Yvals, colour=Population, linetype=Treatment), size=.6) +
    scale_linetype_manual(values=c("dashed", "solid")) +
    theme(
      plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"),
      strip.background = element_rect(fill="grey"),
      strip.text = element_text(size=15, colour="white", face="bold"))+
    facet_grid(.~tempvar)+
    geom_point(data=lines, aes(x=DPI, y=Yvals, colour=Population)) + geom_errorbar(data=lines, aes(ymax=Yvals+SE, ymin=Yvals-SE), width=2) + theme(legend.position="none")
}  

##define a second function without dividing by population
timecourseplot_nopop <- function(genename, name,Metadata, readcounts, names, linescol = c(2,1)){
  ### Plot time course for a focal population
  # Pick a gene
  y <- readcounts[,genename]  
  y
  Totalreads <- rowSums(readcounts)
  Y <- cbind(focalgene = y, Totalreads = Totalreads)
  Y <- as.data.frame(Y)
  Y <- Y %>%
    mutate(pop = focalgene/Totalreads)
  names=as.matrix(names)
  colnames(names)[1] <- "Fish_ID"
  Y=cbind(names,Y)
  final=Y[,c(1,4)]
  final$Yvals <- log(final$pop + 0.0000001)
  final=merge(final,Metadata,by="Fish_ID")
  final_raw=final
  means=final[,c(3:6)]
  lines=means %>%
    group_by(Treatment,DPI) %>%
    summarize(SE = sd(Yvals)/n(), Yvals=mean(Yvals))
  
  final$tempvar <- name
  
  ggplot(final, aes(x=DPI, y=Yvals, shape=Treatment)) + 
    theme_classic()+
    geom_line(data=lines, aes(x=DPI, y=Yvals, linetype=Treatment), size=.6) +
    scale_linetype_manual(values=c("dashed", "solid")) +
    theme(
      plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"),
      strip.background = element_rect(fill="grey"),
      strip.text = element_text(size=15, colour="white", face="bold"))+
    facet_grid(.~tempvar)+
    geom_point(data=lines, aes(x=DPI, y=Yvals)) + geom_errorbar(data=lines, aes(ymax=Yvals+SE, ymin=Yvals-SE), width=2) + theme(legend.position="none")
}  

##define our fibrosis plot function (split by population)##
fibplot <- function(genename, name,Metadata, readcounts, names, linescol = c(2,1)){
  ### Plot time course for a focal population
  # Pick a gene
  y <- readcounts[,genename]  
  y
  Totalreads <- rowSums(readcounts)
  Y <- cbind(focalgene = y, Totalreads = Totalreads)
  Y <- as.data.frame(Y)
  Y <- Y %>%
    mutate(pop = focalgene/Totalreads)
  names=as.matrix(names)
  colnames(names)[1] <- "Fish_ID"
  Y=cbind(names,Y)
  final=Y[,c(1,4)]
  final$Yvals <- log(final$pop + 0.0000001)
  final=merge(final,Metadata,by="Fish_ID")
  final_raw=final
  means=final[,c(3:4,11)]
  
  
  final$tempvar <- name
  
  ggplot(final, aes(x=Fib_L_count, y=Yvals, colour=Population, shape=Population)) +
    scale_shape_manual(values=c(1, 16, 17))+ 
    scale_color_manual(values=c('#C3A016FF','#58A787FF', '#0C1F4BFF')) +
    geom_point(alpha=0.2) + geom_smooth(method = "lm", aes(fill=Population))+
    scale_fill_manual(values=c('#C3A016FF','#58A787FF', '#0C1F4BFF')) +
    theme_classic() +theme(
      plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"),
      strip.background = element_rect(fill="grey"),
      strip.text = element_text(size=15, colour="white", face="bold"))+
    facet_grid(.~tempvar) + theme(legend.position="none")}  


##define our pop comparison plot
popcomp <- function(genename, name,Metadata, readcounts, names, linescol = c(2,1)){
  ### Plot time course for a focal population
  # Pick a gene
  y <- readcounts[,genename] 
  y
  Totalreads <- rowSums(readcounts)
  Y <- cbind(focalgene = y, Totalreads = Totalreads)
  Y <- as.data.frame(Y)
  Y <- Y %>%
    mutate(pop = focalgene/Totalreads)
  names=as.matrix(names)
  colnames(names)[1] <- "Fish_ID"
  Y=cbind(names,Y)
  final=Y[,c(1,4)]
  final$Yvals <- log(final$pop + 0.0000001)
  final=merge(final,Metadata,by="Fish_ID")
  final_raw=final
  means=final[,c(3:6)]
  lines=means %>%
    group_by(Treatment) %>%
    summarize(SE = sd(Yvals)/n(), Yvals=mean(Yvals))
  
  final$tempvar <- name
  
  ggplot(final, aes(x=Treatment, y=Yvals)) +geom_violin()+geom_boxplot(width=0.1)+
    theme_classic()+theme(
      plot.title = element_text(hjust = 0.5, color="black", size=14, face="bold"),
      strip.background = element_rect(fill="grey"),
      strip.text = element_text(size=15, colour="white", face="bold"))+
    facet_grid(.~tempvar)}  



##let's start with Supp Fig 2, wherein we create a multi panel of our correlations between cp/alum using ALL genes.

  ##start with treatment main effects##
  ##read in and process data##
  al.lfc=read.csv("Alum_GLM_matrix_results_fixed.csv")  
  al.lfc=al.lfc[,c(1,4,12)]  
  colnames(al.lfc)[3] <- "al.lfc"
  al.lfc=al.lfc[,c(1,3)]
  cp.lfc=read.csv("CestodeProt_GLM_matrix_results_fixed.csv")  
  cp.lfc=cp.lfc[,c(1,4,12)]  
  colnames(cp.lfc)[3] <- "cp.lfc"
  cp.lfc=cp.lfc[,c(1,3)]
  all=merge(al.lfc, cp.lfc, by="gene", all.x=TRUE)
  all=all[!duplicated(all), ] ##make sure to remove duplicates##
  shared.treat.plot <- all[,-1]
  rownames(shared.treat.plot) <- all[,1]
  ##here we're creating a column to color points based on which treatment they're more responsive to
  shared.treat.plot$effect=shared.treat.plot$al.lfc*shared.treat.plot$cp.lfc
  shared.treat.plot=shared.treat.plot %>% mutate(color = ifelse(effect > 0, "blue", "red"))
  shared.treat.plot=shared.treat.plot[,c(1:2,4)]
  
  ##and now we make that first plot, save it for multipanels
  treatcomp=ggplot(shared.treat.plot, aes(y=cp.lfc, x=al.lfc,color=color)) + 
    geom_point() + 
    geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
    theme_classic() +
    scale_color_manual(values=c("#004f7a", "#cc3d24")) +
    theme(legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black")) +
    xlab("Alum Coefficient") + ylab("Cestode Protein Coefficient") +
    theme(legend.position="none") 
  
  ##last we want the stats for display##
  cor.test(all$al.lfc,all$cp.lfc, method="pearson")
  #P<0.001, cor coefficient= 0.373
  
  
  ##same for treat*DPI##
  ##read in and process data##
  al.lfc=read.csv("Alum_GLM_matrix_results_fixed.csv")  
  al.lfc=al.lfc[,c(1,8,18)]  
  colnames(al.lfc)[3] <- "al.lfc"
  al.lfc=al.lfc[,c(1,3)]
  cp.lfc=read.csv("CestodeProt_GLM_matrix_results_fixed.csv")  
  cp.lfc=cp.lfc[,c(1,8,18)]  
  colnames(cp.lfc)[3] <- "cp.lfc"
  cp.lfc=cp.lfc[,c(1,3)]
  all=merge(al.lfc, cp.lfc, by="gene", all.x=TRUE)
  all=all[!duplicated(all), ] ##make sure to remove duplicates##
  shared.treatdpi.plot <- all[,-1]
  rownames(shared.treatdpi.plot) <- all[,1]
  ##here we're creating a column to color points based on which treatment they're more responsive to
  shared.treatdpi.plot$effect=shared.treatdpi.plot$al.lfc*shared.treatdpi.plot$cp.lfc
  shared.treatdpi.plot=shared.treatdpi.plot %>% mutate(color = ifelse(effect > 0, "blue", "red"))
  shared.treatdpi.plot=shared.treatdpi.plot[,c(1:2,4)]
  
  ##and now we make that first plot, save it for multipanels
  treatdpicomp=ggplot(shared.treatdpi.plot, aes(y=cp.lfc, x=al.lfc,color=color)) + 
    geom_point() + 
    geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
    theme_classic() +
    scale_color_manual(values=c("#004f7a", "#cc3d24")) +
    theme(legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black")) +
    xlab("Alum*Time Coefficient") + ylab("Cestode Protein*Time Coefficient") +
    theme(legend.position="none") 
  ##last we want the stats for display##
  cor.test(all$al.lfc,all$cp.lfc, method="pearson")
  #P<0.001, cor coefficient= 0.569
  
  
  ##same for treat*Geno##
  ##read in and process data##
  al.lfc=read.csv("Alum_GLM_matrix_results_fixed.csv")  
  al.lfc=al.lfc[,c(1,6,14:15)]  
  colnames(al.lfc)[3] <- "al.lfc.GS"
  colnames(al.lfc)[4] <- "al.lfc.RS"
  al.lfc=al.lfc[,c(1,3:4)]
  cp.lfc=read.csv("CestodeProt_GLM_matrix_results_fixed.csv")  
  cp.lfc=cp.lfc[,c(1,6,14:15)]   
  colnames(cp.lfc)[3] <-  "cp.lfc.GS"
  colnames(cp.lfc)[4] <- "cp.lfc.RS"
  cp.lfc=cp.lfc[,c(1,3:4)]
  all=merge(al.lfc, cp.lfc, by="gene", all.x=TRUE)
  all=all[!duplicated(all), ] ##make sure to remove duplicates!
  shared.treatpop.plot <- all[,-1]
  rownames(shared.treatpop.plot) <- all[,1]
  ##here we're creating a column to color points based on which treatment they're more responsive to
  shared.treatpop.plot$effectGS=shared.treatpop.plot$al.lfc.GS*shared.treatpop.plot$cp.lfc.GS
  shared.treatpop.plot=shared.treatpop.plot %>% mutate(colorGS = ifelse(effectGS > 0, "blue", "red"))
  shared.treatpop.plot$effectRS=shared.treatpop.plot$al.lfc.RS*shared.treatpop.plot$cp.lfc.RS
  shared.treatpop.plot=shared.treatpop.plot %>% mutate(colorRS = ifelse(effectRS > 0, "blue", "red"))
  shared.treatpop.plot=shared.treatpop.plot[,c(1:4,6,8)]
  
  ##and now we make plots, save it for multipanels
  treatpopGS=ggplot(shared.treatpop.plot, aes(y=cp.lfc.GS, x=al.lfc.GS,color=colorGS)) + 
    geom_point() + 
    geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
    theme_classic() +
    scale_color_manual(values=c("#004f7a", "#cc3d24")) +
    theme(legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black")) +
    xlab("Alum*Population (GvS) Coefficient") + ylab("Cestode Protein*Population (GvS) Coefficient") +
    theme(legend.position="none")
  ##last we want the stats for display##
  cor.test(all$al.lfc.GS,all$cp.lfc.GS, method="pearson")
  #P<0.001, cor coefficient= 0.422
  
  treatpopRS=ggplot(shared.treatpop.plot, aes(y=cp.lfc.RS, x=al.lfc.RS,color=colorRS)) + 
    geom_point() + 
    geom_abline(slope=0, intercept=0, linetype="solid") +geom_vline(xintercept=0)+
    theme_classic() +
    scale_color_manual(values=c("#004f7a", "#cc3d24")) +
    theme(legend.background = element_rect(fill="white",
                                           size=0.5, linetype="solid", 
                                           colour ="black")) +
    xlab("Alum*Population (RvS) Coefficient") + ylab("Cestode Protein*Population (RvS) Coefficient") +
    theme(legend.position="none")
  ##last we want the stats for display##
  cor.test(all$al.lfc.RS,all$cp.lfc.RS, method="pearson")
  #P<0.001, cor coefficient= 0.482
  
  ##and here we put all those panels together into one figure##
  ggarrange(treatcomp, treatdpicomp,treatpopGS,treatpopRS,
            widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 2,  common.legend=FALSE)


##supplementary figure 3- line plots over genes which were significantly differentially expressed in both alum and cp models
  ##we make one panel per gene and treatment, then pair them
  gene_i="ENSGACT00000019933"
  name="APP"
  APP_alum=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  APP_cp=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp,names.cp, linescol =c(2,1))
  APP=ggarrange(APP_alum,APP_cp,
                widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=FALSE)

  gene_i="ENSGACT00000010568"
  name="c4b"
  c4b_alum=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  c4b_cp=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp,names.cp, linescol =c(2,1))
  c4b=ggarrange(c4b_alum,c4b_cp,
                widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=FALSE)

  
  gene_i="ENSGACT00000003520"
  name="stat1b"
  stat1b_alum=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  stat1b_cp=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp,names.cp, linescol =c(2,1))
  stat1b=ggarrange(stat1b_alum,stat1b_cp,
                   widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=FALSE)
  
  
  gene_i="ENSGACT00000017488"
  name="SBNO2"
  SBNO2_alum=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  SBNO2_cp=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp,names.cp, linescol =c(2,1))
  SBNO2=ggarrange(SBNO2_alum,SBNO2_cp,
                  widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=FALSE)

  ##lastly we arrange the paired panels into one plot
  ggarrange(c4b,SBNO2,stat1b,APP,
            widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 2,  common.legend=TRUE)
  
 
  
  
   
##Next is Supp Figure 4 which shows the T-cell genes which were significantly differentially expressed as a result of treatment*time*pop interactions
  ##we make one figure per gene
  gene_i="ENSGACT00000017975"
  name="prf1.9"
  prf1=timecourseplot(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  
  gene_i="ENSGACT00000014058"
  name="rtkn2a"
  rtkn2a=timecourseplot(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  
  gene_i="ENSGACT00000027236"
  name="rfxap"
  rfxap=timecourseplot(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  
  ##and arrange them into one figure##
  ggarrange(prf1, rtkn2a, rfxap, 
            widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=3, nrow = 1,  common.legend=TRUE)
  
  

  
##now we make Supp Figure 5 which highlights immune genes which are divergent in response to cestode protein across our population specific models
  
  ### Need to reformat data to split by pop
  readcount.cp <- readcount[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  metadata.cp <- metadata[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  readcount.cp.GOS <- readcount.cp[metadata.cp$Population %in% c("GOS"),]
  metadata.cp.GOS <- metadata.cp[metadata.cp$Population %in% c("GOS"),]
  metadata.cp.GOS$treatmentColor <- as.numeric(factor(metadata.cp.GOS$Treatment))*-1 + 3
  metadata.cp.GOS$Treatment<- relevel(factor(metadata.cp.GOS$Treatment), ref = "PBS")
  names_GOS=readcount.cp.GOS$Fish_ID
  readcount.cp.GOS<- as.data.frame(readcount.cp.GOS)
  readcount.cp.GOS <- readcount.cp.GOS[,-1]
  readcount.cp.GOS=sapply(readcount.cp.GOS,as.numeric)
  
  readcount.cp <- readcount[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  metadata.cp <- metadata[metadata$Treatment %in% c("Worm Protein", "PBS"),]
  readcount.cp.RSL <- readcount.cp[metadata.cp$Population %in% c("RSL"),]
  metadata.cp.RSL <- metadata.cp[metadata.cp$Population %in% c("RSL"),]
  metadata.cp.RSL$treatmentColor <- as.numeric(factor(metadata.cp.RSL$Treatment))*-1 + 3
  metadata.cp.RSL$Treatment<- relevel(factor(metadata.cp.RSL$Treatment), ref = "PBS")
  names_RSL=readcount.cp.RSL$Fish_ID
  readcount.cp.RSL<- as.data.frame(readcount.cp.RSL)
  readcount.cp.RSL <- readcount.cp.RSL[,-1]
  readcount.cp.RSL=sapply(readcount.cp.RSL,as.numeric)
  
  ##and then make individual plots for each gene and population, and pair them
  gene_i="ENSGACT00000007690"
  name="AREL1"
  AREL1_GOS=popcomp(gene_i,name,metadata.cp.GOS,readcount.cp.GOS, names_GOS, linescol =c(2,1))
  AREL1_RSL=popcomp(gene_i,name,metadata.cp.RSL,readcount.cp.RSL, names_RSL, linescol =c(2,1))
  AREL1=ggarrange(AREL1_GOS,AREL1_RSL,
                  widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=TRUE) ##shared with fib##
  
  gene_i="ENSGACT00000023721"
  name="c7b"
  c7b_GOS=popcomp(gene_i,name,metadata.cp.GOS,readcount.cp.GOS, names_GOS, linescol =c(2,1))
  c7b_RSL=popcomp(gene_i,name,metadata.cp.RSL,readcount.cp.RSL, names_RSL, linescol =c(2,1))
  c7b=ggarrange(c7b_GOS,c7b_RSL,
                widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=TRUE) ##shared with fib##
  
  gene_i="ENSGACT00000008175"
  name="cdo1"
  cdo1_GOS=popcomp(gene_i,name,metadata.cp.GOS,readcount.cp.GOS, names_GOS, linescol =c(2,1))
  cdo1_RSL=popcomp(gene_i,name,metadata.cp.RSL,readcount.cp.RSL, names_RSL, linescol =c(2,1))
  cdo1=ggarrange(cdo1_GOS,cdo1_RSL,
                 widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=TRUE) ##shared with fib##
  
  
  gene_i="ENSGACT00000027025"
  name="dpf2"
  dpf2_GOS=popcomp(gene_i,name,metadata.cp.GOS,readcount.cp.GOS, names_GOS, linescol =c(2,1))
  dpf2_RSL=popcomp(gene_i,name,metadata.cp.RSL,readcount.cp.RSL, names_RSL, linescol =c(2,1))
  dpf2=ggarrange(dpf2_GOS,dpf2_RSL,
                 widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=TRUE) ##shared with fib##
  
  gene_i="ENSGACT00000012159"
  name="mxb"
  mxb_GOS=popcomp(gene_i,name,metadata.cp.GOS,readcount.cp.GOS, names_GOS, linescol =c(2,1))
  mxb_RSL=popcomp(gene_i,name,metadata.cp.RSL,readcount.cp.RSL, names_RSL, linescol =c(2,1))
  mxb=ggarrange(mxb_GOS,mxb_RSL,
                widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=TRUE) ##shared with fib##
  
  ##finally we arrange them
  ggarrange(c7b,mxb, AREL1, cdo1, dpf2,
            widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow =3,  common.legend=TRUE) ##shared with fib##

  
##last we do supplementary figure 6, which is our trajectory plot##
  
  # First choose which genes to focus on, picking all with significant LFCs
  # Get a list of file names
  Contrasts <- list.files("GLM_output")
  Nfiles <- length(Contrasts)
  
  sig.genes <- c()
  for(i in 1:Nfiles){
    file1 <- Contrasts[i]
    file1 <- fread(paste("GLM_Output/",file1, sep = ""))
    file1.sig <- file1[file1$padj < 0.01,] ##make sure to set the correct alpha##
    dim(file1.sig)      
    sig.genes <- c(sig.genes, file1.sig$gene)
  }
  length(sig.genes)
  sig.genes <- unique(sig.genes)
  
  #Read in raw read data
  metadata2 <- read.csv("ExpDesign_noWA.csv")
  dim(metadata2)
  metadata2 <- metadata2[c(1:312),]
  readcounts <- fread("reads_noWA.csv", data.table = F)
  dim(readcounts)
  
  readcounts <- readcounts[, names(readcounts) %in% sig.genes]
  dim(readcounts)
  
  ##run PCA##
  PCA <- prcomp(readcounts, scale = T, center= T)
  
  ######## DAPC just for PBS and  Alum
  ### Try DAPC
  head(metadata2)
  metasubset <- metadata2[metadata2$Treatment %in% c("PBS", "Alum"),]
  PCAsubset <- PCA$x[ metadata2$Treatment %in% c("PBS", "Alum") , 1:100]
  groups <- paste(metasubset$Treatment, metasubset$Population, metasubset$DPI, sep = "_")
  DFA <- lda(groups ~ PCAsubset)
  DFA.p <- predict(DFA)
  {
    x1means <- tapply(DFA.p$x[,1], groups, mean)
    x2means <- tapply(DFA.p$x[,2], groups, mean)
    x3means <- tapply(DFA.p$x[,3], groups, mean)
    x4means <- tapply(DFA.p$x[,4], groups, mean)
    x5means <- tapply(DFA.p$x[,5], groups, mean)
    
    x1sds <- tapply(DFA.p$x[,1], groups, sd)
    x2sds <- tapply(DFA.p$x[,2], groups, sd)
    x3sds <- tapply(DFA.p$x[,3], groups, sd)
    x4sds <- tapply(DFA.p$x[,4], groups, sd)
    x5sds <- tapply(DFA.p$x[,5], groups, sd)
    
    Ns <- summary(factor(groups))
    x1ses <- x1sds/sqrt(Ns)
    x2ses <- x2sds/sqrt(Ns)
    x3ses <- x3sds/sqrt(Ns)
    x4ses <- x4sds/sqrt(Ns)
    x5ses <- x5sds/sqrt(Ns)
    
    lowerCI1 <- x1means - x1ses
    lowerCI2 <- x2means - x2ses
    lowerCI3 <- x3means - x3ses
    lowerCI4 <- x4means - x4ses
    lowerCI5 <- x5means - x5ses
    
    upperCI1 <- x1means + x1ses
    upperCI2 <- x2means + x2ses
    upperCI3 <- x3means + x3ses
    upperCI4 <- x4means + x4ses
    upperCI5 <- x5means + x5ses
    
    summarytable <- data.frame(x1means, x2means, x3means, x4means, x5means, lowerCI1, lowerCI2, lowerCI3, lowerCI4, lowerCI5, upperCI1,upperCI2, upperCI3, upperCI4, upperCI5)
    summarytable$Treatment <- c(  rep("Alum", 12), rep("PBS", 12))
    summarytable$Day <- c(rep(c(1,10,42,90), 6))
    summarytable$Pop <- rep(c("GOS", "GOS", "GOS","GOS","RSL", "RSL", "RSL", "RSL", "SAY","SAY", "SAY" ,"SAY"), 2)
    summarytable$colorplot <- as.numeric(factor(summarytable$Pop))
    summarytable$colorplot[summarytable$colorplot == 3] <- 4
  }
  
  
  # Make the multi-panel plot (axes may be slightly different than those shown in pub, but same general code)
  {
    par(mfrow = c(1,3), mar = c(5,5,2,1))
    tempDFA <- DFA.p$x[metasubset$Population == "SAY",]
    plot(tempDFA[,2] ~ tempDFA[,1], pch = 1, cex = 0.5, col = rgb(0,0,0,1), xlim = c(-4.5,5), ylim = c(-12,-2), xlab = "DFA axis 1", ylab = "DFA axis 2", cex.lab = 1.3, main = c("SAY"))
    summarytable2 <- summarytable[-c(1:8,13:20,22:24),] # remove Saline on later days
    summarytemp <- summarytable2[summarytable2$Pop == "SAY",]
    points(summarytemp$x1means, summarytemp$x2means, pch = c(16,16,16,16,1), col = "#C3A016FF", cex = 2.5, lwd = 3)
    arrows(x0 = summarytemp$x1means[c(5,1:3)], y0 = summarytemp$x2means[c(5,1:3)], x1 = summarytemp$x1means[c(1:4)], y1 = summarytemp$x2means[c(1:4)], lwd = 3, length = 0.2, col = "#C3A016FF")
    text(summarytemp$x1means+0.2, summarytemp$x2means, c( "1", "10","42","90", "PBS"),pos = 4)
    arrows(x0 = summarytemp$lowerCI1[c(1:5)], y0 = summarytemp$x2means[c(1:5)], x1 = summarytemp$upperCI1[c(1:5)], y1 = summarytemp$x2means[c(1:5)], code = 3, angle = 90, length  = 0.02, col = "grey50")
    arrows(x0 = summarytemp$x1means[c(1:5)], y0 = summarytemp$lowerCI2[c(1:5)], x1 = summarytemp$x1means[c(1:5)], y1 = summarytemp$upperCI2[c(1:5)], code = 3, angle = 90, length  = 0.02, col = "grey50")
    
    tempDFA <- DFA.p$x[metasubset$Population == "GOS",]
    plot(tempDFA[,2] ~ tempDFA[,1], pch = 1, cex = 0.5, col = rgb(0,0,0,1),  xlim = c(-13,-1), ylim = c(0.5,6.5), xlab = "DFA axis 1", ylab = "DFA axis 2", cex.lab = 1.3, main = c("GOS"))
    summarytable2 <- summarytable[-c(5:12,14:24),] # remove Saline on later days
    summarytemp <- summarytable2[summarytable2$Pop == "GOS",]
    points(summarytemp$x1means, summarytemp$x2means, pch = c(16,16,16,16, 1), col = "#58A787FF", cex = 2.5, lwd = 3)
    arrows(x0 = summarytemp$x1means[c(5,1:3)], y0 = summarytemp$x2means[c(5,1:3)], x1 = summarytemp$x1means[c(1:4)], y1 = summarytemp$x2means[c(1:4)], lwd = 3, length = 0.2, col = "#58A787FF")
    text(summarytemp$x1means+0.2, summarytemp$x2means, c( "1", "10","42", "90", "PBS"),pos = 4)
    arrows(x0 = summarytemp$lowerCI1[c(1:5)], y0 = summarytemp$x2means[c(1:5)], x1 = summarytemp$upperCI1[c(1:5)], y1 = summarytemp$x2means[c(1:5)], code = 3, angle = 90, length  = 0.02, col = "grey50")
    arrows(x0 = summarytemp$x1means[c(1:5)], y0 = summarytemp$lowerCI2[c(1:5)], x1 = summarytemp$x1means[c(1:5)], y1 = summarytemp$upperCI2[c(1:5)], code = 3, angle = 90, length  = 0.02, col = "grey50")
    
    tempDFA <- DFA.p$x[metasubset$Population == "RSL",]
    plot(tempDFA[,2] ~ tempDFA[,1], pch = 1, cex = 0.5, col = rgb(0,0,0,1), xlim = c(2,15), ylim = c(-2,8), xlab = "DFA axis 1", ylab = "DFA axis 2", cex.lab = 1.3, main = c("RSL"))
    summarytable2 <- summarytable[-c(1:4,9:16,18:24),] # remove Saline on later days
    summarytemp <- summarytable2[summarytable2$Pop == "RSL",]
    points(summarytemp$x1means, summarytemp$x2means, pch = c(16,16,16,16,1), col = "#0C1F4BFF", cex = 2.5, lwd = 3)
    arrows(x0 = summarytemp$x1means[c(5,1:3)], y0 = summarytemp$x2means[c(5,1:3)], x1 = summarytemp$x1means[c(1:4)], y1 = summarytemp$x2means[c(1:4)], lwd = 3, length = 0.2, col = "#0C1F4BFF")
    text(summarytemp$x1means+0.2, summarytemp$x2means, c( "1", "10","42","90", "PBS"),pos = 4)
    arrows(x0 = summarytemp$lowerCI1[c(1:5)], y0 = summarytemp$x2means[c(1:5)], x1 = summarytemp$upperCI1[c(1:5)], y1 = summarytemp$x2means[c(1:5)], code = 3, angle = 90, length  = 0.02, col = "grey50")
    arrows(x0 = summarytemp$x1means[c(1:5)], y0 = summarytemp$lowerCI2[c(1:5)], x1 = summarytemp$x1means[c(1:5)], y1 = summarytemp$upperCI2[c(1:5)], code = 3, angle = 90, length  = 0.02, col = "grey50")
    
  }
  
  
  
  
  ######## DAPC just for PBS and Worm protein
  ### Try DAPC
  head(metadata2)
  metasubset <- metadata2[metadata2$Treatment %in% c("PBS", "Worm Protein"),]
  PCAsubset <- PCA$x[ metadata2$Treatment %in% c("PBS", "Worm Protein") , 1:100]
  groups <- paste(metasubset$Treatment, metasubset$Population, metasubset$DPI, sep = "_")
  DFA <- lda(groups ~ PCAsubset)
  DFA.p <- predict(DFA)
  {
    x1means <- tapply(DFA.p$x[,1], groups, mean)
    x2means <- tapply(DFA.p$x[,2], groups, mean)
    x3means <- tapply(DFA.p$x[,3], groups, mean)
    x4means <- tapply(DFA.p$x[,4], groups, mean)
    x5means <- tapply(DFA.p$x[,5], groups, mean)
    
    x1sds <- tapply(DFA.p$x[,1], groups, sd)
    x2sds <- tapply(DFA.p$x[,2], groups, sd)
    x3sds <- tapply(DFA.p$x[,3], groups, sd)
    x4sds <- tapply(DFA.p$x[,4], groups, sd)
    x5sds <- tapply(DFA.p$x[,5], groups, sd)
    
    Ns <- summary(factor(groups))
    x1ses <- x1sds/sqrt(Ns)
    x2ses <- x2sds/sqrt(Ns)
    x3ses <- x3sds/sqrt(Ns)
    x4ses <- x4sds/sqrt(Ns)
    x5ses <- x5sds/sqrt(Ns)
    
    lowerCI1 <- x1means - x1ses
    lowerCI2 <- x2means - x2ses
    lowerCI3 <- x3means - x3ses
    lowerCI4 <- x4means - x4ses
    lowerCI5 <- x5means - x5ses
    
    upperCI1 <- x1means + x1ses
    upperCI2 <- x2means + x2ses
    upperCI3 <- x3means + x3ses
    upperCI4 <- x4means + x4ses
    upperCI5 <- x5means + x5ses
    
    summarytable <- data.frame(x1means, x2means, x3means, x4means, x5means, lowerCI1, lowerCI2, lowerCI3, lowerCI4, lowerCI5, upperCI1,upperCI2, upperCI3, upperCI4, upperCI5)
    summarytable$Treatment <- c( rep("PBS", 12), rep("WP", 12))
    summarytable$Day <- c(rep(c(1,10,42,90), 6))
    summarytable$Pop <- rep(c("GOS", "GOS", "GOS","GOS","RSL", "RSL", "RSL", "RSL", "SAY", "SAY", "SAY", "SAY"), 2)
    summarytable$colorplot <- as.numeric(factor(summarytable$Pop))
    summarytable$colorplot[summarytable$colorplot == 3] <- 4
  }
  
  
  
  # Make the plot
  {
    par(mfrow = c(1,3), mar = c(5,5,2,1))
    
    tempDFA <- DFA.p$x[metasubset$Population == "SAY",]
    plot(tempDFA[,2] ~ tempDFA[,1], pch = 1, cex = 0.5, col = rgb(0,0,0,1), xlim = c(-10, -1.5), ylim = c(-10, -1.5), xlab = "DFA axis 1", ylab = "DFA axis 2", cex.lab = 1.3, main = c("SAY"))
    summarytable2 <- summarytable[-c(1:8,10:12),] # remove Saline on later days
    summarytemp <- summarytable2[summarytable2$Pop == "SAY",]
    points(summarytemp$x1means, summarytemp$x2means, pch = c(1,16,16,16,16), col = "#C3A016FF", cex = 2.5, lwd = 3)
    arrows(x0 = summarytemp$x1means[c(1:4)], y0 = summarytemp$x2means[c(1:4)], x1 = summarytemp$x1means[c(2:5)], y1 = summarytemp$x2means[c(2:5)], lwd = 3, length = 0.2, col = "#C3A016FF")
    text(summarytemp$x1means+0.2, summarytemp$x2means, c("PBS", "1", "10","42","90"),pos = 4)
    arrows(x0 = summarytemp$lowerCI1[c(1:5)], y0 = summarytemp$x2means[c(1:5)], x1 = summarytemp$upperCI1[c(1:5)], y1 = summarytemp$x2means[c(1:5)], code = 3, angle = 90, length  = 0.02, col = "grey50")
    arrows(x0 = summarytemp$x1means[c(1:5)], y0 = summarytemp$lowerCI2[c(1:5)], x1 = summarytemp$x1means[c(1:5)], y1 = summarytemp$upperCI2[c(1:5)], code = 3, angle = 90, length  = 0.02, col = "grey50")
    
    tempDFA <- DFA.p$x[metasubset$Population == "GOS",]
    plot(tempDFA[,2] ~ tempDFA[,1], pch = 1, cex = 0.5, col = rgb(0,0,0,1), xlim = c(-8, 0), ylim = c(1.5, 14), xlab = "DFA axis 1", ylab = "DFA axis 2", cex.lab = 1.3, main = c("GOS"))
    summarytable2 <- summarytable[-c(2:12),] # remove PBS
    summarytemp <- summarytable2[summarytable2$Pop == "GOS",]
    points(summarytemp$x1means, summarytemp$x2means, pch = c(1,16,16,16,16),  col = "#58A787FF", cex = 2.5, lwd = 3)
    arrows(x0 = summarytemp$x1means[c(1:4)], y0 = summarytemp$x2means[c(1:4)], x1 = summarytemp$x1means[c(2:5)], y1 = summarytemp$x2means[c(2:5)], lwd = 3, length = 0.2, col = "#58A787FF")
    text(summarytemp$x1means+0.2, summarytemp$x2means, c("PBS", "1", "10","42","90"),pos = 4)
    arrows(x0 = summarytemp$lowerCI1[c(1:5)], y0 = summarytemp$x2means[c(1:5)], x1 = summarytemp$upperCI1[c(1:5)], y1 = summarytemp$x2means[c(1:5)], code = 3, angle = 90, length  = 0.02, col = "grey50")
    arrows(x0 = summarytemp$x1means[c(1:5)], y0 = summarytemp$lowerCI2[c(1:5)], x1 = summarytemp$x1means[c(1:45)], y1 = summarytemp$upperCI2[c(1:5)], code = 3, angle = 90, length  = 0.02, col = "grey50")
    
    tempDFA <- DFA.p$x[metasubset$Population == "RSL",]
    plot(tempDFA[,2] ~ tempDFA[,1], pch = 1, cex = 0.5, col = rgb(0,0,0,1), xlim = c(4, 15), ylim = c(-4, 3), xlab = "DFA axis 1", ylab = "DFA axis 2", cex.lab = 1.3, main = c("RSL"))
    summarytable2 <- summarytable[-c(1:4,6:12),] # remove Saline on later days
    summarytemp <- summarytable2[summarytable2$Pop == "RSL",]
    points(summarytemp$x1means, summarytemp$x2means, pch = c(1,16,16,16,16), col = "#0C1F4BFF", cex = 2.5, lwd = 3)
    text(summarytemp$x1means+0.2, summarytemp$x2means, c("PBS", "1", "10","42","90"),pos = 4)
    arrows(x0 = summarytemp$x1means[c(1:4)], y0 = summarytemp$x2means[c(1:4)], x1 = summarytemp$x1means[c(2:5)], y1 = summarytemp$x2means[c(2:5)], lwd = 3, length = 0.2, col = "#0C1F4BFF")
    arrows(x0 = summarytemp$lowerCI1[c(1:5)], y0 = summarytemp$x2means[c(1:5)], x1 = summarytemp$upperCI1[c(1:5)], y1 = summarytemp$x2means[c(1:5)], code = 3, angle = 90, length  = 0.02, col = "grey50")
    arrows(x0 = summarytemp$x1means[c(1:5)], y0 = summarytemp$lowerCI2[c(1:5)], x1 = summarytemp$x1means[c(1:5)], y1 = summarytemp$upperCI2[c(1:5)], code = 3, angle = 90, length  = 0.02, col = "grey50")
  }
  
  
  
  