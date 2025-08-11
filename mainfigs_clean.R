##This is the code to create our main figures for the injection experiment pub (minus the venn diagram, which was hand made)##
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
  

  
  
  
##Ok let's dive into our figures, starting with Fig 2
  ##here we're going to make correlation plots of response to alum versus cestode protein
  ##we'll do this for main treatment effects, treatment*time effects, and treatment*population effects
  
  ##start with treatment main effects##
  ##read in and process data##
  al.lfc=read.csv("Alum_GLM_matrix_results_fixed.csv")  
  al.lfc=al.lfc[,c(1,4,12)]  
  colnames(al.lfc)[3] <- "al.lfc"
  al.lfc.sig=al.lfc[al.lfc$Treat_P < .01, ] ##selecting only sig genes##
  al.lfc.sig=al.lfc.sig[,c(1,3)]
  cp.lfc=read.csv("CestodeProt_GLM_matrix_results_fixed.csv")  
  cp.lfc=cp.lfc[,c(1,4,12)]  
  colnames(cp.lfc)[3] <- "cp.lfc"
  cp.lfc.sig=cp.lfc[cp.lfc$Treat_P < .01, ] ##selecting only sig genes##
  cp.lfc.sig=cp.lfc.sig[,c(1,3)]
  al.sig=merge(al.lfc.sig, cp.lfc, by="gene", all.x=TRUE)
  al.sig=al.sig[,c(1:2,4)]
  cp.sig=merge(cp.lfc.sig, al.lfc, by = "gene",all.x=TRUE)
  cp.sig=cp.sig[,c(1,4,2)]
  all.sig=rbind(al.sig,cp.sig)
  all.sig=all.sig[!duplicated(all.sig), ] ##make sure to remove duplicates##
  shared.treat.plot <- all.sig[,-1]
  rownames(shared.treat.plot) <- all.sig[,1]
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
  cor.test(all.sig$al.lfc,all.sig$cp.lfc, method="pearson")
  #P<0.001, cor coefficient= 0.411
  
  ##same for treat*DPI##
  ##read in and process data##
  al.lfc=read.csv("Alum_GLM_matrix_results_fixed.csv")  
  al.lfc=al.lfc[,c(1,8,18)]  
  colnames(al.lfc)[3] <- "al.lfc"
  al.lfc.sig=al.lfc[al.lfc$Treat_DPI_P < .01, ]  ##selecting only sig genes##
  al.lfc.sig=al.lfc.sig[,c(1,3)]
  cp.lfc=read.csv("CestodeProt_GLM_matrix_results_fixed.csv")  
  cp.lfc=cp.lfc[,c(1,8,18)]  
  colnames(cp.lfc)[3] <- "cp.lfc"
  cp.lfc.sig=cp.lfc[cp.lfc$Treat_DPI_P < .01, ] ##selecting only sig genes##
  cp.lfc.sig=cp.lfc.sig[,c(1,3)]
  al.sig=merge(al.lfc.sig, cp.lfc, by="gene", all.x=TRUE)
  al.sig=al.sig[,c(1:2,4)]
  cp.sig=merge(cp.lfc.sig, al.lfc, by = "gene",all.x=TRUE)
  cp.sig=cp.sig[,c(1,4,2)]
  all.sig=rbind(al.sig,cp.sig)
  all.sig=all.sig[!duplicated(all.sig), ] ##make sure to remove duplicates##
  shared.treatdpi.plot <- all.sig[,-1]
  rownames(shared.treatdpi.plot) <- all.sig[,1]
  ##here we're creating a column to color points based on which treatment they're more responsive to
  shared.treatdpi.plot$effect=shared.treatdpi.plot$al.lfc*shared.treatdpi.plot$cp.lfc
  shared.treatdpi.plot=shared.treatdpi.plot %>% mutate(color = ifelse(effect > 0, "blue", "red"))
  shared.treatdpi.plot=shared.treatdpi.plot[,c(1:2,4)]
  
  ##and now we make that plot, save it for multipanels
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
  cor.test(all.sig$al.lfc,all.sig$cp.lfc, method="pearson")
  #P<0.001, cor coefficient= 0.702
  
  
  ##same for treat*Geno##
  ##read in and process data##
  al.lfc=read.csv("Alum_GLM_matrix_results_fixed.csv")  
  al.lfc=al.lfc[,c(1,6,14:15)]  
  colnames(al.lfc)[3] <- "al.lfc.GS"
  colnames(al.lfc)[4] <- "al.lfc.RS"
  al.lfc.sig=al.lfc[al.lfc$Pop_Trt_P < .01, ] ##selecting only sig genes##
  al.lfc.sig=al.lfc.sig[,c(1,3:4)]
  cp.lfc=read.csv("CestodeProt_GLM_matrix_results_fixed.csv")  
  cp.lfc=cp.lfc[,c(1,6,14:15)]   
  colnames(cp.lfc)[3] <-  "cp.lfc.GS"
  colnames(cp.lfc)[4] <- "cp.lfc.RS"
  cp.lfc.sig=cp.lfc[cp.lfc$Pop_Trt_P < .01, ] ##selecting only sig genes##
  cp.lfc.sig=cp.lfc.sig[,c(1,3:4)]
  al.sig=merge(al.lfc.sig, cp.lfc, by="gene", all.x=TRUE)
  al.sig=al.sig[,c(1:3,5:6)]
  cp.sig=merge(cp.lfc.sig, al.lfc, by = "gene",all.x=TRUE)
  cp.sig=cp.sig[,c(1,5:6,2:3)]
  all.sig=rbind(al.sig,cp.sig)
  all.sig=all.sig[!duplicated(all.sig), ] ##make sure to remove duplicates!
  shared.treatpop.plot <- all.sig[,-1]
  rownames(shared.treatpop.plot) <- all.sig[,1]
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
  cor.test(all.sig$al.lfc.GS,all.sig$cp.lfc.GS, method="pearson")
  #P<0.001, cor coefficient= 0.604
  
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
  cor.test(all.sig$al.lfc.RS,all.sig$cp.lfc.RS, method="pearson")
  #P<0.001, cor coefficient= 0.733
  
  ##and here we put all those panels together into one figure##
  ggarrange(treatcomp, treatdpicomp,treatpopGS,treatpopRS,
            widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 2,  common.legend=FALSE)
  
  
  
  
  
##Great, now we'll move on to figure 3 where we look at genes which are shared between alum main effect and cp*time (statistically overrepresented overlap)
  ##we make two panels for each gene, alum and cestode protein
  ##ANXA4
  gene_i="ENSGACT00000022600"
  name="anxa4"
  anxa4_alum=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  anxa4_cp=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp,names.cp, linescol =c(2,1))
  anxa4=ggarrange(anxa4_alum,anxa4_cp,
                  widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=FALSE) ##shared with fib##
  
  ##BCAP31
  gene_i="ENSGACT00000000357"
  name="bcap31"
  bcap31_alum=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  bcap31_cp=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp,names.cp, linescol =c(2,1))
  bcap31=ggarrange(bcap31_alum,bcap31_cp,
                   widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=FALSE)
  
  ##CCDC22
  gene_i="ENSGACT00000017674"
  name="ccdc22"
  ccdc22_alum=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  ccdc22_cp=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp,names.cp, linescol =c(2,1))
  ccdc22=ggarrange(ccdc22_alum,ccdc22_cp,
                   widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=FALSE) ##shared with fib##
  
  ##LAMP1
  gene_i="ENSGACT00000020286"
  name="LAMP1"
  LAMP1_alum=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  LAMP1_cp=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp,names.cp, linescol =c(2,1))
  LAMP1=ggarrange(LAMP1_alum,LAMP1_cp,
                  widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=FALSE) ##shared with fib##
  
  ##MRC1B
  gene_i="ENSGACT00000022556"
  name="mrc1b"
  mrc1b_alum=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  mrc1b_cp=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp,names.cp, linescol =c(2,1))
  mrc1b=ggarrange(mrc1b_alum,mrc1b_cp,
                  widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=FALSE) ##shared with fib/fib*pop##
  
  ##UFD1L
  gene_i="ENSGACT00000005239"
  name="ufd1l"
  ufd1l_alum=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  ufd1l_cp=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp,names.cp, linescol =c(2,1))
  ufd1l=ggarrange(ufd1l_alum,ufd1l_cp,
                  widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=FALSE) ##shared with fib##
  
  ##KLF2B
  gene_i="ENSGACT00000022354"
  name="klf2b"
  klf2b_alum=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  klf2b_cp=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp,names.cp, linescol =c(2,1))
  klf2b=ggarrange(klf2b_alum,klf2b_cp,
                  widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=FALSE) ##shared with fib##
  
  ##CXCL19
  gene_i="ENSGACT00000025278"
  name="cxcl19"
  cxcl19_alum=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
  cxcl19_cp=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp,names.cp, linescol =c(2,1))
  cxcl19=ggarrange(cxcl19_alum,cxcl19_cp,
                   widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=FALSE) ##shared with fib##
  
  ##then put them all together
  ggarrange(anxa4,bcap31,ccdc22,cxcl19, klf2b,LAMP1,mrc1b,ufd1l,
            widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 4,  common.legend=FALSE) ##all but bcap31 are shared with fib##
  

  
  
##Next we'll do Figure 4 where we look at how those genes plotted in figure 3 relate to fibrosis##
  ##plot some of our overlapping fibrosis genes which are alum main/cp*time##
  ##note all of these are generally positively associated with fibrosis
  ##and all show upregulation in alum, but more downregulation in CP
  
  ##here we make one panel per gene because our fibrosis model did not include a treatment term##
  gene_i="ENSGACT00000022600"
  name="anxa4"
  anxa4=fibplot(gene_i,name,metadata_fib,readcount_fib,names_fib, linescol =c(2,1))
  
  gene_i="ENSGACT00000017674"
  name="ccdc22"
  ccdc22=fibplot(gene_i,name,metadata_fib,readcount_fib,names_fib, linescol =c(2,1))
  
  gene_i="ENSGACT00000025278"
  name="cxcl19"
  cxcl19=fibplot(gene_i,name,metadata_fib,readcount_fib,names_fib, linescol =c(2,1))
  
  gene_i="ENSGACT00000022556"
  name="mrc1b"
  mrc1b=fibplot(gene_i,name,metadata_fib,readcount_fib,names_fib, linescol =c(2,1))
  
  gene_i="ENSGACT00000005239"
  name="ufd1l"
  ufd1l=fibplot(gene_i,name,metadata_fib,readcount_fib,names_fib, linescol =c(2,1))
  
  gene_i="ENSGACT00000020286"
  name="LAMP1"
  LAMP1=fibplot(gene_i,name,metadata_fib,readcount_fib,names_fib, linescol =c(2,1))
  
  gene_i="ENSGACT00000022354"
  name="klf2b"
  klf2b=fibplot(gene_i,name,metadata_fib,readcount_fib,names_fib, linescol =c(2,1))
  
  ##we then combine figures into one plot##
  ggarrange(anxa4, ccdc22, cxcl19, klf2b, LAMP1, mrc1b, ufd1l,
            widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=4, nrow = 2,  common.legend=TRUE)
  
##Now we'll move on to treatment specific effects, starting with alum and Figure 5 (Tcells) ##
  ##with alum we see clear picture of tcells
  ##and this overlaps some with fibrosis
    ##start with heatmap of alum and fibrosis main effects on tcells
    ##our heatmap will display model coefficients for alum and fibrosis main effects for our t-cell genes which are sig for alum
    ##we will order it manually by alum coefficients##
    
    ##start by reading in and processing data
    tcell_main=read.csv("tcell_main_GOI.csv") ##list of tcell genes hand curated; also includes our annotation related to how genes regulate Tcells##
    alum=read.csv("Alum_GLM_matrix_results_fixed.csv")  ##alum model results
    alum_C=alum[,c(1,12)]  
    fib=read.csv("Fibrosis_GLM_matrix_results_fixed.csv") ##fibrosis model results##
    fib_C=fib[,c(1,8)]
    ##merge
    tcell_main=merge(tcell_main,alum_C,by.x="Transcript",by.y="gene")  
    tcell_main=merge(tcell_main,fib_C,by.x="Transcript",by.y="gene")  
    ##order them and process##
    tcell_main=tcell_main[order(-tcell_main$Treat_C),] 
    tcell_main$Regulate=as.factor(tcell_main$Regulate)
    tcell_heatmap=tcell_main[,c(2,4:5,3)]  
    tcell <- tcell_heatmap[,-1]
    rownames(tcell) <- tcell_heatmap[,1]
    tcell=as.matrix(tcell[,c(1:2)])
    rownames(tcell) <- tcell_heatmap[,1]
    tcell_annos=tcell_main[,c(1,3)]
    tcell_annos$blank="y" ##doing this cause pheatmap doesn't like having just one annotation column##
    tcell_annos=tcell_annos[,c(2:3)]
    rownames(tcell_annos) <- tcell_heatmap[,1]
    
    
    ##set a color palatte
    paletteLength <- 50
    myColor <- colorRampPalette(c("#762a83", "white", "#1b7837"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    ##basically we're setting white as 0##
    myBreaks <- c(seq(min(tcell), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(tcell)/paletteLength, max(tcell), length.out=floor(paletteLength/2)))
    
    # make our heatmap (will be modified in illustrator)
    pheatmap(tcell,cluster_cols=FALSE,cluster_rows=FALSE, annotation_row = tcell_annos, 
                  annotation_colors = list(
                    Regulate = c(NEG = "#5977ff", POS = "#f74747", AP = "#edf8b1",REG = "grey"),
                    blank = c(y="grey")),
                  color=myColor, breaks=myBreaks)
    
    ##let's test whether alum and fibrosis effects on these genes are congruent statistically##
    tcell=as.data.frame(tcell)
    cor.test(tcell$Treat_C,tcell$Fib_C, method="pearson") #p<0.001, r=.666
    ##ok so we see mainly congruence! (+ about 65% of alum responsive genes are also sig associated with fibrosis/fibrosis*pop)
    
  ##let's go ahead and plot out our ones where pop*fib is sig too for a multi-panel figure##
    ##we do this because the fibrosis main coefficient is already plotted, so we don't need to plot those again##
    ##(here we're just trying to show population trends)##
    
    ##like above we make one panel per gene##
    gene_i="ENSGACT00000017961"
    name="tnfaip8l2b"
    tnfaip8l2b=fibplot(gene_i,name,metadata_fib,readcount_fib,names_fib, linescol =c(2,1))
    
    gene_i="ENSGACT00000007096"
    name="cyldl"
    cyldl=fibplot(gene_i,name,metadata_fib,readcount_fib,names_fib, linescol =c(2,1))
    
    gene_i="ENSGACT00000021776"
    name="irf4a"
    irf4a=fibplot(gene_i,name,metadata_fib,readcount_fib,names_fib, linescol =c(2,1))
    
    gene_i="ENSGACT00000010177"
    name="rbx1"
    rbx1=fibplot(gene_i,name,metadata_fib,readcount_fib,names_fib, linescol =c(2,1))
    
    gene_i="ENSGACT00000012461"
    name="prkd2"
    prkd2=fibplot(gene_i,name,metadata_fib,readcount_fib,names_fib, linescol =c(2,1))
    
    
    ##and then arrange them together##
    ggarrange(cyldl,irf4a, prkd2, rbx1, tnfaip8l2b,
              widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 3,  common.legend=TRUE)
    
    
    
    
    
    
    
##Ok we have one more alum/t-cell figure to make, figure 6
##here we see tcell genes that change over time in response to alum, they all start up, then end up suppressed
    ##as above, we make one panel per gene##
    gene_i="ENSGACT00000009199"
    name="cxcl12a"
    cxcl12a=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
    
    gene_i="ENSGACT00000027215"
    name="HMGB1"
    HMGB1=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
    
    gene_i="ENSGACT00000019527"
    name="ANPEP"
    ANPEP=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
    
    gene_i="ENSGACT00000014797"
    name="dbnlb"
    dbnlb=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
    
    gene_i="ENSGACT00000017149"
    name="efhd2"
    efhd2=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
    
    gene_i="ENSGACT00000002459"
    name="gpr183a"
    gpr183a=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
    
    gene_i="ENSGACT00000009721"
    name="pik3cd"
    pik3cd=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
    
    gene_i="ENSGACT00000021067"
    name="ptger4b"
    ptger4b=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
    
    gene_i="ENSGACT00000003533"
    name="stat4"
    stat4=timecourseplot_nopop(gene_i,name,metadata.alum,readcount.alum, names.al, linescol =c(2,1))
    
    ##and then combine for our figure
    ggarrange(ANPEP, cxcl12a, dbnlb, efhd2, gpr183a, HMGB1, pik3cd, ptger4b, stat4,
              widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=3, nrow = 3,  common.legend=FALSE)
    ##only two of these (cxcl12a and HMGB1) overlap with fibrosis, so no need for additional panels##
  
##now we'll move on to cestode protein unique genes, starting with Fig 7 (nfkb responsive over time##)
    ##once again we use our functions to make one panel per gene.
    gene_i="ENSGACT00000016117"
    name="abhd12"
    abhd12=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp, names.cp, linescol =c(2,1))
    
    gene_i="ENSGACT00000013616"
    name="commd7"
    commd7=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp, names.cp, linescol =c(2,1))
    
    gene_i="ENSGACT00000019916"
    name="mul1a"
    mul1a=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp, names.cp, linescol =c(2,1))
    
    gene_i="ENSGACT00000017902"
    name="nfkbie"
    nfkbie=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp, names.cp, linescol =c(2,1))
    
    gene_i="ENSGACT00000009520"
    name="rab10"
    rab10=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp, names.cp, linescol =c(2,1))
    
    gene_i="ENSGACT00000015374"
    name="sash1a"
    sash1a=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp, names.cp, linescol =c(2,1))
    
    gene_i="ENSGACT00000003341"
    name="shrprbck1r"
    shrprbck1r=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp, names.cp, linescol =c(2,1))
    
    ##and then arrange them into one figure
    ggarrange(abhd12, commd7, mul1a, nfkbie, rab10, sash1a, shrprbck1r,
              widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=4, nrow = 2,  common.legend=FALSE)
    
##Next for Figure 8 we consider CP genes which are differentially responsive across populations, all inflammatory##
    ##once again we use our functions to make one panel per gene.
    gene_i="ENSGACT00000020499"
    name="apc"
    apc=timecourseplot(gene_i,name,metadata.cp,readcount.cp, names.cp, linescol =c(2,1))
    
    gene_i="ENSGACT00000013313"
    name="dusp10"
    dusp10=timecourseplot(gene_i,name,metadata.cp,readcount.cp, names.cp, linescol =c(2,1))
    
    gene_i="ENSGACT00000019348"
    name="grna"
    grna=timecourseplot(gene_i,name,metadata.cp,readcount.cp, names.cp, linescol =c(2,1))
    
    gene_i="ENSGACT00000013594"
    name="mef2cb"
    mef2cb=timecourseplot(gene_i,name,metadata.cp,readcount.cp, names.cp, linescol =c(2,1))
    
    gene_i="ENSGACT00000009525"
    name="tax1bp1a"
    tax1bp1a=timecourseplot(gene_i,name,metadata.cp,readcount.cp, names.cp, linescol =c(2,1))
    
    ##and combine into one figure##
    ggarrange(dusp10, grna, tax1bp1a, apc, mef2cb, 
              widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=3, nrow = 2,  common.legend=TRUE)

    
    
    
        
##Next up we'll make figure 9 which highlights spi1b, a putative fibrosis gene
    ##this one is a combination of two types of plots made individually
    gene_i="ENSGACT00000020522"
    name="spi1b"
    spi1b_cp=timecourseplot_nopop(gene_i,name,metadata.cp,readcount.cp, names.cp, linescol =c(2,1))
    spi1b_fib=fibplot(gene_i,name,metadata_fib,readcount_fib,names_fib, linescol =c(2,1))
    ##then combined into a multipanel figure##
    ggarrange(spi1b_cp,spi1b_fib, 
              widths = c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3), ncol=2, nrow = 1,  common.legend=FALSE)
    
    
##End with our trajectory analyses and plot for Figure 10###
    # PCA trajectory analysis
    
    # First choose which genes to focus on, picking all with significant LFCs
    # Get a list of file names
    Contrasts <- list.files("GLM_output")
    Nfiles <- length(Contrasts)
    
    sig.genes <- c()
    for(i in 1:Nfiles){
      file1 <- Contrasts[i]
      file1 <- fread(paste("GLM_Output/",file1, sep = ""))
      file1.sig <- file1[file1$padj < 0.01,] ##make sure to use correct cut-off
      dim(file1.sig)      
      sig.genes <- c(sig.genes, file1.sig$gene)
    }
    length(sig.genes)
    sig.genes <- unique(sig.genes)
    
    #Read in raw read data
    metadata2 <- read.csv("ExpDesign_noWA_no90.csv")
    dim(metadata2)
    metadata2 <- metadata2[c(1:276),]
    readcounts <- fread("reduced_matrix_HK.csv", data.table = F)
    dim(readcounts)
    readcounts <- readcounts[, names(readcounts) %in% sig.genes]
    dim(readcounts)
    
    ##run the PCA##
    PCA <- prcomp(readcounts, scale = T, center= T)
    
    ################DAPC just for PBS and  Alum###############
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
      summarytable$Treatment <- c(  rep("Alum", 9), rep("PBS", 9))
      summarytable$Day <- c(rep(c(1,10,42), 6))
      summarytable$Pop <- rep(c("GOS", "GOS", "GOS","RSL", "RSL", "RSL", "SAY", "SAY" ,"SAY"), 2)
      summarytable$colorplot <- as.numeric(factor(summarytable$Pop))
      summarytable$colorplot[summarytable$colorplot == 3] <- 4
    }
    
    
    
    # Make a plot with three panels, one per population
    {
      par(mfrow = c(1,3), mar = c(5,5,2,1))
      
      tempDFA <- DFA.p$x[metasubset$Population == "SAY",]
      plot(tempDFA[,2] ~ tempDFA[,1], pch = 1, cex = 0.5, col = rgb(0,0,0,1), xlim = c(-6.5,8), ylim = c(-6,6.5), xlab = "DFA axis 1", ylab = "DFA axis 2", cex.lab = 1.3, main = c("SAY"))
      summarytable2 <- summarytable[-c(17:18),] # remove Saline on later days
      summarytemp <- summarytable2[summarytable2$Pop == "SAY",]
      points(summarytemp$x1means, summarytemp$x2means, pch = c(16,16,16,1), col = "#C3A016FF", cex = 2.5, lwd = 3)
      arrows(x0 = summarytemp$x1means[c(4,1:2)], y0 = summarytemp$x2means[c(4,1:2)], x1 = summarytemp$x1means[c(1:3)], y1 = summarytemp$x2means[c(1:3)], lwd = 3, length = 0.2, col = "#C3A016FF")
      text(summarytemp$x1means+0.2, summarytemp$x2means, c( "1", "10","42", "PBS"),pos = 4)
      arrows(x0 = summarytemp$lowerCI1[c(1:4)], y0 = summarytemp$x2means[c(1:4)], x1 = summarytemp$upperCI1[c(1:4)], y1 = summarytemp$x2means[c(1:4)], code = 3, angle = 90, length  = 0.02, col = "grey50")
      arrows(x0 = summarytemp$x1means[c(1:4)], y0 = summarytemp$lowerCI2[c(1:4)], x1 = summarytemp$x1means[c(1:4)], y1 = summarytemp$upperCI2[c(1:4)], code = 3, angle = 90, length  = 0.02, col = "grey50")
      
      tempDFA <- DFA.p$x[metasubset$Population == "GOS",]
      plot(tempDFA[,2] ~ tempDFA[,1], pch = 1, cex = 0.5, col = rgb(0,0,0,1),  xlim = c(-6.5,8), ylim = c(-6,6.5), xlab = "DFA axis 1", ylab = "DFA axis 2", cex.lab = 1.3, main = c("GOS"))
      summarytable2 <- summarytable[-c(11:12),] # remove Saline on later days
      summarytemp <- summarytable2[summarytable2$Pop == "GOS",]
      points(summarytemp$x1means, summarytemp$x2means, pch = c(16,16,16, 1), col = "#58A787FF", cex = 2.5, lwd = 3)
      arrows(x0 = summarytemp$x1means[c(4,1:2)], y0 = summarytemp$x2means[c(4,1:2)], x1 = summarytemp$x1means[c(1:3)], y1 = summarytemp$x2means[c(1:3)], lwd = 3, length = 0.2, col = "#58A787FF")
      text(summarytemp$x1means+0.2, summarytemp$x2means, c( "1", "10","42", "PBS"),pos = 4)
      arrows(x0 = summarytemp$lowerCI1[c(1:4)], y0 = summarytemp$x2means[c(1:4)], x1 = summarytemp$upperCI1[c(1:4)], y1 = summarytemp$x2means[c(1:4)], code = 3, angle = 90, length  = 0.02, col = "grey50")
      arrows(x0 = summarytemp$x1means[c(1:4)], y0 = summarytemp$lowerCI2[c(1:4)], x1 = summarytemp$x1means[c(1:4)], y1 = summarytemp$upperCI2[c(1:4)], code = 3, angle = 90, length  = 0.02, col = "grey50")
      
      tempDFA <- DFA.p$x[metasubset$Population == "RSL",]
      plot(tempDFA[,2] ~ tempDFA[,1], pch = 1, cex = 0.5, col = rgb(0,0,0,1), xlim = c(-6.5,8), ylim = c(-6,6.5), xlab = "DFA axis 1", ylab = "DFA axis 2", cex.lab = 1.3, main = c("RSL"))
      summarytable2 <- summarytable[-c(14:15),] # remove Saline on later days
      summarytemp <- summarytable2[summarytable2$Pop == "RSL",]
      points(summarytemp$x1means, summarytemp$x2means, pch = c(16,16,16,1), col = "#0C1F4BFF", cex = 2.5, lwd = 3)
      arrows(x0 = summarytemp$x1means[c(4,1:2)], y0 = summarytemp$x2means[c(4,1:2)], x1 = summarytemp$x1means[c(1:3)], y1 = summarytemp$x2means[c(1:3)], lwd = 3, length = 0.2, col = "#0C1F4BFF")
      text(summarytemp$x1means+0.2, summarytemp$x2means, c( "1", "10","42", "PBS"),pos = 4)
      arrows(x0 = summarytemp$lowerCI1[c(1:4)], y0 = summarytemp$x2means[c(1:4)], x1 = summarytemp$upperCI1[c(1:4)], y1 = summarytemp$x2means[c(1:4)], code = 3, angle = 90, length  = 0.02, col = "grey50")
      arrows(x0 = summarytemp$x1means[c(1:4)], y0 = summarytemp$lowerCI2[c(1:4)], x1 = summarytemp$x1means[c(1:4)], y1 = summarytemp$upperCI2[c(1:4)], code = 3, angle = 90, length  = 0.02, col = "grey50")
    }
    
    
    
    
    
    ######## DAPC just for PBS and Worm protein################
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
      summarytable$Treatment <- c( rep("PBS", 9), rep("WP", 9))
      summarytable$Day <- c(rep(c(1,10,42), 6))
      summarytable$Pop <- rep(c("GOS", "GOS", "GOS","RSL", "RSL", "RSL", "SAY", "SAY", "SAY"), 2)
      summarytable$colorplot <- as.numeric(factor(summarytable$Pop))
      summarytable$colorplot[summarytable$colorplot == 3] <- 4
    }
    
    # Make a plot with three panels, one per population
    {  
      par(mfrow = c(1,3), mar = c(5,5,2,1))
      tempDFA <- DFA.p$x[metasubset$Population == "SAY",]
      plot(tempDFA[,2] ~ tempDFA[,1], pch = 1, cex = 0.5, col = rgb(0,0,0,1), xlim = c(-6, 9), ylim = c(-6, 6), xlab = "DFA axis 1", ylab = "DFA axis 2", cex.lab = 1.3, main = c("SAY"))
      summarytable2 <- summarytable[-c(8:9),] # remove Saline on later days
      summarytemp <- summarytable2[summarytable2$Pop == "SAY",]
      points(summarytemp$x1means, summarytemp$x2means, pch = c(1,16,16,16), col = "#C3A016FF", cex = 2.5, lwd = 3)
      arrows(x0 = summarytemp$x1means[c(1:3)], y0 = summarytemp$x2means[c(1:3)], x1 = summarytemp$x1means[c(2:4)], y1 = summarytemp$x2means[c(2:4)], lwd = 3, length = 0.2, col = "#C3A016FF")
      text(summarytemp$x1means+0.2, summarytemp$x2means, c("PBS", "1", "10","42"),pos = 4)
      arrows(x0 = summarytemp$lowerCI1[c(1:4)], y0 = summarytemp$x2means[c(1:4)], x1 = summarytemp$upperCI1[c(1:4)], y1 = summarytemp$x2means[c(1:4)], code = 3, angle = 90, length  = 0.02, col = "grey50")
      arrows(x0 = summarytemp$x1means[c(1:4)], y0 = summarytemp$lowerCI2[c(1:4)], x1 = summarytemp$x1means[c(1:4)], y1 = summarytemp$upperCI2[c(1:4)], code = 3, angle = 90, length  = 0.02, col = "grey50")
      
      tempDFA <- DFA.p$x[metasubset$Population == "GOS",]
      plot(tempDFA[,2] ~ tempDFA[,1], pch = 1, cex = 0.5, col = rgb(0,0,0,1), xlim = c(-6, 9), ylim = c(-6,6), xlab = "DFA axis 1", ylab = "DFA axis 2", cex.lab = 1.3, main = c("GOS"))
      summarytable2 <- summarytable[-c(2:3),] # remove PBS
      summarytemp <- summarytable2[summarytable2$Pop == "GOS",]
      points(summarytemp$x1means, summarytemp$x2means, pch = c(1,16,16,16), col = "#58A787FF", cex = 2.5, lwd = 3)
      arrows(x0 = summarytemp$x1means[c(1:3)], y0 = summarytemp$x2means[c(1:3)], x1 = summarytemp$x1means[c(2:4)], y1 = summarytemp$x2means[c(2:4)], lwd = 3, length = 0.2, col = "#58A787FF")
      text(summarytemp$x1means+0.2, summarytemp$x2means, c("PBS", "1", "10","42"),pos = 4)
      arrows(x0 = summarytemp$lowerCI1[c(1:4)], y0 = summarytemp$x2means[c(1:4)], x1 = summarytemp$upperCI1[c(1:4)], y1 = summarytemp$x2means[c(1:4)], code = 3, angle = 90, length  = 0.02, col = "grey50")
      arrows(x0 = summarytemp$x1means[c(1:4)], y0 = summarytemp$lowerCI2[c(1:4)], x1 = summarytemp$x1means[c(1:4)], y1 = summarytemp$upperCI2[c(1:4)], code = 3, angle = 90, length  = 0.02, col = "grey50")
      
      tempDFA <- DFA.p$x[metasubset$Population == "RSL",]
      plot(tempDFA[,2] ~ tempDFA[,1], pch = 1, cex = 0.5, col = rgb(0,0,0,1), xlim = c(-6, 9), ylim = c(-6, 6), xlab = "DFA axis 1", ylab = "DFA axis 2", cex.lab = 1.3, main = c("RSL"))
      summarytable2 <- summarytable[-c(5:6),] # remove Saline on later days
      summarytemp <- summarytable2[summarytable2$Pop == "RSL",]
      points(summarytemp$x1means, summarytemp$x2means, pch = c(1,16,16,16), col = "#0C1F4BFF", cex = 2.5, lwd = 3)
      text(summarytemp$x1means+0.2, summarytemp$x2means, c("PBS", "1", "10","42"),pos = 4)
      arrows(x0 = summarytemp$x1means[c(1:3)], y0 = summarytemp$x2means[c(1:3)], x1 = summarytemp$x1means[c(2:4)], y1 = summarytemp$x2means[c(2:4)], lwd = 3, length = 0.2, col = "#0C1F4BFF")
      arrows(x0 = summarytemp$lowerCI1[c(1:4)], y0 = summarytemp$x2means[c(1:4)], x1 = summarytemp$upperCI1[c(1:4)], y1 = summarytemp$x2means[c(1:4)], code = 3, angle = 90, length  = 0.02, col = "grey50")
      arrows(x0 = summarytemp$x1means[c(1:4)], y0 = summarytemp$lowerCI2[c(1:4)], x1 = summarytemp$x1means[c(1:4)], y1 = summarytemp$upperCI2[c(1:4)], code = 3, angle = 90, length  = 0.02, col = "grey50")}
    
