

Run_Caller1v2_Intersect <- function(SubID, CNV_Caller1_data_table, CNV_Caller2_data_table, overlap_threshold=0.5, merged_conf_threshold=15, merged_numsnp_threshold=10, merged_length_threshold=20000) {
  CNV_Caller1_ID_tmp <- CNV_Caller1_data_table[CNV_Caller1_data_table$ID == SubID]
  CNV_Caller2_ID_tmp <- CNV_Caller2_data_table[CNV_Caller2_data_table$ID == SubID]
  Caller1vs2.byID.Overlaps <- findOverlaps(CNV_Caller1_ID_tmp, CNV_Caller2_ID_tmp,
                                           maxgap=-1L, minoverlap=1L,
                                           type=c('any', 'start', 'end', 'within', 'equal'),
                                           select= 'all', 
                                           ignore.strand=FALSE)
  Caller1vs2.byID.Intersect <- pintersect(CNV_Caller1_ID_tmp[queryHits(Caller1vs2.byID.Overlaps)], CNV_Caller2_ID_tmp[subjectHits(Caller1vs2.byID.Overlaps)])
  #print(Caller1vs2.byID.Intersect)
  
  if(isEmpty(Caller1vs2.byID.Intersect)){
    pb$tick()
    return(GRanges(seqnames=NULL,ranges=NULL,strand=NULL, ID=SubID))
  } 
  else {
    CNV_Caller1.byID.AnyOverlap <- CNV_Caller1_ID_tmp[queryHits(Caller1vs2.byID.Overlaps)]
    CNV_Caller1.byID.AnyOverlap$PercentOverlap.OrigCall <- width(Caller1vs2.byID.Intersect) / width(CNV_Caller1.byID.AnyOverlap)
    CNV_Caller2.byID.AnyOverlap <- CNV_Caller2_ID_tmp[subjectHits(Caller1vs2.byID.Overlaps)]
    CNV_Caller2.byID.AnyOverlap$PercentOverlap.OrigCall <- width(Caller1vs2.byID.Intersect) / width(CNV_Caller2.byID.AnyOverlap)
    Caller1vs2.Overlap.byID.Append <- sort(sortSeqlevels(c(CNV_Caller1.byID.AnyOverlap, CNV_Caller2.byID.AnyOverlap)))
    Caller1vs2.Overlap.byID.Overlap.Filtered <- Caller1vs2.Overlap.byID.Append[Caller1vs2.Overlap.byID.Append$PercentOverlap.OrigCall>=overlap_threshold,]
    Caller1vs2.Overlap.byID.Overlap.Conf.SNPNUM.Length.Filtered <- Caller1vs2.Overlap.byID.Overlap.Filtered[Caller1vs2.Overlap.byID.Overlap.Filtered$numsnp >= merged_numsnp_threshold & Caller1vs2.Overlap.byID.Overlap.Filtered$length >= merged_length_threshold & Caller1vs2.Overlap.byID.Overlap.Filtered$confidencescore >= merged_conf_threshold,]
    pb$tick()
    return(Caller1vs2.Overlap.byID.Overlap.Conf.SNPNUM.Length.Filtered)
    
  }
} 


Find_overlap_with_referenceSets<-function(CNVCall.GR, referenceData.GR,overlap_threshold = 0.5){
  
  CNVCall.GR.split <- split(CNVCall.GR, seq_along(CNVCall.GR))
  CNVCall.NoHits <- list()
  CNVCall.Hits <- list()
  
  pb <- progress_bar$new(total = length(CNVCall.GR.split), format = "[:bar] :percent")
  for(name in names(CNVCall.GR.split)){
    # print(name)
    
    CNVCall.tmp <- CNVCall.GR.split[[name]]
    CNV.MinSeg.Overlap <-(CNVCall.tmp$length * overlap_threshold)
    # CNV.MinSeg.Overlap <- CNVCall.tmp$HalfCNV
    overlapstodrop.byCNV <- findOverlaps(CNVCall.tmp, referenceData.GR, minoverlap = CNV.MinSeg.Overlap, ignore.strand=TRUE)
    CNVstoDrop.hits <- S4Vectors::queryHits(overlapstodrop.byCNV)
    if(isEmpty(CNVstoDrop.hits)){
      CNVCall.NoHits <- c(CNVCall.NoHits, CNVCall.tmp)
    } else {  
      
      CNVCall.Hits <- c(CNVCall.Hits, CNVCall.tmp)
    }
    pb$tick()
  }  
  
  # print(CNVCall.Hits)
  CNVCall.NoHits.df <- lapply(CNVCall.NoHits, as.data.frame)
  CNVCall.NoHits.df.complete <- bind_rows(CNVCall.NoHits.df)
  
  CNVCall.NoHits.df.complete.distinct <- distinct(CNVCall.NoHits.df.complete)
  
  CNVCall.Hits.df <- lapply(CNVCall.Hits, as.data.frame)
  CNVCall.Hits.df.complete <- bind_rows(CNVCall.Hits.df)
  CNVCall.Hits.df.complete.distinct <- distinct(CNVCall.Hits.df.complete)
  
  return(list(CNVCall.NoHits.df.complete.distinct, CNVCall.Hits.df.complete.distinct))
}


Find_overlap_with_KnownRiskLoci<-function(referenceData.GR, CNVCall.GR){

  referenceData.GR.split <- split(referenceData.GR, seq_along(referenceData.GR))
  KnownRiskLoci.Hits <- list()
  
  pb <- progress_bar$new(total = length(referenceData.GR.split), format = "[:bar] :percent")
  for(name in names(referenceData.GR.split)){
    
    referenceData.GR.tmp <- referenceData.GR.split[[name]]
    CNV.MinSeg.Overlap <-floor(referenceData.GR.tmp$FortyPercLoci)
    overlapstodrop.byCNV <- suppressWarnings(findOverlaps(CNVCall.GR, referenceData.GR.tmp, minoverlap = CNV.MinSeg.Overlap, ignore.strand=TRUE))
    CNVstoDrop.hits <- S4Vectors::queryHits(overlapstodrop.byCNV)

    if(!isEmpty(CNVstoDrop.hits)){
      KnowRiskLoci.Overlap <- CNVCall.GR[CNVstoDrop.hits]
      KnowRiskLoci.Overlap$Syndrome <- referenceData.GR.tmp$Syndrome
      KnowRiskLoci.Overlap$SyndromeLocation <- referenceData.GR.tmp$FullLocation_hg38
      KnowRiskLoci.Overlap$DiseaseAssoc <- referenceData.GR.tmp$DiseaseAssociation_SCZMarshall_ASDSanders_DDD_NDDCoe
      KnowRiskLoci.Overlap$Ref<-referenceData.GR.tmp$Ref
      KnownRiskLoci.Hits <- c(KnownRiskLoci.Hits, KnowRiskLoci.Overlap)
      # CNVCall.CentTeloSegDupHits <- c(CNVCall.CentTeloSegDupHits, CNVCall.tmp)
    }
    pb$tick()
  }
  KnownRiskLoci.Hits.df <- lapply(KnownRiskLoci.Hits, as.data.frame)
  KnownRiskLoci.Hits.df.complete <- bind_rows(KnownRiskLoci.Hits.df)

  return(KnownRiskLoci.Hits.df.complete)
}

# with NRXN1 Exons
Find_overlap_with_NRXN1_Exons<-function(referenceData.GR, CNVCall.GR){
  
  referenceData.GR.split <- split(referenceData.GR, seq_along(referenceData.GR))
  KnownRiskLoci.Hits <- list()
  
  pb <- progress_bar$new(total = length(referenceData.GR.split), format = "[:bar] :percent")
  for(name in names(referenceData.GR.split)){
  
    referenceData.GR.tmp <- referenceData.GR.split[[name]]
    CNV.MinSeg.Overlap <-floor(referenceData.GR.tmp$length)
    overlapstodrop.byCNV <- suppressWarnings(findOverlaps(CNVCall.GR, referenceData.GR.tmp, minoverlap = CNV.MinSeg.Overlap, ignore.strand=TRUE))
    CNVstoDrop.hits <- S4Vectors::queryHits(overlapstodrop.byCNV)
    
    if(!isEmpty(CNVstoDrop.hits)){
      KnowRiskLoci.Overlap <- CNVCall.GR[CNVstoDrop.hits]
      KnowRiskLoci.Overlap$Syndrome <- 'NRXN1'
      KnowRiskLoci.Overlap$DiseaseAssoc <- 'SCZ'
      KnownRiskLoci.Hits <- c(KnownRiskLoci.Hits, KnowRiskLoci.Overlap)
    }
    pb$tick()
  }
  KnownRiskLoci.Hits.df <- lapply(KnownRiskLoci.Hits, as.data.frame)
  KnownRiskLoci.Hits.df.complete <- bind_rows(KnownRiskLoci.Hits.df)
  
  return(KnownRiskLoci.Hits.df.complete)
}







Split_GRObjects_bySubject<-function(CNVCalls.GR){
  CNVCalls.GR.bySub <- split(CNVCalls.GR, CNVCalls.GR$ID)
  
  CNVCalls.Reduced.AllSub <- list()
  pb <- progress_bar$new(total = length(CNVCalls.GR.bySub), format = "[:bar] :percent")
  #reduce genomic ranges to unique non-overlapping ranges
  for(name in names(CNVCalls.GR.bySub)){
    # print(name)
    CNVCalls.Sub.GR <- CNVCalls.GR.bySub[[name]]
    # print(head(CNVCalls.Sub.GR))
    CNVCalls.Sub.GR.reduced <- GenomicRanges::reduce(CNVCalls.Sub.GR)
    CNVCalls.Sub.GR.reduced$ID <- name
    # print(CNVCalls.Sub.GR.reduced)
    CNVCalls.Reduced.AllSub <- c(CNVCalls.Reduced.AllSub, CNVCalls.Sub.GR.reduced)
    pb$tick()
  }
  
  CNVCalls.Reduced.AllSub.df <- lapply(CNVCalls.Reduced.AllSub, as.data.frame)
  CNVCalls.Reduced.AllSub.df.complete <- bind_rows(CNVCalls.Reduced.AllSub.df)
  CNVCalls.Reduced.AllSub.df.complete$ID.Short <- CNVCalls.Reduced.AllSub$ID
  CNVCalls.Reduced.AllSub.df.complete$HalfWidth <- floor(CNVCalls.Reduced.AllSub.df.complete$width/2)
  
  return(CNVCalls.Reduced.AllSub.df.complete)
  
}



Find_overlap_with_refSeq<-function(CNVCall.GR, referenceData.GR){
  CNVCall.GR.AllSubj <- split(CNVCall.GR, seq_along(CNVCall.GR))
  RareCNV <- data.frame()
  pb <- progress_bar$new(total = length(CNVCall.GR.AllSubj), format = "[:bar] :percent")
  for(name in names(CNVCall.GR.AllSubj)){
    
    CNVCall.GR.oneSubject <- CNVCall.GR.AllSubj[[name]]
    Geneoverlaps.byCNV <- findOverlaps(referenceData.GR, CNVCall.GR.oneSubject, ignore.strand=TRUE)
    CNV.genehits <- S4Vectors::queryHits(Geneoverlaps.byCNV)
    CNV.cnvgenes <- referenceData.GR[CNV.genehits]
    CNV.cnvgenes.df <- CNV.cnvgenes %>% as.data.frame()
    if(!isEmpty(CNV.cnvgenes.df)){
      # print('Gene Hits')
      CNV.cnvgenes.df$ID <- unique(CNVCall.GR.oneSubject$ID)
      CNV.cnvgenes.df$FullLocation_hg38 <- CNVCall.GR.oneSubject$FullLocation_hg38
      CNV.cnvgenes.df$FullLocation_numsnp <- CNVCall.GR.oneSubject$numsnp
      CNV.cnvgenes.df$Algorithm <- CNVCall.GR.oneSubject$Algorithm
      CNV.cnvgenes.df$confidencescore <- CNVCall.GR.oneSubject$confidencescore
      RareCNV <- rbind(RareCNV, CNV.cnvgenes.df)
    }
    pb$tick()

  }
  return(RareCNV)
}


Run.LOEUFSum_byModule <- function(Modules, LOEUFScore_df, SubjectIDs, DelDup_Status) {
  All.CNVSubset.LOEUFScore.Sum.byModule.bySubj <- data.frame(ID = SubjectIDs)
  for(Module in Modules){
    # print(Module)
    LOEUFScore_tmp = LOEUFScore_df[LOEUFScore_df$ModuleNumName == Module,]
    LOEUFScore_tmp <- LOEUFScore_tmp[!is.na(LOEUFScore_tmp$LOEUF.perc),]
    # print(head(LOEUFScore_tmp))
    LOEUFScore.Sum.bySubj <- group_by(LOEUFScore_tmp, ID) %>% summarise(sum(LOEUF.perc)) %>% as.data.frame()
    LOEUFScore.col.name <- paste(Module, DelDup_Status,'LOEUFScoreSum_min10SNPs20kbCNVs', sep = '_')
    colnames(LOEUFScore.Sum.bySubj) <- c('ID', LOEUFScore.col.name)
    All.CNVSubset.LOEUFScore.Sum.byModule.bySubj <-left_join(All.CNVSubset.LOEUFScore.Sum.byModule.bySubj, LOEUFScore.Sum.bySubj)
  }
  return(All.CNVSubset.LOEUFScore.Sum.byModule.bySubj)
  
}

Make_Summary_Stats_for_RareCNV_bySubj<-function(RareCNV.wBrainSpan, DelDup_Status, Batch12.PennCNV.QC.GSAMD.Outlier.distinct){
  # Function 
  RareCNV.LOEUF.top10perc <- dplyr::filter(RareCNV.wBrainSpan, LOEUF.top10perc == 'LOEUF.top10perc')
  # Count for each subject, how many genes that are top10 perc LOEUF are affected by CNVs
  RareCNV.LOEUF.top10perc.GeneCount <- group_by(RareCNV.LOEUF.top10perc, ID) %>% summarise(length(unique(HGNCUpdatedSymbol)))
  colnames(RareCNV.LOEUF.top10perc.GeneCount)[2] <- paste0('Total_FullExon',DelDup_Status,'_LOEUF.top10perc.GeneCount_min10SNPs20kbCNVs') 
  
  # Count for each subject, how many genes that are top20 perc LOEUF are affected by CNVs
  RareCNV.LOEUF.top20perc <- dplyr::filter(RareCNV.wBrainSpan, LOEUF.top20perc == 'LOEUF.top20perc')
  RareCNV.LOEUF.top20perc.GeneCount <- group_by(RareCNV.LOEUF.top20perc, ID) %>% summarise(length(unique(HGNCUpdatedSymbol)))
  colnames(RareCNV.LOEUF.top20perc.GeneCount)[2] <-  paste0('Total_FullExon',DelDup_Status,'_LOEUF.top20perc.GeneCount_min10SNPs20kbCNVs')
  
  # Count for each subject, how many genes in total are affected by CNVs
  RareCNV.GeneCount <- group_by(RareCNV.wBrainSpan, ID) %>% summarise(length(unique(HGNCUpdatedSymbol)))
  colnames(RareCNV.GeneCount)[2] <- paste0('Total_FullExon',DelDup_Status,'_GeneCount_min10SNPs20kbCNVs')
  
  
  RareCNV.CountPerModule <- table(RareCNV.wBrainSpan$ModuleNumName, RareCNV.wBrainSpan$ID) %>% as.data.frame()# summarise(length(unique(HGNCUpdatedSymbol)))
  RareCNV.CountPerModule$Var1 <- str_c(RareCNV.CountPerModule$Var1, paste0(DelDup_Status,"_GeneCount_min10SNPs20kbCNVs"), sep = '_')
  RareCNV.CountPerModule.wide <- spread(RareCNV.CountPerModule, Var1, Freq)
  colnames(RareCNV.CountPerModule.wide)[1] <- 'ID'
  
  RareCNV.LOEUF.top10perc.CountPerModule <- table(RareCNV.LOEUF.top10perc$ModuleNumName, RareCNV.LOEUF.top10perc$ID) %>% as.data.frame()# summarise(length(unique(HGNCUpdatedSymbol)))
  RareCNV.LOEUF.top10perc.CountPerModule$Var1 <- str_c(RareCNV.LOEUF.top10perc.CountPerModule$Var1, paste0(DelDup_Status,"_LOEUF.top10perc.GeneCount_min10SNPs20kbCNVs") , sep = '_')
  RareCNV.LOEUF.top10perc.CountPerModule.wide <- spread(RareCNV.LOEUF.top10perc.CountPerModule, Var1, Freq)
  colnames(RareCNV.LOEUF.top10perc.CountPerModule.wide)[1] <- 'ID'

  RareCNV.LOEUF.top20perc.CountPerModule <- table(RareCNV.LOEUF.top20perc$ModuleNumName, RareCNV.LOEUF.top20perc$ID) %>% as.data.frame()# summarise(length(unique(HGNCUpdatedSymbol)))
  RareCNV.LOEUF.top20perc.CountPerModule$Var1 <- str_c(RareCNV.LOEUF.top20perc.CountPerModule$Var1, paste0(DelDup_Status,"_LOEUF.top20perc.GeneCount_min10SNPs20kbCNVs") , sep = '_')
  RareCNV.LOEUF.top20perc.CountPerModule.wide <- spread(RareCNV.LOEUF.top20perc.CountPerModule, Var1, Freq)
  colnames(RareCNV.LOEUF.top20perc.CountPerModule.wide)[1] <- 'ID'

  All.CNVSubset.LOEUFScore.Sum.byModule.bySubj<-Run.LOEUFSum_byModule(unique(Kangmodules.minmod30.DS2.ModuleNumNames), RareCNV.wBrainSpan, unique(Batch12.PennCNV.QC.GSAMD.Outlier.distinct$ID), DelDup_Status)
  
  RareCNV.wLOEUF <- RareCNV.wBrainSpan[!is.na(RareCNV.wBrainSpan$LOEUF.perc),]
  RareCNV.LOEUFScore.Sum.bySubj <- group_by(RareCNV.wLOEUF, as.factor(ID)) %>% summarise(sum(LOEUF.perc)) %>% as.data.frame()
  colnames(RareCNV.LOEUFScore.Sum.bySubj) <- c('ID', paste0("Total_FullExon",DelDup_Status,"_LOEUFScoreSum_min10SNPs20kbCNVs"))
  
  All.CNVSubset.LOEUFScore.Sum.byModule.bySubj <- left_join(All.CNVSubset.LOEUFScore.Sum.byModule.bySubj, RareCNV.LOEUFScore.Sum.bySubj)
  All.CNVSubset.LOEUFScore.Sum.byModule.bySubj[is.na(All.CNVSubset.LOEUFScore.Sum.byModule.bySubj)] <- 0 
  
  AllLargeCNV_bySubj.Summarystats <- (left_join(Batch12.PennCNV.QC.GSAMD.Outlier.distinct, RareCNV.GeneCount) %>% left_join(RareCNV.CountPerModule.wide) 
                                                                    %>% left_join(RareCNV.LOEUF.top10perc.CountPerModule.wide) %>% left_join(RareCNV.LOEUF.top20perc.CountPerModule.wide) %>% left_join(RareCNV.LOEUFScore.Sum.bySubj) 
                                                                    %>% left_join(All.CNVSubset.LOEUFScore.Sum.byModule.bySubj)) 
  
  
  
  
  AllLargeCNV_bySubj.Summarystats$ASD_Gene <- 0
  AllLargeCNV_bySubj.Summarystats$ASD_Gene[(AllLargeCNV_bySubj.Summarystats$ID %in% RareCNV.wBrainSpan$ID[RareCNV.wBrainSpan$ASD == 'ASD'])] <- 1
  
  AllLargeCNV_bySubj.Summarystats$NDD_Gene <- 0
  AllLargeCNV_bySubj.Summarystats$NDD_Gene[(AllLargeCNV_bySubj.Summarystats$ID %in% RareCNV.wBrainSpan$ID[RareCNV.wBrainSpan$NDD == 'NDD'])] <- 1
  
  NDD.Gene.IDs <- RareCNV.wBrainSpan[!is.na(RareCNV.wBrainSpan$NDD),]
  NDD.Gene.IDs.unique <- unique(NDD.Gene.IDs$ID) %>% as.data.frame()
  
  AllLargeCNV_bySubj.Summarystats$SCZ_Gene<- 0
  AllLargeCNV_bySubj.Summarystats$SCZ_Gene[(AllLargeCNV_bySubj.Summarystats$ID %in% RareCNV.wBrainSpan$ID[RareCNV.wBrainSpan$SCZ == 'SCZ'])] <- 1
  
  col_lengths<-length(colnames(AllLargeCNV_bySubj.Summarystats))
  
  colnames(AllLargeCNV_bySubj.Summarystats)[col_lengths-2]<-paste0("ASD_Gene_FullExon",DelDup_Status,"_min10SNPs20kbCNVs")
  colnames(AllLargeCNV_bySubj.Summarystats)[col_lengths-1]<-paste0("NDD_Gene_FullExon",DelDup_Status,"_min10SNPs20kbCNVs")
  colnames(AllLargeCNV_bySubj.Summarystats)[col_lengths]<-paste0("SCZ_Gene_FullExon",DelDup_Status,"_min10SNPs20kbCNVs")
  
  return(AllLargeCNV_bySubj.Summarystats)
  
}






