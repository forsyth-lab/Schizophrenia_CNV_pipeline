# Header

cat("
 _____ ______      _   _______   _       ___  ______
|  __ \\| ___ \\    | \\ | |  _  \\ | |     / _ \\ | ___ \\
| |  \\/| |_/ /__ _|  \\| | | | | | |    / /_\\ \\| |_/ /
| | __ |    // _` | . ` | | | | | |    |  _  || ___ \\
| |_\\ \\| |\\ \\ (_| | |\\  | |/ /  | |____| | | || |_/ /
 \\____/\\_| \\_\\__,_\\_| \\_/___/    \\_____\\_| |_/\\____/

")

options(warn = -1)

cat("\nCNV Preprocessing Pipline, developed by Dr.Jennifer Forsyth, Jinhan Zhu\n")

# load libraries:
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Cairo))

#Reading and writing tables
suppressPackageStartupMessages(library(openxlsx))

#String operations
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tools))

#Functional programming
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(broom))

#analysis
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
# UI
suppressPackageStartupMessages(library(progress))

source("GRange_Overlap_Functions.R")
#Read in PennCNV Calls for Batch 1 in hg38 coords
Batch1.PennCNV.Min1kb3SNPConf10.hg38.Coords.AllAvail <- read.xlsx('../input_files/Batch1.PennCNV.xlsx')


#Read in QuantiSNP Calls for Batch 1 in hg38 coords
Batch1.QuantiSNP.Min1kb3SNPConf10.hg38.Coords.AllAvail <- read.xlsx('../input_files/Batch1.QuantiSNP.xlsx')

#Read in QuantiSNP Calls for Batch 2
Batch2.QuantiSNP.Min1kb3SNPConf10<-fread("../input_files/Batch2.QuantiSNP.txt", header = TRUE)
#Read in PennCNV Calls for Batch 2
Batch2.PennCNV.Min1kb3SNPConf10<-fread("../input_files/Batch2.PennCNV.txt", header = TRUE)

cat("\n################################################### Find Overlaps Between CNV Call Algorithms ##############################################\n")
Batch2.Penn.GR.Min1kb3SNPConf10 <- makeGRangesFromDataFrame(Batch2.PennCNV.Min1kb3SNPConf10,
                                                            keep.extra.columns=TRUE,
                                                            ignore.strand=FALSE,
                                                            seqinfo=NULL,
                                                            seqnames.field=c('seqnames', 'seqname',
                                                                             'chromosome', 'chrom',
                                                                             'Chr', 'chromosome_name',
                                                                             'seqid'),
                                                            start.field='start',
                                                            end.field=c('end', 'stop'),
                                                            strand.field='strand',
                                                            starts.in.df.are.0based=FALSE)
Batch2.QSNP.GR.Min1kb3SNPConf10 <- makeGRangesFromDataFrame(Batch2.QuantiSNP.Min1kb3SNPConf10,
                                                            keep.extra.columns=TRUE,
                                                            ignore.strand=FALSE,
                                                            seqinfo=NULL,
                                                            seqnames.field=c('seqnames', 'seqname',
                                                                             'chromosome', 'chrom',
                                                                             'Chr', 'chromosome_name',
                                                                             'seqid'),
                                                            start.field='start',
                                                            end.field=c('end', 'stop'),
                                                            strand.field='strand',
                                                            starts.in.df.are.0based=FALSE)


Batch1.Penn.GR.hg38.Min1kb3SNPConf10 <- makeGRangesFromDataFrame(Batch1.PennCNV.Min1kb3SNPConf10.hg38.Coords.AllAvail,
                                                                 keep.extra.columns=TRUE,
                                                                 ignore.strand=FALSE,
                                                                 seqinfo=NULL,
                                                                 seqnames.field=c('seqnames', 'seqname',
                                                                                  'chromosome', 'chrom',
                                                                                  'Chr', 'chromosome_name', 'chr',
                                                                                  'seqid'),
                                                                 start.field='start',
                                                                 end.field=c('end', 'stop'),
                                                                 strand.field='strand',
                                                                 starts.in.df.are.0based=FALSE)

Batch1.QSNP.GR.hg38.Min1kb3SNPConf10 <- makeGRangesFromDataFrame(Batch1.QuantiSNP.Min1kb3SNPConf10.hg38.Coords.AllAvail,
                                                                 keep.extra.columns=TRUE,
                                                                 ignore.strand=FALSE,
                                                                 seqinfo=NULL,
                                                                 seqnames.field=c('seqnames', 'seqname',
                                                                                  'chromosome', 'chrom',
                                                                                  'Chr', 'chromosome_name',
                                                                                  'seqid'),
                                                                 start.field='start',
                                                                 end.field=c('end', 'stop'),
                                                                 strand.field='strand',
                                                                 starts.in.df.are.0based=FALSE)

####################################
#Subset Batch 2 Penn GR Object for Dels vs Dup
Batch2.Penn.GR.Min1kb3SNPConf10.Del <- Batch2.Penn.GR.Min1kb3SNPConf10[Batch2.Penn.GR.Min1kb3SNPConf10$LossvGain == 'Loss']
Batch2.Penn.GR.Min1kb3SNPConf10.Dup <- Batch2.Penn.GR.Min1kb3SNPConf10[Batch2.Penn.GR.Min1kb3SNPConf10$LossvGain == 'Gain']

#Subset Batch 2 QuantiSNP GR Object for Dels vs Dup
Batch2.QSNP.GR.Min1kb3SNPConf10.Del <- Batch2.QSNP.GR.Min1kb3SNPConf10[Batch2.QSNP.GR.Min1kb3SNPConf10$LossvGain == 'Loss']
Batch2.QSNP.GR.Min1kb3SNPConf10.Dup <- Batch2.QSNP.GR.Min1kb3SNPConf10[Batch2.QSNP.GR.Min1kb3SNPConf10$LossvGain == 'Gain']

#Get Batch 2 IDs with Del Calls in both algorithms
Batch2.CNV.CrossAlgorithm.Del.IDs <- Batch2.Penn.GR.Min1kb3SNPConf10.Del$ID[Batch2.Penn.GR.Min1kb3SNPConf10.Del$ID %in% Batch2.QSNP.GR.Min1kb3SNPConf10.Del$ID]
Batch2.CNV.CrossAlgorithm.Del.UniqueIDs <- Batch2.CNV.CrossAlgorithm.Del.IDs[!duplicated(Batch2.CNV.CrossAlgorithm.Del.IDs)]

#Get Batch 2 IDs with Dup Calls in both algorithms
Batch2.CNV.CrossAlgorithm.Dup.IDs <- Batch2.Penn.GR.Min1kb3SNPConf10.Dup$ID[Batch2.Penn.GR.Min1kb3SNPConf10.Dup$ID %in% Batch2.QSNP.GR.Min1kb3SNPConf10.Dup$ID]
Batch2.CNV.CrossAlgorithm.Dup.UniqueIDs <- Batch2.CNV.CrossAlgorithm.Dup.IDs[!duplicated(Batch2.CNV.CrossAlgorithm.Dup.IDs)]


#Subset Batch 1 hg38 Penn GR Object for Dels vs Dup
Batch1.Penn.GR.hg38.Min1kb3SNPConf10.Del <- Batch1.Penn.GR.hg38.Min1kb3SNPConf10[Batch1.Penn.GR.hg38.Min1kb3SNPConf10$LossvGain == 'Loss']
Batch1.Penn.GR.hg38.Min1kb3SNPConf10.Dup <- Batch1.Penn.GR.hg38.Min1kb3SNPConf10[Batch1.Penn.GR.hg38.Min1kb3SNPConf10$LossvGain == 'Gain']

#Subset Batch 1 hg38 QuantiSNP GR Object for Dels vs Dup
Batch1.QSNP.GR.hg38.Min1kb3SNPConf10.Del <- Batch1.QSNP.GR.hg38.Min1kb3SNPConf10[Batch1.QSNP.GR.hg38.Min1kb3SNPConf10$LossvGain == 'Loss']
Batch1.QSNP.GR.hg38.Min1kb3SNPConf10.Dup <- Batch1.QSNP.GR.hg38.Min1kb3SNPConf10[Batch1.QSNP.GR.hg38.Min1kb3SNPConf10$LossvGain == 'Gain']



#Get Batch 1 IDs with Del Calls in both algorithms
Batch1.CNV.hg38.CrossAlgorithm.Del.IDs <- Batch1.Penn.GR.hg38.Min1kb3SNPConf10.Del$ID[Batch1.Penn.GR.hg38.Min1kb3SNPConf10.Del$ID %in% Batch1.QSNP.GR.hg38.Min1kb3SNPConf10.Del$ID]
Batch1.CNV.hg38.CrossAlgorithm.Del.UniqueIDs <- Batch1.CNV.hg38.CrossAlgorithm.Del.IDs[!duplicated(Batch1.CNV.hg38.CrossAlgorithm.Del.IDs)]

#Get Batch 1 IDs with Dup Calls in both algorithms
Batch1.CNV.hg38.CrossAlgorithm.Dup.IDs <- Batch1.Penn.GR.hg38.Min1kb3SNPConf10.Dup$ID[Batch1.Penn.GR.hg38.Min1kb3SNPConf10.Dup$ID %in% Batch1.QSNP.GR.hg38.Min1kb3SNPConf10.Dup$ID]
Batch1.CNV.hg38.CrossAlgorithm.Dup.UniqueIDs <- Batch1.CNV.hg38.CrossAlgorithm.Dup.IDs[!duplicated(Batch1.CNV.hg38.CrossAlgorithm.Dup.IDs)]


#All Batch2 Del Call Subs
cat("Run for Batch 2 Deletion","\n")
pb <- progress_bar$new(total = length(Batch2.CNV.CrossAlgorithm.Del.UniqueIDs), format = "[:bar] :percent")
Batch2.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All <- map(Batch2.CNV.CrossAlgorithm.Del.UniqueIDs, Run_Caller1v2_Intersect, Batch2.Penn.GR.Min1kb3SNPConf10.Del, Batch2.QSNP.GR.Min1kb3SNPConf10.Del) #%>% reduce(cbind) %>% as.data.frame()
names(Batch2.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All) <- Batch2.CNV.CrossAlgorithm.Del.UniqueIDs

#All Batch1 Del Call Subs
cat("Run for Batch 1 Deletion","\n")
pb <- progress_bar$new(total = length(Batch1.CNV.hg38.CrossAlgorithm.Del.UniqueIDs), format = "[:bar] :percent")
Batch1.hg38.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All <- map(Batch1.CNV.hg38.CrossAlgorithm.Del.UniqueIDs, Run_Caller1v2_Intersect, Batch1.Penn.GR.hg38.Min1kb3SNPConf10.Del, Batch1.QSNP.GR.hg38.Min1kb3SNPConf10.Del) #%>% reduce(cbind) %>% as.data.frame()
names(Batch1.hg38.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All) <- Batch1.CNV.hg38.CrossAlgorithm.Del.UniqueIDs


#All Batch2 Dup Call Subs
cat("Run for Batch 2 Duplication","\n")
pb <- progress_bar$new(total = length(Batch2.CNV.CrossAlgorithm.Dup.UniqueIDs), format = "[:bar] :percent")
Batch2.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All <- map(Batch2.CNV.CrossAlgorithm.Dup.UniqueIDs, Run_Caller1v2_Intersect, Batch2.Penn.GR.Min1kb3SNPConf10.Dup, Batch2.QSNP.GR.Min1kb3SNPConf10.Dup) #%>% reduce(cbind) %>% as.data.frame()
names(Batch2.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All) <- Batch2.CNV.CrossAlgorithm.Dup.UniqueIDs


#All Batch1 Dup Call Subs
cat("Run for Batch 1 Duplication","\n")
pb <- progress_bar$new(total = length(Batch1.CNV.hg38.CrossAlgorithm.Dup.UniqueIDs), format = "[:bar] :percent")
Batch1.hg38.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All <- map(Batch1.CNV.hg38.CrossAlgorithm.Dup.UniqueIDs, Run_Caller1v2_Intersect, Batch1.Penn.GR.hg38.Min1kb3SNPConf10.Dup, Batch1.QSNP.GR.hg38.Min1kb3SNPConf10.Dup) #%>% reduce(cbind) %>% as.data.frame()
names(Batch1.hg38.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All) <- Batch1.CNV.hg38.CrossAlgorithm.Dup.UniqueIDs

cat("\n########################### Drop CNVs overlapping with Centromere, Telomere, Seg Dup Reference Genome Datasets #############################\n")
SegDup.min10kb.Cento.Telo.hg38 <- read.xlsx('../GeneticReferenceData/SegDup.min10kb.Cento.Telo.hg38.forCNVAnalysis.042223.xlsx')
SegDup.min10kb.Cento.Telo.hg38$chrom <- str_replace_all(SegDup.min10kb.Cento.Telo.hg38$chrom, 'chr', '')
#Create GR Objects for Centromeres, Telomeres, Segmental Duplication regions to be dropped 
Cent.Telo.SegDup.min10kb.hg38.GR <- makeGRangesFromDataFrame(SegDup.min10kb.Cento.Telo.hg38,
                                                             keep.extra.columns=TRUE,
                                                             ignore.strand=FALSE,
                                                             seqinfo=NULL,
                                                             seqnames.field=c('seqnames', 'seqname',
                                                                              'chromosome', 'chrom',
                                                                              'Chr', 'chromosome_name',
                                                                              'seqid'),
                                                             start.field='start',
                                                             end.field=c('end','End', 'stop'),
                                                             strand.field='strand',
                                                             starts.in.df.are.0based=FALSE)
# Batch 2 
cat("Run for Batch 2 Deletion","\n")
Batch2.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact <- plyr::compact(Batch2.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All)
Batch2.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact.df <- lapply(Batch2.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact, as_tibble)
Batch2.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.df <- bind_rows(Batch2.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact.df)
Batch2.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.GR <- GRanges(Batch2.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.df)


Drop_Telo_Centro_SegDup.del.results<-suppressWarnings(Find_overlap_with_referenceSets(Batch2.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.GR, Cent.Telo.SegDup.min10kb.hg38.GR))

Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct <- Drop_Telo_Centro_SegDup.del.results[[1]]
Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Del.CentTeloSegDup.df.complete.distinct <- Drop_Telo_Centro_SegDup.del.results[[2]]


############# Run for Batch 1
cat("Run for Batch 1 Deletion","\n")
Batch1.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact <- plyr::compact(Batch1.hg38.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All)
Batch1.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact.df <- lapply(Batch1.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact, as_tibble)
Batch1.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.df <- bind_rows(Batch1.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact.df)
Batch1.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.GR <- GRanges(Batch1.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.df)

Drop_Telo_Centro_SegDup.del.results<-suppressWarnings(Find_overlap_with_referenceSets(Batch1.PennvsQSNP.Del.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.GR, Cent.Telo.SegDup.min10kb.hg38.GR))

Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct <- Drop_Telo_Centro_SegDup.del.results[[1]]
Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Del.CentTeloSegDup.df.complete.distinct <- Drop_Telo_Centro_SegDup.del.results[[2]]

########### Run for Batch 2
cat("Run for Batch 2 Duplication","\n")
Batch2.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact <- plyr::compact(Batch2.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All)
Batch2.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact.df <- lapply(Batch2.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact, as_tibble)
Batch2.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.df <- bind_rows(Batch2.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact.df)
Batch2.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.GR <- GRanges(Batch2.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.df)

Drop_Telo_Centro_SegDup.dup.results<-suppressWarnings(Find_overlap_with_referenceSets(Batch2.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.GR, Cent.Telo.SegDup.min10kb.hg38.GR))

Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct <- Drop_Telo_Centro_SegDup.dup.results[[1]]
Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.CentTeloSegDup.df.complete.distinct <- Drop_Telo_Centro_SegDup.dup.results[[2]]


cat("Run for Batch 1 Duplication","\n")
Batch1.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact <- plyr::compact(Batch1.hg38.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All)
Batch1.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact.df <- lapply(Batch1.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact, as_tibble)
Batch1.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.df <- bind_rows(Batch1.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.compact.df)
Batch1.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.GR <- GRanges(Batch1.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.df)

Drop_Telo_Centro_SegDup.dup.results<-suppressWarnings(Find_overlap_with_referenceSets(Batch1.PennvsQSNP.Dup.Min50Perc.Min10SNPs20kb.Conf15.Intersect.All.GR, Cent.Telo.SegDup.min10kb.hg38.GR))

Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct <- Drop_Telo_Centro_SegDup.dup.results[[1]]
Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.CentTeloSegDup.df.complete.distinct <- Drop_Telo_Centro_SegDup.dup.results[[2]]

#Prep Cleaned Batch2 Penn vs QSNP Ranges without Seg Dup, Cenotromere, Telomere Overlap for Intersect with gnomad common CNVs 
#Reformat Full Location Info so can intersect with details of Neuropsych CNV Annot
Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct$FullLocation_hg38 <- Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct$FullLocation 
Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct$FullLocation_hg38 <- str_replace_all(Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct$FullLocation_hg38, 'chr', '') 

Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct$FullLocation_hg38 <- Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct$FullLocation 
Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct$FullLocation_hg38 <- str_replace_all(Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct$FullLocation_hg38, 'chr', '') 

## Make Cleaned up CNV Calls into GRanges Object Again to Intersect with Known Neuropsychiatric CNVs
Batch2.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.GR <- GRanges(Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct)
Batch2.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.GR <- GRanges(Batch2.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct)


# Prep Cleaned Batch1 Penn vs QSNP Ranges without Seg Dup, Cenotromere, Telomere Overlap for Intersect with gnomad common CNVs
#Reformat Full Location Info so can intersect with details of Neuropsych CNV Annot
Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct$FullLocation_hg38 <- paste(Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct$seqnames, Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct$start, sep = ':')
Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct$FullLocation_hg38 <- paste(Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct$FullLocation_hg38, Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct$end, sep = '-')

Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct$FullLocation_hg38 <- paste(Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct$seqnames, Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct$start, sep = ':')
Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct$FullLocation_hg38 <- paste(Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct$FullLocation_hg38, Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct$end, sep = '-')

## Make Cleaned up CNV Calls into GRanges Object Again to Intersect with gnomad common CNVs
Batch1.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.GR <- GRanges(Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Del.NoCentTeloSegDup.df.complete.distinct)
Batch1.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.GR <- GRanges(Batch1.PennvsQSNP.MinConf15.Min10SNP20kb.Dup.NoCentTeloSegDup.df.complete.distinct)

cat("\n############################################## Dropping common CNVs based on gnomad frequency data #########################################\n")
#Read in gnomad common dels + dups (Pop Max Freq < .001) with hg38 coordinates
gnomad.v2.CNVs.common001.complete.min10kb.DelorDup.LiftOverAvail.hg38 <- read.xlsx('../GeneticReferenceData/gnomad.v2.CNVs.common001.complete.min10kb.DelorDup.LiftOverAvail.hg38.042323.xlsx')
#Extract subset of common CNVs to check for overlap with CNV Calls to drop
gnomad.v2.CNVs.common01.complete.min10kb.DelorDup.LiftOverAvail.hg38 <- gnomad.v2.CNVs.common001.complete.min10kb.DelorDup.LiftOverAvail.hg38[gnomad.v2.CNVs.common001.complete.min10kb.DelorDup.LiftOverAvail.hg38$POPMAX_AF >= 0.01,]
gnomad.v2.CNVs.common01.Del.min10kb.hg38 <- gnomad.v2.CNVs.common01.complete.min10kb.DelorDup.LiftOverAvail.hg38[gnomad.v2.CNVs.common01.complete.min10kb.DelorDup.LiftOverAvail.hg38$SVTYPE == 'DEL',]
gnomad.v2.CNVs.common01.Del.min10kb.hg38.reduce <- dplyr::select(gnomad.v2.CNVs.common01.Del.min10kb.hg38, chrom, start, end, old.chrom, old.start, old.end, POPMAX_AF, SVLEN, length, SVTYPE, PROTEIN_CODING__LOF, PROTEIN_CODING__INTRONIC, PROTEIN_CODING__INTERGENIC)
gnomad.v2.CNVs.common01.Dup.min10kb.hg38 <- gnomad.v2.CNVs.common01.complete.min10kb.DelorDup.LiftOverAvail.hg38[gnomad.v2.CNVs.common01.complete.min10kb.DelorDup.LiftOverAvail.hg38$SVTYPE == 'DUP',]
gnomad.v2.CNVs.common01.Dup.min10kb.hg38.reduce <- dplyr::select(gnomad.v2.CNVs.common01.Dup.min10kb.hg38, chrom, start, end, old.chrom, old.start, old.end, POPMAX_AF, SVLEN, length, SVTYPE, PROTEIN_CODING__LOF, PROTEIN_CODING__INTRONIC, PROTEIN_CODING__INTERGENIC)



gnomad.common01.Del.min10kb.hg38.GR <- makeGRangesFromDataFrame(gnomad.v2.CNVs.common01.Del.min10kb.hg38.reduce,
                                                                keep.extra.columns=TRUE,
                                                                ignore.strand=FALSE,
                                                                seqinfo=NULL,
                                                                seqnames.field=c( 'chrom',
                                                                                  'chr'),
                                                                start.field='start',
                                                                end.field=c('end','End'),
                                                                strand.field='strand',
                                                                starts.in.df.are.0based=FALSE)

gnomad.common01.Dup.min10kb.hg38.GR <- makeGRangesFromDataFrame(gnomad.v2.CNVs.common01.Dup.min10kb.hg38.reduce,
                                                                keep.extra.columns=TRUE,
                                                                ignore.strand=FALSE,
                                                                seqinfo=NULL,
                                                                seqnames.field=c( 'chrom',
                                                                                  'chr'),
                                                                start.field='start',
                                                                end.field=c('end','End'),
                                                                strand.field='strand',
                                                                starts.in.df.are.0based=FALSE)

cat("Run for Batch 2 Deletion","\n")
Del.overlapWithGnomad.results<-suppressWarnings(Find_overlap_with_referenceSets(Batch2.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.GR, gnomad.common01.Del.min10kb.hg38.GR))
# saveRDS(Del.overlapWithGnomad.results, './intermediate_results/Del.overlapWithGnomad.results.rda')
Batch2.PennvsQSNP.Del.Intersect.AllSubj.RareCNV.df.complete.distinct<-Del.overlapWithGnomad.results[[1]]
Batch2.PennvsQSNP.Del.Intersect.AllSubj.Common01.CNVHits.df.complete<-Del.overlapWithGnomad.results[[2]]

Batch2.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.GR <- GRanges(Batch2.PennvsQSNP.Del.Intersect.AllSubj.RareCNV.df.complete.distinct)


cat("Run for Batch 1 Deletion","\n")
Del.overlapWithGnomad.results<-suppressWarnings(Find_overlap_with_referenceSets(Batch1.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.GR, gnomad.common01.Del.min10kb.hg38.GR))
# saveRDS(Del.overlapWithGnomad.results, './intermediate_results/Del.overlapWithGnomad.results.rda')
Batch1.PennvsQSNP.Del.Intersect.AllSubj.RareCNV.df.complete.distinct<-Del.overlapWithGnomad.results[[1]]
Batch1.PennvsQSNP.Del.Intersect.AllSubj.Common01.CNVHits.df.complete<-Del.overlapWithGnomad.results[[2]]
Batch1.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.GR <- GRanges(Batch1.PennvsQSNP.Del.Intersect.AllSubj.RareCNV.df.complete.distinct)



cat("Run for Batch 2 Duplication","\n")
Dup.overlapWithGnomad.results<-suppressWarnings(Find_overlap_with_referenceSets(Batch2.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.GR, gnomad.common01.Dup.min10kb.hg38.GR))
# saveRDS(Del.overlapWithGnomad.results, './intermediate_results/Del.overlapWithGnomad.results.rda')
Batch2.PennvsQSNP.Dup.Intersect.AllSubj.RareCNV.df.complete.distinct<-Dup.overlapWithGnomad.results[[1]]
Batch2.PennvsQSNP.Dup.Intersect.AllSubj.Common01.CNVHits.df.complete<-Dup.overlapWithGnomad.results[[2]]
Batch2.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.GR <- GRanges(Batch2.PennvsQSNP.Dup.Intersect.AllSubj.RareCNV.df.complete.distinct)




cat("Run for Batch 1 Duplication","\n")
Dup.overlapWithGnomad.results<-suppressWarnings(Find_overlap_with_referenceSets(Batch1.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.GR, gnomad.common01.Dup.min10kb.hg38.GR))
# saveRDS(Del.overlapWithGnomad.results, './intermediate_results/Del.overlapWithGnomad.results.rda')
Batch1.PennvsQSNP.Dup.Intersect.AllSubj.RareCNV.df.complete.distinct<-Dup.overlapWithGnomad.results[[1]]
Batch1.PennvsQSNP.Dup.Intersect.AllSubj.Common01.CNVHits.df.complete<-Dup.overlapWithGnomad.results[[2]]
Batch1.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.GR <- GRanges(Batch1.PennvsQSNP.Dup.Intersect.AllSubj.RareCNV.df.complete.distinct)

cat("\n#################################################### Prep to Intersect with Know CNV Loci ##################################################\n")
DD.SZ.ASD.CNVs <- read.xlsx('../GeneticReferenceData/DDD_Expert66_Coe2014FixedNDDDirection_CNV_GenomicDisorders_DelDup_PlusBrainSpanPaperLoci_hg19_Liftover_hg38_Fix1q21DupSCZ_062524.xlsx')
#Drop regions that couldn't be lifted over
DD.SZ.ASD.CNVs <- DD.SZ.ASD.CNVs[!is.na(DD.SZ.ASD.CNVs$UCSC_LiftOver_hg38_min0.8_Overlap),]
DD.SZ.ASD.CNVs <- DD.SZ.ASD.CNVs[!DD.SZ.ASD.CNVs$LossvGain == "NOTSIG",]
DD.SZ.ASD.CNVs <- DD.SZ.ASD.CNVs[!DD.SZ.ASD.CNVs$`Significance.in.Orig.Paper/Resource` == "NOT SIG",]


DD.SZ.ASD.CNVs$width <- DD.SZ.ASD.CNVs$hg38_end - DD.SZ.ASD.CNVs$hg38_start
DD.SZ.ASD.CNVs$HalfLoci <- round(DD.SZ.ASD.CNVs$width/2, 0)
DD.SZ.ASD.CNVs$FortyPercLoci <- round(DD.SZ.ASD.CNVs$width*0.4, 0)
DD.SZ.ASD.CNVs$FullLocation_hg38 <- DD.SZ.ASD.CNVs$UCSC_LiftOver_hg38_min0.8_Overlap

DD.SZ.ASD.CNVs <- DD.SZ.ASD.CNVs[!is.na(DD.SZ.ASD.CNVs$LossvGain),]

DD.SZ.ASD.CNVs$chr <- DD.SZ.ASD.CNVs$hg38_chr
DD.SZ.ASD.CNVs$start <- DD.SZ.ASD.CNVs$hg38_start
DD.SZ.ASD.CNVs$end <- DD.SZ.ASD.CNVs$hg38_end

DD.SZ.ASD.CNVs.Del <- DD.SZ.ASD.CNVs[DD.SZ.ASD.CNVs$LossvGain == 'Loss',]
SZ.CNVs.Del <- DD.SZ.ASD.CNVs.Del[DD.SZ.ASD.CNVs.Del$DiseaseAssociation_SCZMarshall_ASDSanders_DDD_NDDCoe == 'SCZ',]

DD.SZ.ASD.CNVs.Dup <- DD.SZ.ASD.CNVs[DD.SZ.ASD.CNVs$LossvGain == 'Gain',]
SZ.CNVs.Dup <- DD.SZ.ASD.CNVs.Dup[DD.SZ.ASD.CNVs.Dup$DiseaseAssociation_SCZMarshall_ASDSanders_DDD_NDDCoe == 'SCZ',]

NRXN1.refseq.hg38.exons <- read.xlsx('../GeneticReferenceData/NRXN1_Exons_hg38_Full_092723.xlsx')
colnames(NRXN1.refseq.hg38.exons) <- c("chr", "start", "end", "transcript")
NRXN1.refseq.hg38.exons$length <- NRXN1.refseq.hg38.exons$end - NRXN1.refseq.hg38.exons$start
NRXN1.refseq.hg38.exons$FiftyPercExon <- NRXN1.refseq.hg38.exons$length/2

#Create GR Objects for Known Neuropsychiatric CNVs
KnownNDD.Del.GR.p1 <- makeGRangesFromDataFrame(DD.SZ.ASD.CNVs.Del,
                                               keep.extra.columns=TRUE,
                                               ignore.strand=FALSE,
                                               seqinfo=NULL,
                                               seqnames.field=c('seqnames', 'seqname',
                                                                'chromosome', 'chrom',
                                                                'Chr', 'chromosome_name',
                                                                'seqid'),
                                               start.field='start',
                                               end.field=c('end','End', 'stop'),
                                               strand.field='strand',
                                               starts.in.df.are.0based=FALSE)

KnownNDD.Dup.GR.p1 <- makeGRangesFromDataFrame(DD.SZ.ASD.CNVs.Dup,
                                               keep.extra.columns=TRUE,
                                               ignore.strand=FALSE,
                                               seqinfo=NULL,
                                               seqnames.field=c('seqnames', 'seqname',
                                                                'chromosome', 'chrom',
                                                                'Chr', 'chromosome_name',
                                                                'seqid'),
                                               start.field='start',
                                               end.field=c('end', 'End','stop'),
                                               strand.field='strand',
                                               starts.in.df.are.0based=FALSE)


NRXN1.exons.hg38.GR <- makeGRangesFromDataFrame(NRXN1.refseq.hg38.exons,
                                                keep.extra.columns=TRUE,
                                                ignore.strand=FALSE,
                                                seqinfo=NULL,
                                                seqnames.field=c('seqnames', 'seqname',
                                                                 'chromosome', 'chrom',
                                                                 'Chr', 'chromosome_name',
                                                                 'seqid'),
                                                start.field='start',
                                                end.field=c('end','End', 'stop'),
                                                strand.field='strand',
                                                starts.in.df.are.0based=FALSE)


cat("Run for Batch 2 Deletion","\n")
Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete<-Find_overlap_with_KnownRiskLoci(KnownNDD.Del.GR.p1, Batch2.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.GR)

Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete$ID.Short <- Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete$ID
Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete$ID.Short <- str_replace_all(Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete$ID.Short, '_rep', '')
Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete.distinct <- distinct(Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete, ID, seqnames, start, end, .keep_all = TRUE)
Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.df.SCZ <- Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete[Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete$DiseaseAssoc == 'SCZ',]
Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.df.ASD <- Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete[Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete$DiseaseAssoc == 'ASD',]



cat("Run for Batch 1 Deletion","\n")
Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete<-Find_overlap_with_KnownRiskLoci(KnownNDD.Del.GR.p1, Batch1.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.GR)

Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete$ID.Short <- Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete$ID
Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete$ID.Short <- str_replace_all(Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete$ID, '_RPT', '')
Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete.distinct <- distinct(Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete, ID, seqnames, start, end, .keep_all = TRUE)
Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.df.SCZ <- Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete[Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete$DiseaseAssoc == 'SCZ',]
Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.df.ASD <- Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete[Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete$DiseaseAssoc == 'ASD',]



cat("Run for Batch 2 Duplication","\n")

Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete<-Find_overlap_with_KnownRiskLoci(KnownNDD.Dup.GR.p1, Batch2.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.GR)
Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete$ID.Short <- Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete$ID
Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete$ID.Short <- str_replace_all(Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete$ID.Short, '_rep', '')
Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete.distinct <- distinct(Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete, ID, seqnames, start, end, .keep_all = TRUE)
Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.df.SCZ <- Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete[Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete$DiseaseAssoc == 'SCZ',]
Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.df.ASD <- Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete[Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete$DiseaseAssoc == 'ASD',]

cat("Run for Batch 1 Duplication","\n")
Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete<-Find_overlap_with_KnownRiskLoci(KnownNDD.Dup.GR.p1, Batch1.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.GR)
Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete$ID.Short <- Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete$ID
Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete$ID.Short <- str_replace_all(Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete$ID.Short, '_rep', '')
Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete.distinct <- distinct(Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete, ID, seqnames, start, end, .keep_all = TRUE)
Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.df.SCZ <- Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete[Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete$DiseaseAssoc == 'SCZ',]
Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.df.ASD <- Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete[Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete$DiseaseAssoc == 'ASD',]

cat("#################################################### Check Deletion Overlap with NRXN1 Exons ###############################################\n")
Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NRXN1.Del.Loci.Hits.df.complete<-Find_overlap_with_NRXN1_Exons(NRXN1.exons.hg38.GR, Batch2.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.GR)


Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NRXN1.Del.Loci.Hits.df.complete<-Find_overlap_with_NRXN1_Exons(NRXN1.exons.hg38.GR, Batch1.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.GR)

Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NRXN1.Del.Loci.Hits.df.complete$ID.Short <- Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NRXN1.Del.Loci.Hits.df.complete$ID
Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NRXN1.Del.Loci.Hits.df.complete$SyndromeLocation <- '2:50147488-51259674'

cat("################################################### Concatenate Known CNV Loci #############################################################\n")

Batch12.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete <- bind_rows(Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete, Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete)

Batch12.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.NDDSCZASD.DelDup.Loci.Hits.df.complete <- bind_rows(Batch12.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete, Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete) %>% bind_rows(Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete)

Batch12.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.NDDSCZASD.DelDup.Loci.Hits.df.complete <- bind_rows(Batch12.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.NDDSCZASD.DelDup.Loci.Hits.df.complete, Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NRXN1.Del.Loci.Hits.df.complete) %>% bind_rows(Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NRXN1.Del.Loci.Hits.df.complete)
Batch12.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.NDDSCZASD.DelDup.Loci.Hits.df.complete.distinct <- distinct(Batch12.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.NDDSCZASD.DelDup.Loci.Hits.df.complete)


cat("############################### Split GR Objects by Subject, Reduce Genomic Ranges and Get Length Per Subj #################################\n")

cat("Run for Batch 2 Deletion","\n")
Batch2.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.Reduced.AllSub.df.complete<-Split_GRObjects_bySubject(Batch2.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.GR)
Batch2.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj <- group_by(Batch2.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.Reduced.AllSub.df.complete, ID) %>% summarise(sum(width))
colnames(Batch2.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj)[2] <- 'Total_Del_Length_min10SNP20kb_NoCentTeloSegDup'

cat("Run for Batch 1 Deletion","\n")
Batch1.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.Reduced.AllSub.df.complete<-Split_GRObjects_bySubject(Batch1.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.GR)
Batch1.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj <- group_by(Batch1.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.Reduced.AllSub.df.complete, ID) %>% summarise(sum(width))
colnames(Batch1.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj)[2] <- 'Total_Del_Length_min10SNP20kb_NoCentTeloSegDup'

cat("Run for Batch 2 Duplication")

Batch2.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.Reduced.AllSub.df.complete<-Split_GRObjects_bySubject(Batch2.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.GR)
Batch2.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj <- group_by(Batch2.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.Reduced.AllSub.df.complete, ID) %>% summarise(sum(width))
colnames(Batch2.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj)[2] <- 'Total_Dup_Length_min10SNP20kb_NoCentTeloSegDup'


cat("Run for Batch 1 Duplication","\n")
Batch1.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.Reduced.AllSub.df.complete<-Split_GRObjects_bySubject(Batch1.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.GR)
Batch1.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj <- group_by(Batch1.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.Reduced.AllSub.df.complete, ID) %>% summarise(sum(width))
colnames(Batch1.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj)[2] <- 'Total_Dup_Length_min10SNP20kb_NoCentTeloSegDup'


Batch12.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj<-rbind(Batch2.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj,Batch1.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj )

Batch12.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj<-rbind(Batch2.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj,Batch1.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj )

cat("############################################ Annotate Rare CNVs by RefSeq Protein-Coding Genes #############################################\n")
MainChroms<-1:22
MainChroms<-c(MainChroms, "X","Y")

hg38.refseq.exons.distinctbyHGNC.coding<-read_tsv('../GeneticReferenceData/hg38_ncbiRefSeq_Exons_UCSC_distinctbyHGNC_coding.txt',show_col_types = FALSE)
hg38.refseq.exons.distinctbyHGNC.coding.MainChroms <- hg38.refseq.exons.distinctbyHGNC.coding[hg38.refseq.exons.distinctbyHGNC.coding$chrom %in% MainChroms,]
hg38.refseq.exons.distinctbyHGNC.coding.MainChroms$chrom <- str_replace_all(hg38.refseq.exons.distinctbyHGNC.coding.MainChroms$chrom, 'chr', '')
colnames(hg38.refseq.exons.distinctbyHGNC.coding.MainChroms)[4] <- 'HGNCUpdatedSymbol'

#Make RefSeq Genomic Range
Refseq.ProtCoding.GR <- makeGRangesFromDataFrame(hg38.refseq.exons.distinctbyHGNC.coding.MainChroms,
                                                 keep.extra.columns=TRUE,
                                                 ignore.strand=FALSE,
                                                 seqinfo=NULL,
                                                 seqnames.field=c('seqnames', 'seqname',
                                                                  'chromosome', 'chrom',
                                                                  'Chr', 'chromosome_name',
                                                                  'seqid'),
                                                 start.field='exonStarts',
                                                 end.field=c('end', 'End','stop', 'exonEnds'),
                                                 strand.field='strand',
                                                 starts.in.df.are.0based=FALSE)

############################################ Prep for Kang Modules, PLI, ASD/SCZ Gene Annotations  #######################################

Kangmodules.minmod30.DS2.Genes <- read.xlsx( '../GeneticReferenceData/Kangmodules.minmod30.DS2.wBDPPI_Updated2021.xlsx')
Kangmodules.minmod30.DS2.Genes.NoDupe <- Kangmodules.minmod30.DS2.Genes[!duplicated(Kangmodules.minmod30.DS2.Genes$HGNCUpdatedSymbol),]

Kangmodules.minmod30.DS2.ModuleNumNames <- unique(Kangmodules.minmod30.DS2.Genes.NoDupe$ModuleNumName)
Kangmodules.minmod30.DS2.ModuleNumNames <- Kangmodules.minmod30.DS2.ModuleNumNames[-9]

Kangmodules.minmod30.DS2.wBDPPI <- Kangmodules.minmod30.DS2.Genes
Kangmodules.minmod30.DS2.wBDPPI.nodupe <- Kangmodules.minmod30.DS2.wBDPPI[!duplicated(Kangmodules.minmod30.DS2.wBDPPI$HGNCUpdatedSymbol),]

Kangmodules.minmod30.DS2.key <- select(Kangmodules.minmod30.DS2.wBDPPI.nodupe, OrigSymbol, HGNCUpdatedSymbol)

NewpLIList.canonical.reduce <- read.xlsx('../GeneticReferenceData/Karczewski2020.pLI.canonical.UpdatedHGNCSymbol.110723.xlsx')


NewpLIList.canonical.forannot <- dplyr::select(NewpLIList.canonical.reduce, HGNCUpdatedSymbol, ensembl_gene_id, oe_lof_upper, LOEUF.perc, LOEUF.top10perc, LOEUF.top20perc)

#Check for dupe HGNC Symbols and select highest LOEUF scored gene
NewpLIList.dupegenes <- NewpLIList.canonical.forannot$HGNCUpdatedSymbol[duplicated(NewpLIList.canonical.forannot$HGNCUpdatedSymbol)]
NewpLIList.canonical.forannot.dupegenes <- NewpLIList.canonical.forannot[NewpLIList.canonical.forannot$HGNCUpdatedSymbol %in% NewpLIList.dupegenes,]

NewpLIList.canonical.sortLOEUF.forannot <- NewpLIList.canonical.forannot[order(NewpLIList.canonical.forannot$LOEUF.perc, decreasing = TRUE),]
NewpLIList.canonical.mostLOEUF.forannot <- NewpLIList.canonical.sortLOEUF.forannot[!duplicated(NewpLIList.canonical.sortLOEUF.forannot$HGNCUpdatedSymbol),]

DD.ASD.SCZ.genes.complete.wide <- read.xlsx( "../GeneticReferenceData/DD.ASD.SCZ.genes.complete.wide.112723.xlsx")
DD.ASD.SCZ.genes.complete.wide <- select(DD.ASD.SCZ.genes.complete.wide, -OrigSymbol)

cat("################################# Get Genes for Genomic Ranges for Min 5 SNP 10 kb Dup CNV Overlap #########################################\n")
cat("Run for Batch 1 & 2 Deletion","\n")
Batch12.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.GR <- c(Batch1.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.GR, Batch2.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.GR)
Batch12.MinConf15.10SNP20kb.RareFullExonDel<-suppressWarnings(Find_overlap_with_refSeq(Batch12.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.GR, Refseq.ProtCoding.GR))

Batch12.MinConf15.10SNP20kb.RareFullExonDel.DistinctGene <- distinct(Batch12.MinConf15.10SNP20kb.RareFullExonDel, ID, HGNCUpdatedSymbol, .keep_all = TRUE)
Batch12.MinConf15.10SNP20kb.RareFullExonDel.wBrainSpan <- suppressMessages(left_join(Batch12.MinConf15.10SNP20kb.RareFullExonDel.DistinctGene, Kangmodules.minmod30.DS2.wBDPPI.nodupe))
Batch12.MinConf15.10SNP20kb.RareFullExonDel.wBrainSpan.PLI <- suppressMessages(left_join(Batch12.MinConf15.10SNP20kb.RareFullExonDel.wBrainSpan, NewpLIList.canonical.mostLOEUF.forannot))
Batch12.MinConf15.10SNP20kb.RareFullExonDel.wBrainSpan.PLI$LossvGain <- 'Loss'
Batch12.MinConf15.10SNP20kb.RareFullExonDel.wBrainSpan.PLI.DxGenes <- suppressMessages(left_join(Batch12.MinConf15.10SNP20kb.RareFullExonDel.wBrainSpan.PLI, DD.ASD.SCZ.genes.complete.wide))

cat("Run for Batch 1 & 2 Duplication","\n")
Batch12.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.GR <- c(Batch1.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.GR, Batch2.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.GR)

Batch12.MinConf15.10SNP20kb.RareFullExonDup<-suppressWarnings(Find_overlap_with_refSeq(Batch12.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.GR, Refseq.ProtCoding.GR))
Batch12.MinConf15.10SNP20kb.RareFullExonDup.DistinctGene <- distinct(Batch12.MinConf15.10SNP20kb.RareFullExonDup, ID, HGNCUpdatedSymbol, .keep_all = TRUE)
Batch12.MinConf15.10SNP20kb.RareFullExonDup.wBrainSpan <- suppressMessages(left_join(Batch12.MinConf15.10SNP20kb.RareFullExonDup.DistinctGene, Kangmodules.minmod30.DS2.wBDPPI.nodupe))
Batch12.MinConf15.10SNP20kb.RareFullExonDup.wBrainSpan.PLI <- suppressMessages(left_join(Batch12.MinConf15.10SNP20kb.RareFullExonDup.wBrainSpan, NewpLIList.canonical.mostLOEUF.forannot))
Batch12.MinConf15.10SNP20kb.RareFullExonDup.wBrainSpan.PLI$LossvGain <- 'Gain'

cat("############################################ Make Summary Stats for Rare Dels by Subj for Analysis #########################################\n")
Batch12.PennCNV.QC.GSAMD.Outlier.distinct<-read.xlsx("../input_files/Batch12.PennCNV.QC.GSAMD.Outlier.distinct.xlsx")

AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDel.Summarystats<-suppressMessages(Make_Summary_Stats_for_RareCNV_bySubj(Batch12.MinConf15.10SNP20kb.RareFullExonDel.wBrainSpan.PLI, "Del", Batch12.PennCNV.QC.GSAMD.Outlier.distinct))

AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDup.Summarystats<-suppressMessages(Make_Summary_Stats_for_RareCNV_bySubj(Batch12.MinConf15.10SNP20kb.RareFullExonDup.wBrainSpan.PLI, "Dup", Batch12.PennCNV.QC.GSAMD.Outlier.distinct))
AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats <- suppressMessages(left_join(AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDel.Summarystats, AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDup.Summarystats))

colnames.RareFullExonDelDup <- data.frame(ColName = colnames(AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats))
colnames.RareFullExonDelDup$ColName <- str_replace_all(colnames.RareFullExonDelDup$ColName, "_FullExonDel_", "_Del_")
colnames.RareFullExonDelDup$ColName <- str_replace_all(colnames.RareFullExonDelDup$ColName,  "_FullExonDup_", "_Dup_")
colnames(AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats) <- colnames.RareFullExonDelDup$ColName

AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats<-suppressMessages(left_join(AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats, Batch12.PennvsQSNP.MinConf15.10SNP20kb.Dup.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj)%>%left_join(Batch12.PennvsQSNP.MinConf15.10SNP20kb.Del.NoCentTeloSegDup.RareCNV.ReducedGR.TotalWidth.bySubj))



AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$SCZ_Del_012424 <- 0
AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$SCZ_Del_012424[AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$ID %in% Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.df.SCZ$ID.Short | AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$ID %in% Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.df.SCZ$ID.Short | AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$ID %in% Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NRXN1.Del.Loci.Hits.df.complete$ID |  AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$ID %in% Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NRXN1.Del.Loci.Hits.df.complete$ID] <- 1
AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$NDDSCZASD_Del_012424 <- 0
AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$NDDSCZASD_Del_012424[AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$ID %in% Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete.distinct$ID.Short | AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$ID %in% Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NDDSCZASD.Del.Loci.Hits.df.complete.distinct$ID.Short | AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$ID %in% Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NRXN1.Del.Loci.Hits.df.complete$ID.Short | AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$ID %in% Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Del.NRXN1.Del.Loci.Hits.df.complete$ID] <- 1

AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$SCZ_Dup_012424 <- 0
AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$SCZ_Dup_012424[AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$ID %in% Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.df.SCZ$ID.Short | AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$ID %in% Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.df.SCZ$ID.Short] <- 1
AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$NDDSCZASD_Dup_012424 <- 0
AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$NDDSCZASD_Dup_012424[AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$ID %in% Batch1.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete.distinct$ID.Short | AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$ID %in% Batch2.PennvsQSNP.MinConf15.10SNP20kb.NoCentTeloSegDup.RareCNV.Dup.NDDSCZASD.Dup.Loci.Hits.df.complete.distinct$ID.Short ] <- 1

AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$SCZ_DelorDup_CNV_012424 <- 0 
AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$SCZ_DelorDup_CNV_012424[AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$SCZ_Del_012424 == 1 | AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$SCZ_Dup_012424 == 1] <- 1
AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$NDDSCZASD_DelorDup_CNV_012424 <- 0 
AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$NDDSCZASD_DelorDup_CNV_012424[AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$NDDSCZASD_Del_012424 == 1 | AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats$NDDSCZASD_Dup_012424 == 1] <- 1

AllLargeCNV_bySubj_min10SNPs20kb.DelRareFullExonDup.Summarystats.NAFix<-AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats
AllLargeCNV_bySubj_min10SNPs20kb.DelRareFullExonDup.Summarystats.NAFix[,1:ncol(AllLargeCNV_bySubj_min10SNPs20kb.DelRareFullExonDup.Summarystats.NAFix)][is.na(AllLargeCNV_bySubj_min10SNPs20kb.DelRareFullExonDup.Summarystats.NAFix[,1:ncol(AllLargeCNV_bySubj_min10SNPs20kb.DelRareFullExonDup.Summarystats.NAFix)])] <- 0

write.xlsx(AllLargeCNV_bySubj_min10SNPs20kb.DelRareFullExonDup.Summarystats.NAFix, "../output/AllLargeCNV_bySubj_min10SNPs20kb.RareFullExonDelDup.Summarystats_070524.xlsx")


