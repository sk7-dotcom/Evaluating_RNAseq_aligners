### Counting HiSAT2 data with featurCounts 


str_counts <- featureCounts(bam_files, annot.ext = 'mm10.ref.gtf', isGTFAnnotationFile = TRUE, GTF.featureType = 'gene', GTF.attrType = 'gene_id', countMultiMappingReads = TRUE)

write.csv(str_counts[["counts"]], 'str_counts.csv')

countData <- read.csv('str_count.csv', header = TRUE)

countName <- c("id",
               "day14_immuno_treat_1", 
               "day14_immuno_treat_2", 
               "day14_immuno_treat_3",
               "day14_immuno_treat_4", 
               "day14_immuno_treat_5", 
               "day14_isotype_control_1", 
               "day14_isotype_control_2",
               "day14_isotype_control_3", 
               "day14_isotype_control_4",
               "day14_isotype_control_5",
               "day7_immuno_treat_1", 
               "day7_immuno_treat_2", 
               "day7_immuno_treat_3",
               "day7_immuno_treat_4", 
               "day7_immuno_treat_5", 
               "day7_isotype_control_1", 
               "day7_isotype_control_2",
               "day7_isotype_control_3", 
               "day7_isotype_control_4",
               "day7_isotype_control_5")
names(countData) <- countName
countData$id <- gsub("\\..*", "",countData$id)

# Experimental Setup

Sample <- c("day14_immuno_treat_1", 
            "day14_immuno_treat_2", 
            "day14_immuno_treat_3",
            "day14_immuno_treat_4", 
            "day14_immuno_treat_5", 
            "day14_isotype_control_1", 
            "day14_isotype_control_2",
            "day14_isotype_control_3", 
            "day14_isotype_control_4",
            "day14_isotype_control_5",
            "day7_immuno_treat_1", 
            "day7_immuno_treat_2", 
            "day7_immuno_treat_3",
            "day7_immuno_treat_4", 
            "day7_immuno_treat_5", 
            "day7_isotype_control_1", 
            "day7_isotype_control_2",
            "day7_isotype_control_3", 
            "day7_isotype_control_4",
            "day7_isotype_control_5")

dpi <- c(rep(14, 10), rep(7, 10))
condition <- c(rep('Treatment', 5), rep('Control', 5), rep('Treatment', 5), rep('Control', 5))

group <- c("day14_treat", 
           "day14_treat", 
           "day14_treat",
           "day14_treat", 
           "day14_treat", 
           "day14_control", 
           "day14_control",
           "day14_control", 
           "day14_control",
           "day14_control",
           "day7_treat", 
           "day7_treat", 
           "day7_treat",
           "day7_treat", 
           "day7_treat", 
           "day7_control", 
           "day7_control",
           "day7_control", 
           "day7_control",
           "day7_control")
metaData<- cbind(Sample, dpi, condition, group)

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~group, tidy = TRUE)

# dds$group <- relevel(dds$group, ref = "day7_control")
hisat_final <- DESeq(dds, test="LRT", reduced = ~ 1)

hisat_day7_comp <- results(hisat_final, contrast = list('group_day7_treat_vs_day14_control', 'group_day7_control_vs_day14_control'),
                           independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE)

hisat_day14_comp <- results(hisat_final, contrast = list('group_day14_treat_vs_day14_control'),
                            independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE)

write.csv(hisat_day14_comp, "hisat_day14.csv")
write.csv(hisat_day7_comp, "hisat_day7.csv")
