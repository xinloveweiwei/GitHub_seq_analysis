rm(list = ls())
directory <- '/Users/macbook/Documents/project/DPF2/result/seq_data/RNA_seq/RAW_hjx/1stwt_E7_3_repeats_202_E7_LPS_CRO/counts_fromweiwei'
setwd(directory)
outputPrefix <- 'Gmax_DESeq2' 

getwd()
list.files()

sampleFiles <- list.files(pattern = '*_0-0-*')

sampleNames <- gsub('*-output_basename.counts$', "",sampleFiles)

sampleCondition <- c('ctr','ctr','ctr','ctr','KD','KD','KD','KD')

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

treatments <- c('ctr','KD')

library(DESeq2)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)
View(colData(ddsHTSeq)$condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

dds <- DESeq(ddsHTSeq)

res <- results(dds)
res2 <- subset(res, res$padj<0.05)

#res2 <- as.data.frame(res2) 可以在对话框中看到数据框的内容

# order results by padj value (most significant to least)
res2 <- res2[order(res2$padj),]
#res2中没有count信息，还在dds中。
write.csv(res2, file = 'res2.csv', na='')

tem <- counts(dds, normalized =T)
View(tem)

# save data results and normalized reads to csv
#counts(dds,normalized = T) to get the normalized count data from dds.
resdata <- merge(as.data.frame(res2), as.data.frame(counts(dds, normalized =T)), 
                 by = "row.names", sort = F)
View(resdata)
names(resdata)[1] <- 'gene'   #names() 等于 colnames(). rownames 等于row.names.

write.csv(resdata, file = paste0(outputPrefix, 'of 0_0_WTvsKD.csv'))
