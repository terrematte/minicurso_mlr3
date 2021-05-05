library("DESeq2")

load("../data/tcga_kirc.RData")

codes.metas <- kirc_cli %>% 
  dplyr::filter(metastasis %in% c("M0", "M1")) %>%    
  droplevels()  %>%
  rownames_to_column("samples") %>%
  dplyr::select(samples)  %>%
  pull(.)

data <- round(kirc_rna+1, digits = 0)

countData <- as.matrix(t(data[codes.metas,]))
colData <- droplevels(kirc_cli[codes.metas,])

# create the DESeqDataSet object
ddsObj <- DESeqDataSetFromMatrix(countData = countData,
                                 colData = colData,
                                 design = ~ metastasis)

ddsObj <- DESeq(ddsObj)
res.shr  <- results(ddsObj)
summary(res.shr)

dea.M0.M1 <- as.data.frame(res.shr) %>%
  rownames_to_column("symbol") %>% 
  dplyr::rename(logFC=log2FoldChange, FDR=padj)

df.deseq <- dea.M0.M1  %>% filter(abs(logFC) >=2, pvalue <= 0.01)

df.deseq <- df.deseq[!is.na(df.deseq$logFC),]
dim(df.deseq)


save(dea.M0.M1, file = "../data/dea.M0.M1.rda", compress = T)

dea.M0.M1.lst <- unique(df.deseq$symbol)
write(dea.M0.M1.lst,  file = "../data/dea.M0.M1.lst")


# Keeg Genes
map05211 <- read.csv("../data/map05211.tsv", sep=";", comment.char="#", stringsAsFactors = F)
write(map05211$Gene, file = "../data/genes_keeg.lst")


