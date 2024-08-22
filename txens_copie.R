library(AnnotationHub)
library(ensembldb)

ah <- AnnotationHub()
edb=query(ah, c("EnsDb", "Hsapiens", "105"))[[1]]

tx_lens=transcriptLengths(edb,with.utr5_len = TRUE,with.utr3_len = TRUE)
head(tx_lens)
