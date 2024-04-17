#R 4.3.2 

library(edgeR) #version 3.42.4
library(clusterProfiler) #version 4.8.3 
library(ggplot2)

matrix <- read.table("raw_gene_counts.tsv",sep = "\t",quote = NULL,header = T,
                       row.names = 1,check.names = FALSE)


size <- colSums(matrix)
size

y <- DGEList(counts=matrix)

group <- as.factor(c("control", "control", "control", "KO", "KO", "KO", "em24h", "em24h", "em24h","em48h", "em48h", "em48h"))
y$samples$group <- group
y$samples

#filtering 
table(rowSums(y$counts==0)==12)
keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]
dim(y)

#pre-normalization 
logcpm_before <- cpm(y, log=TRUE)
boxplot(logcpm_before)

#normalization 
y <- calcNormFactors(y, method = "TMM")

logcpm <- cpm(y, log=TRUE, prior.count = 1)
boxplot(logcpm)

#design matrix 
design <- model.matrix(~group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design 


plotMDS(logcpm, labels=group)

y <- estimateDisp(y, design)
plotBCV(y)
y
y$common.dispersion

#differential expression analysis - em24h vs control 
fit <- glmQLFit(y, design)
qlf.2vs1 <- glmQLFTest(fit, coef=2) 
topTags(qlf.2vs1)

FDR <- p.adjust(qlf.2vs1$table$PValue, method="BH")

summary(decideTests(qlf.2vs1, p.value=0.01, lfc=1)) 
results_em24h <- topTags(qlf.2vs1, n=20000, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)$table
results_em24h <-(results_em24h[which(results_em24h$FDR < 0.01 & (results_em24h$logFC >1 | results_em24h$logFC < -1)), ])
write.table(results_em24h, file= 'results_em24h.txt',sep='\t')

#differential expression analysis - em48h vs control
fit <- glmQLFit(y, design)
qlf.2vs1 <- glmQLFTest(fit, coef=3)
topTags(qlf.2vs1)

FDR <- p.adjust(qlf.2vs1$table$PValue, method="BH")

summary(decideTests(qlf.2vs1, p.value=0.01, lfc=1)) 
results_em48h <- topTags(qlf.2vs1, n=20000, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)$table
results_em48h <-(results_em48h[which(results_em48h$PValue < 0.01 & (results_em48h$logFC >1 | results_em48h$logFC < -1)), ])
write.table(results_em48h, file= 'results_em48h.txt',sep='\t')

#differential expression analysis - KO vs control
fit <- glmQLFit(y, design)
qlf.2vs1 <- glmQLFTest(fit, coef=4) 

FDR <- p.adjust(qlf.2vs1$table$PValue, method="BH")

summary(decideTests(qlf.2vs1, p.value=0.01, lfc=1)) 
results_ko <- topTags(qlf.2vs1, n=20000, adjust.method = "BH", sort.by = "PValue", p.value = 0.01)$table
results_ko <-(results_ko[which(results_ko$FDR < 0.01 & (results_ko$logFC >1 | results_ko$logFC < -1)), ])
write.table(results, file= 'results_ko.txt',sep='\t')


#enrichment analysis - em24h vs control
em24h <- rownames(results_em24h)

GO_bp_em24h <-enrichGO(em24h, keyType = "ENSEMBL", 'org.Hs.eg.db', ont="BP",
                     pvalueCutoff=0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH", universe = matrix$gene_id)

dotplot(GO_bp_em24h, showCategory=20, title="EM127 24h vs WT NT", font.size = 10)

results_GO_em24h <- GO_bp_em24h@result[which(GO_bp_em24h@result$p.adjust < 0.05),]
write.table(results_GO_em24h, file="enrichment_em24h.txt")


#enrichment analysis - KO vs control
ko <- rownames(results_ko)

GO_bp_ko <-enrichGO(ko, keyType = "ENSEMBL", 'org.Hs.eg.db', ont="BP",
                       pvalueCutoff=0.05, qvalueCutoff = 0.05, pAdjustMethod = "BH", universe = matrix$gene_id)

dotplot(GO_bp_ko, showCategory=20, title="SMYD KO vs WT NT", font.size = 10)

results_GO_ko <- GO_bp_ko@result[which(GO_bp_ko@result$p.adjust < 0.05),]
write.table(results_GO_ko, file="enrichment_KO.txt")
