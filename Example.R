packageVersion("hgu133plus2hsentrezgcdf")
library(affy)

CELs.entrezgene = ReadAffy(cdfname="hgu133plus2hsentrezgcdf",
                           filenames=c('../data/GSM1154244_2Problema.CEL',
                                       '../data/GSM1154245_2Problema_2.CEL',
                                       '../data/GSM1154246_2Problema_3.CEL',
                                       '../data/GSM1154247_1Control.CEL', 
                                       '../data/GSM1154248_1Control_2.CEL',
                                       '../data/GSM1154249_1Control_3.CEL'))

eset.entrezgene = rma(CELs.entrezgene)
rma.exprs <- exprs(eset.entrezgene)

library(multtest)

groups=c(1, 1, 1, 0, 0, 0)
stats=mt.teststat(rma.exprs, groups, test="t")
rawp=2*(1-pnorm(abs(stats)))

# Permutation adjusted p-values for simple multiple testing procedures
procs<-c("Bonferroni","BH")
res2<-mt.rawp2adjp(rawp, procs)

# Plot results from all multiple testing procedures
allp <- cbind(res2$adjp[order(res2$index),])

sum(adjp < 0.05)