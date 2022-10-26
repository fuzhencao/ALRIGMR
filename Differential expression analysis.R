df<-read.csv('GSE144269_voom.csv',header = T)#GSE144269 data
m=as.matrix(df[, 1])
v=as.vector(m)
X<-df[, -1]
row.names(X) <- v
X1normal<-X[,seq(1,140,2)]
X2cancer<-X[,seq(0,140,2)]
X<-cbind(X1normal,X2cancer)

library(limma)
group <- c(rep(c("a","b"),c(70,70)))#group information
design <- model.matrix(~0+factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(X)
contrast.matrix <- makeContrasts("a-b", levels = design)
fit <- lmFit(X,design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)
DEG_ot <- topTable(fit2, coef =1, n = Inf,p.value = 1)
plot(DEG_ot$logFC,-log10(DEG_ot$P.Value))
write.csv(DEG_ot, "DEG.csv",quote = FALSE)
#ÌôÑ¡²îÒì»ùÒò
choose<-subset(DEG_ot,DEG_ot$logFC>=1&DEG_ot$P.Value<=0.05 | DEG_ot$logFC<=-1&DEG_ot$P.Value<=0.05)
row1<-rownames(choose)
row2<-rownames(X)
shaixuangene<-match(row1,row2)
shaixuangene<-X[c(shaixuangene),]
write.csv(shaixuangene, "DEG_genes.csv",quote = FALSE)




