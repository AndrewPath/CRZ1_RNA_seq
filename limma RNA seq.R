install.packages("gplots")
library(gplots)
#take the countadat from deseq script(as a matrix)
#filter per rownames to keep out low reads count
#countdatafilt <-countdata %>% filter_all(all_vars(.> 300))
colnames(countdata) <- c("wt4b","wt7b","wt10b","wt14b", "wt4c", "wt7c","wt10c","wt14c","d4b","d7b","d10b","d14b", "d4c", "d7c","d10c","d14c" )
countdata <- read.table("/Users/andreacacciotti/Desktop/countdata.txt" )

#in logCPM prior count is used here to damp down the variances of logarithms of low counts.
#logCPM <- cpm(dge, log=TRUE, prior.count=3)
#Create DGEList object

dge <- DGEList(counts=countdata)
#Calculate normalization factors

dge <- calcNormFactors(dge)
#Filter low-expressed genes

cutoff <- 1
drop <- which(apply(cpm(dge), 1, max) < cutoff)
d <- dge[-drop,] 
dim(d) # number of genes left

snames <- colnames(countdata) 
snames
genotype <- c("w","w","w","w","w","w","w","w","d","d","d","d","d","d","d","d") 
time <- c("4","7","10","14","4","7","10","14","4","7","10","14","4","7","10","14")
group <- interaction(genotype, time)
group

mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
y
#lmFit fits a linear model using weighted least squares for each gene:

fit <- lmFit(y, mm)
head(coef(fit))

#Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models:

contr <- makeContrasts(groupw.14 - groupd.14, levels = colnames(coef(fit)))
contr

#Estimate contrast for each gene

tmp <- contrasts.fit(fit, contr)

#shrink the error

tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))

#filter and take DEG in new df
countdatafilt <-top.table%>%filter(top.table$adj.P.Val < 0.05)

#take the top gene of the delta (<) and wt (>)
#topexprd14 <- countdatafilt%>%filter(countdatafilt$logFC < 0)
#take the matrix rownames as list(for each condition) 
#degd14 <- rownames(topexprd14)
#x<-list(degwt7,degwt14,degd7,degd14) 
#saveRDS(x, file="listconf.RData")

#per venn diagram metti i confronti direttamente tra cond e temp
#filter and take DEG in new df
countdatafilt <-top.table%>%filter(top.table$adj.P.Val < 0.05)
#take the matrix rownames as list(for each condition) 
deg4 <- rownames(countdatafilt)
x<-list(deg4,deg7,deg10,deg14) 
names(x) <- c("4","7","10","14")
saveRDS(x, file="listconf.RData")

#volcano plot DEG
#transform matrix in dataframe 
df1 <- top.table[3:4]<- NULL
library(tibble)
df <- tibble::rownames_to_column(top.table, "gene")
df$diffexpressed <- "NO"
#if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
df$diffexpressed[df$logFC > 0.6 & df$P.Value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$logFC < -0.6 & df$P.Value < 0.05] <- "DOWN"
p <- ggplot(data=df, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + geom_point() + theme_minimal() 
print(p)
p + geom_vline(xintercept=c(-0.6, 0.6), col="red")

#heatmap with  the expression value (corretto)
logCPM <- cpm(d, log=TRUE, prior.count=1)
#trasforma matrix logcpm in dataframe
logCPMdf <- as.data.frame(logCPM)
#hm
wt14.vs.d14.topgenes <- df$gene[1:100]
logCPMdf$name <- rownames(logCPMdf)
i <- which(logCPMdf$name %in% wt14.vs.d14.topgenes)
<<<<<<< HEAD
mycol <- colorpanel(1000,"red","black","green")
heatmap.2(logCPM[i,], scale="row",
          labRow=logCPMdf$name[i], labCol=group, 
          col=viridis(n = 64,option = "A",direction = 1), trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")
dev.off()
=======
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(logCPM[i,], scale="row",
          labRow=logCPMdf$name[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

>>>>>>> eccf4c55f9e610bf3bb8d4a0158f37bf1714c5e8

#prove per hm 
#hm tot delle 4 condizioni
#Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models:

contr <- makeContrasts(groupw.14 - groupd.14, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))

countdatafilt14 <-top.table%>%filter(top.table$adj.P.Val < 0.05)
#unire tutti i dataframe
totDGE<-rbind(countdatafilt4,countdatafilt7,countdatafilt10,countdatafilt14)
#sort per significatività
totDGE <- totDGE[order(totDGE$adj.P.Val),]

#transform matrix in dataframe and put rownames as a variable
library(tibble)
df1 <- tibble::rownames_to_column(totDGE, "gene")
#togliamo i 4 e 7 dai test
test<-logCPMdf
test <- test[,-9]
test <- mutate_all(test, function(x) as.numeric(as.character(x)))
test <- as.matrix(test, header=T)
is.numeric(test)

#hm top gene senza 4 e 7 (not work) 
wt14.vs.d14.topgenes <- df1$gene[1:500]
logCPMdf$name <- rownames(logCPMdf)
i <- which(logCPMdf$name %in% wt14.vs.d14.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(test[i,], scale="row",
          labRow=logCPMdf$name[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          dendrogram="column")

#togliamo i 4 e 7 dai test
test<-logCPMdf
test <- test[,-9]
test <- mutate_all(test, function(x) as.numeric(as.character(x)))
test <- as.matrix(test, header=T)
is.numeric(test)


fit <- lmFit(y, mm)
head(coef(fit))

#Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models:
#contrast for each condition, change the number after the point
contr <- makeContrasts(groupw.7 - groupd.7, levels = colnames(coef(fit)))
contr

#Estimate contrast for each gene

tmp <- contrasts.fit(fit, contr)

#shrink the error

tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
length(which(top.table$adj.P.Val < 0.05))

#filter and take DEG in new df
countdatafilt <-top.table%>%filter(top.table$adj.P.Val < 0.05)
#take the top gene of the delta (<) and wt (>),for each condition create new matrix
topexprd7 <- countdatafilt%>%filter(countdatafilt$logFC < 0)
#take the matrix rownames as list(for each condition),
degd7 <- rownames(topexprd7)
#create the x list with the 4 condition for venn
x<-list(degwt7,degwt14,degd7,degd14) 

#save the list and pass to mac                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          

saveRDS(x, file="fname.RData")

#Load with

readRDS("fname.RData")
#rename factor inside  with 
names(x) <- c("WT7","WT14","D7","D14")

install.packages("ggVennDiagram")
library(ggVennDiagram)