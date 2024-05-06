##### load packages #####
library(DESeq2)
library(tidyverse)
library(airway)
library(sva)
library(devtools)
library(WGCNA) # may ask install impure and preprocessCore packages
library(ComplexHeatmap)
library(circlize)
library(gridtext)
library(factoextra)

devtools::install_github("zhangyuqing/sva-devel")

##### DEGs with DESeq2 #####
# The example script used Arabidopsis heat stress dataset from 5 individaul publishes.
# There are sample metadata file (coldata.txt) and raw count data file (counts.txt) we should prepare.
# Follow this part, we would get DEG list from each cluster.
# We can also send the DEG list to PRIDICT.com for getting kmer list in each cluster.

###### load data ######
path <- '/Volumes/R425/Lavakau/kmers_data/GRN/script/arabidopsis'
setwd(path)
coldata <- read.delim('coldata.txt') # read sample metadata
counts.data <- read.delim('counts.txt') # read raw count
#######################

###### preprocess ######
rname <- counts.data[,1] # gene name to the rawname
counts.data <- counts.data[,7:ncol(counts.data)]
row.names(counts.data) <- rname

colname <- coldata[,1]
colnames(counts.data) <- colname
head(counts.data)

rname2 <- coldata[,1] # sample name to raw name
coldata <- coldata[,2:ncol(coldata)]
row.names(coldata) <- rname2
coldata$genotype <- 'col'
coldata$temperature <- paste(coldata$treatment, coldata$temperature, sep = '_')
head(coldata)
#######################

###### ComBat-seq batch adjusting ######
# batch effect adjusting with ComBat-seq
# first, check raw data quality
gsg <- goodSamplesGenes(t(counts.data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes) # 1031/26385 False/True as good genes
table(gsg$goodSamples)

data <- counts.data[gsg$goodGenes == TRUE,]
htree <- hclust(dist(t(data)), method = "average")
plot(htree) # check the outlier groups in your input data

# Accordign to the htree, s1 to s6 are outliers in input data.
# Aherefore, we remove s1:s6 to keep the credibility of downstreaming analysis. 
coldata <- coldata[7:28,]
data <- data[,7:28]


batch <- coldata[, 4] # assign batch infromation
count_data_mat <- as.matrix(data) 
covar_mat <- data.frame(coldata[,c(1,3)]) # assign multiple covariates, here are treatment (control/heat) and treatment time (00/30/60 minutes)
adjusted.count.data <- ComBat_seq(count_data_mat, batch = batch, covar_mod = covar_mat)


# check the sample name in the column of adjusted count data are the same with 
# sample name in the row of coldata
all(colnames(adjusted.count.data) %in% row.names(coldata))
all(colnames(adjusted.count.data) == row.names(coldata))


# force time column as factor
coldata$time <- factor(coldata$time)
#######################

###### DESeq ######
dds <- DESeqDataSetFromMatrix(countData = adjusted.count.data,
                              colData = coldata,
                              design = ~ 1)
dds

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

ddsMF <- dds
ddsMF$group <- factor(paste(ddsMF$treatment, ddsMF$time, sep = '_'))
design(ddsMF) <- ~ group
head(ddsMF)

ddsMF$group
ddsMF$group <- relevel(ddsMF$group, ref = "CT_0")
ddsMF <- DESeq(ddsMF)
#######################

####### quality control ######
#transform counts data into newly DESeqDataSets
vsdMF <- vst(ddsMF, blind = FALSE)
head(assay(vsdMF), 5) #assay is to extract the matrux of normalized values



plotPCA(vsdMF, intgroup = c("treatment", "time"))  

pcaData <- plotPCA(vsdMF, intgroup = c("treatment", "time"), returnData = TRUE)
percenVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = as.factor(time), shape = treatment)) + 
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percenVar[1], "% variance")) +
  ylab(paste0("PC2: ", percenVar[2], "% variance")) +
  labs(color = 'treatment time (minutes)') +
  coord_fixed()

# exprot vst-transformed count data
setwd(path2)
write.table(assay(vsdMF), 'vst_norm_count.txt', quote = F, sep = '\t')
#######################

###### export DEG lists in each cluster ######
resultsNames(ddsMF)
range <- data.frame(group = c('group'),
                    treatment = c('HS_30', 'HS_60'),
                    base = 'CT_0')
dir.create(paste(path, 'deg', sep = '/'))
dir <- paste(path, 'deg', sep = '/')


for (i in 1:nrow(range)) {
  res <- results(ddsMF, contrast = c(range$group[i], range$treatment[i], range$base[i])) 
  res0.05 <- results(ddsMF, alpha = 0.05, contrast = c(range$group[i], range$treatment[i], range$base[i]))
  summary(res)
  summary(res0.05)
  
  write.csv(as.data.frame(res0.05), file = paste0(dir, '/', range$treatment[i], '_', range$base[i], '.csv'))
  res0.05 <- read.csv(paste0(dir, '/', range$treatment[i], '_', range$base[i], '.csv'))
  res0.05$Up_regulated <- 'NO' 
  res0.05$Up_regulated[res0.05$log2FoldChange >= 1.0 & res0.05$padj < 0.05] <- 'YES'
  res0.05$Down_regulated <- 'NO'
  res0.05$Down_regulated[res0.05$log2FoldChange <= -1.0 & res0.05$padj < 0.05] <- 'YES'
  res0.05$Negative <- 'NO' 
  res0.05$Negative[res0.05$log2FoldChange > -0.8 & res0.05$padj > 0.05] <- 'YES'
  res0.05$Negative[res0.05$log2FoldChange < 0.8 & res0.05$padj > 0.05] <- 'YES'
  write.csv(res0.05, paste0(dir, '/', range$treatment[i], '_', range$base[i], '.csv'), row.names = FALSE)
  
  
  # output gene list
  library(pgirmess)
  res0.05 <- read.csv(paste0(dir, '/', range$treatment[i], '_', range$base[i], '.csv'))
  
  up <- res0.05 %>% filter(Up_regulated =='YES') %>% select(X)
  write.table(up$X, file = paste0(dir, '/', 'ara_', range$treatment[i], '_', range$base[i], '_u', '.txt'),
              row.names = F, quote = F, sep = '\t')
  
  dwn <- res0.05 %>% filter(Down_regulated =='YES') %>% select(X)
  write.table(dwn$X, file = paste0(dir, '/', 'ara_', range$treatment[i], '_', range$base[i], '_d', '.txt'),
              row.names = F, quote = F, sep = '\t')
  
  neg <- res0.05 %>% filter(Negative =='YES') %>% select(X)
  write.table(neg$X, file = paste0(dir, '/', 'ara_', range$treatment[i], '_', range$base[i], '_n', '.txt'),
              row.names = F, quote = F, sep = '\t')
  
}
#######################

###### use Upset plot to determine cluster ###### 
path2 <- paste(getwd(), 'deg', sep = '/')
u.30 <- read.table(paste(path2, 'ara_HS_30_CT_0_u.txt', sep = '/'), header = T)[,1]
d.30 <- read.table(paste(path2, 'ara_HS_30_CT_0_d.txt', sep = '/'), header = T)[,1]
n.30 <- read.table(paste(path2, 'ara_HS_30_CT_0_n.txt', sep = '/'), header = T)[,1]

u.60 <- read.table(paste(path2, 'ara_HS_60_CT_0_u.txt', sep = '/'), header = T)[,1]
d.60 <- read.table(paste(path2, 'ara_HS_60_CT_0_d.txt', sep = '/'), header = T)[,1]
n.60 <- read.table(paste(path2, 'ara_HS_60_CT_0_n.txt', sep = '/'), header = T)[,1]

lt <- list(u.30, d.30, n.30, u.60, d.60, n.60)
list_to_matrix(lt)

names(lt) <- c('Up_0.5HR', 'Dwn_0.5HR', 'Neg_0.5HR', 'Up_1HR', 'Dwn_1HR', 'Neg_1HR')
m1 <- make_comb_mat(lt)
comb_name(m1)
comb_size(m1)


UpSet(t(m1[comb_degree(m1) == 2]))

ss = set_size(m1)
cs = comb_size(m1)
ht = UpSet(m1, 
           set_order = order(ss),
           comb_order = order(comb_degree(m1), -cs),
           top_annotation = HeatmapAnnotation(
             "Gene Intersections" = anno_barplot(cs, 
                                                 ylim = c(0, max(cs)*1.1),
                                                 border = FALSE, 
                                                 gp = gpar(fill = "black"), 
                                                 height = unit(4, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           left_annotation = rowAnnotation(
             "Gene Per Dataset" = anno_barplot(-ss, 
                                               baseline = 0,
                                               axis_param = list(
                                                 at = c(0, -1500, -5000),
                                                 labels = c(0, 1500, 5000),
                                                 labels_rot = 0),
                                               border = FALSE, 
                                               gp = gpar(fill = "black"), 
                                               width = unit(4, "cm")
             ),
             set_name = anno_text(set_name(m1), 
                                  location = 0.5, 
                                  just = "center",
                                  width = max_text_width(set_name(m1)) + unit(4, "mm"))
           ), 
           right_annotation = NULL,
           show_row_names = FALSE)


ht = draw(ht)
od = column_order(ht)
decorate_annotation("Gene Intersections", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
})

comb_degree(m1)
UU <- extract_comb(m1, "100100")
DD <- extract_comb(m1, "010010")
NN <- extract_comb(m1, "001001")


ma <- t(data.frame(UU = c('1','1'),
                   DD = c('-1', '-1')))
colnames(ma) <- c('0.5HR', '1HR')


count <- c('1298','1294')

colors = structure(c('blue','grey60','red'), names = c("-1", "0", "1"))
ha = rowAnnotation('# of genes' = anno_text(count, location = 0.5, just = "center",
                                            gp = gpar(fill = rep(2:4, each = 4),
                                                      col = "black", border = "black", fill = "orange"),
                                            width = max_text_width(count)*1.2))

Heatmap(ma,
        name = 'log2(FC)',
        col = colors,
        show_row_dend =FALSE,
        show_row_names = TRUE,
        row_names_side = 'left',
        cluster_columns = FALSE,
        column_names_side = 'top',
        column_names_gp = gpar(fontsize = 15),
        column_names_rot = 0,
        row_title = 'HS responsive clusters',
        right_annotation = ha,
        heatmap_legend_param = list(direction = "horizontal",
                                    nrow = 1),
        column_split = c(1,2)) # check the clusters


setwd(path2)
write.table(UU, 'UU.txt', col.names = F, row.names = F, quote = F, sep = '\t')
write.table(DD, 'DD.txt', col.names = F, row.names = F, quote = F, sep = '\t')
write.table(NN, 'NN.txt', col.names = F, row.names = F, quote = F, sep = '\t')
#######################

###### attribution file for network construction ######
attr <- rbind(data.frame(gene = UU, cluster = 'UU'),
              data.frame(gene = DD, cluster = 'DD'))
setwd(path2)
write.table(attr, 'attr.txt', row.names = F,quote = F, sep = '\t')
#######################