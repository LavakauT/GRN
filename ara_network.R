##### load packages #####
library(dplyr)
library(tibble)
library(ggplot2)
library(pgirmess)
library(randomForest)
library(GENIE3)
#######################

###### load corresponding data ######
path <- '~/absolute/path/to/GRN-main'
setwd(path)
attr <- read.table('attr.txt', header = T)
norm.counts <- read.delim('file/vst_norm_count.txt')

# plantTFDB database for transcription factors annotation:
# https://planttfdb.gao-lab.org/download.php#tfext
at.tf <- read.delim('Ath_TF_list.txt') %>%  # protein, Gene_ID ,and TF-family
  distinct(.keep_all = TRUE) %>%
  rename(protein = TF_ID,
         gene = Gene_ID)
  
# TAIR database for general genes annotation:
# https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release
at.de <- read.delim('TAIR10_functional_descriptions.txt',
                    header = T) %>% 
  select(1,3) %>% # Gene_ID and short description
  distinct() %>% 
  rename(protein = Model_name)

at <- left_join(at.de, at.tf, by = 'protein')
attr.anno <- left_join(attr, at %>% 
                         select(gene, Family, Short_description) %>% 
                         distinct(), by = 'gene')
# check duplication deu to the similar TF families
attr.anno$gene %>% duplicated() %>% table()
# keep first match and remove the others
attr.anno <- attr.anno[!duplicated(attr.anno$gene),]

# exprot DEG attribution with more annotation
path3 <- paste(path, 'network', sep = '/')
dir.create(path2)
setwd(path2)
write.table(attr.anno, 'attr.anno.txt', row.names = F, quote = F, sep = '\t')
#######################

###### GENIE3 prediction ######
# RESTRICT THE CANDIDATE REGULATORS (TF)
input.genes <- attr.anno %>% filter(Family != '') %>% pull(1) %>% unique()
expr.matrix <- as.matrix(norm.counts[attr.anno$gene,])

weight.matrix.tf <- GENIE3(expr.matrix, regulators = input.genes)

# Get all the regulatory links
link.list <- getLinkList(weight.matrix.tf, threshold = 0.01)
# dim(link.list)
sum.link.list <- data.frame(quantile(link.list$weight))

setwd(path2)
range <- data.frame(confidence = c('low', 'medium', 'high'),
                    quantile = c('25%', '50%', '75%'))
for (i in 1:3) {
  sub.link.list <- getLinkList(weight.matrix.tf, threshold = sum.link.list[range[i,]$quantile,])
  write.table(sub.link.list, paste0('ara_hs_genie3_tf_', range[i,]$confidence, '.txt'),
              row.names = F, quote = F, sep = '\t')
}
#######################

##### ARACNe-AP #####
matrix <- data.frame(expr.matrix) %>% 
  rownames_to_column(., var = 'gene')

# export the matrix file and TF list for ARACNe-AP
setwd(path2)
write.table(matrix, 'matrix.txt', row.names = F, quote = F, sep = '\t')

write.table(input.genes, 'tfs.txt',
            row.names = F, col.names = F, quote = F, sep = '\t')
# Please follow the aracne.txt to infer regulatory links
#######################



###### Input files to Cytoscape for merged network ######
#######################



###### double ranking (MI vs weight) ######
setwd(path2)
df <- read.csv('Merged Network default edge.csv') %>% 
  select(weight, MI)

df$rank1 <- rank(-df$weight)
df$rank2 <- rank(-df$MI)
df$rank_all <- apply(df[,3:4], 1, median, na.rm=T)


df <- df[order(df$rank_all),]


columns <- c('source', 'weight', 'MI')
filter <- data.frame(matrix(nrow = 0, ncol = length(columns))) 
colnames(filter) <-  columns

range <- c(round(nrow(df)*0.03),
           round(nrow(df)*0.05),
           round(nrow(df)*0.1),
           round(nrow(df)*0.15),
           round(nrow(df)*0.2))

for (i in 1:5) {
  i <- range[i]
  z <- paste0('top', i)
  df2 <- df[1:i,]
  weight <- min(df2$weight)
  MI <- min(df2$MI)
  df2 <- data.frame('source' = z,
                    'weight' = weight,
                    'MI' = MI)
  filter <- rbind(filter, df2)
}

write.table(filter, 'filter.txt',
            row.names = F, quote = F, sep = '\t')
#######################


###### Use the filter file to subset the top20 sub-network ######
#######################
