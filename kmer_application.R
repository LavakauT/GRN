##### load packages #####
library(dplyr)
library(tidyverse)
library(stringr)

#######################

####### kmer nodes edges application ######
path <- '/Volumes/R425/Lavakau/kmers_data/GRN/script/arabidopsis/deg'
setwd(paste(path, 'deg', sep = '/'))

sub.dir <- list.dirs(paste(path, 'motif', sep = '/'),
                     full.names = FALSE,
                     recursive = FALSE)

dir <- paste(path, 'motif', sep = '/')
dir2 <- paste(path, 'training/tss_5utr_promoter_pcc', sep = '/')

for (z in 1:length(sub.dir)) {
  attr <- read.csv('Merged Network default node.csv') %>% 
    rename(ID = name)
  kmers <- read.table(paste(dir, sub.dir[z], paste0(sub.dir[z], '_top_sim_dap.txt'), sep = '/'),
                      header = T) %>% 
    filter(PCC >= 0.9) %>% 
    filter(motif != '')
  imp <- read.table(paste(dir2, sub.dir[z], 'RF', 'imp.txt', sep = '/'), header = T) %>% 
    column_to_rownames(., var = 'file')
  
  
  imp <- data.frame(t(imp)) %>% 
    filter(count != '') %>% 
    select(-count)
  imp <- data.frame(t(imp))
  imp <- imp[,kmers$k.mers]
  
  
  dir3 <- paste0(paste(dir2, sub.dir[z], sep = '/'))
  tf <- unique(kmers$motif)
  
  columns2 = c("ID", "Class") 
  df.com2 = data.frame(matrix(nrow = 0, ncol = length(columns2))) 
  colnames(df.com2) = columns2
  df.com2$ID <- as.character()
  df.com2$Class <- as.character()
  
  for (i in 1:length(tf)) {
    kmer <- kmers %>% 
      filter(motif == tf[i]) %>% 
      pull(k.mers)
    
    columns = c("ID", "Class") 
    df.com = data.frame(matrix(nrow = 0, ncol = length(columns))) 
    colnames(df.com) = columns
    df.com$ID <- as.character()
    df.com$Class <- as.character()
    
    for (x in 1:length(kmer)) {
      col <- imp %>% 
        select(all_of(kmer[x]))
      col[is.na(col)] <- 0
      col$max <- apply(col, 1, max, na.rm=TRUE)
      col <- row.names(col[order(-col$max),][1,])
      col <- str_replace_all(col, '_RF_imp', '')
      
      df <- read.delim(paste(dir3, col, sep = '/'))
      df <- cbind(df[,1:2],
                  df %>%
                    select(all_of(kmer[x])))
      df <- df %>% filter(Class == 'pos')
      colnames(df) <- c('ID', 'Class', kmer[x])
      df.com <- full_join(df.com, df, by = c('ID', 'Class'))
    }
    
    df.com[is.na(df.com)] <- 0
    
    if(length(df.com) >3){
      df.com <- df.com %>% 
        mutate(tf = rowSums(df.com[,-c(1:2)]))
    } else{
      colnames(df.com) <- c('ID', 'Class', 'tf')
    }
    
    df.com <- df.com %>% 
      select(ID, Class, tf)
    colnames(df.com) <- c('ID', 'Class', paste(tf[i], sub.dir[z], sep = '_'))
    
    df.com2 <- full_join(df.com2, df.com, by = c('ID', 'Class'))
  }
  attr.sub <- left_join(attr, df.com2, by = 'ID')
  write.csv(attr.sub, paste0('Merged Network default node_', sub.dir[z], '.csv'))
}



attr.edge <- read.csv('Merged Network default edge.csv')
attr.edge$name <- str_replace_all(attr.edge$name, ' \\(interacts with\\) ', '_')
attr.edge <- separate(attr.edge, name, into = c('node_X', 'node_Y'), sep = '_')

attr.anno.x <- attr[c('ID', 'Family')] %>% 
  rename(node_X = ID,
         TF_X = Family)
attr.anno.y <- attr[c('ID', 'Family')] %>% 
  rename(node_Y = ID,
         TF_Y = Family)
attr.edge <- merge(attr.edge, attr.anno.x, by = 'node_X')
attr.edge <- merge(attr.edge, attr.anno.y, by = 'node_Y')
attr.edge$TF <- paste(attr.edge$TF_X, attr.edge$TF_Y, sep = '_')

write.csv(attr.edge,
          'Merged Network default edge2.csv',
          row.names = FALSE)

attr.edge$kmer <- 'Nodata'
for (i in 1:nrow(attr.edge)) {
  target <- attr.edge[i, 'node_Y']
  source.family <- attr.edge[i, 'TF_X']
  cluster <- attr %>% 
    filter(ID == target) %>% 
    pull(cluster)
  edge <- read.csv(paste0('Merged Network default node_', cluster, '.csv')) %>% 
    column_to_rownames(., var = 'ID')
  edge.colname <- colnames(edge)
  edge.colname <- str_remove_all(edge.colname, paste0('_', cluster))
  colnames(edge) <- edge.colname
  
  if(length(edge[target,source.family]) == 0){
    next
  }
  
  if(is.na(edge[target,source.family])){
    next
  }else if(edge[target,source.family] >= 1){
    attr.edge[i,]$kmer <- 'Support'
  }else if (edge[target,source.family] == 0){
    attr.edge[i,]$kmer <- 'Reject'
  }
  
}


write.csv(attr.edge,
          'Merged Network default edge_kmer.csv',
          row.names = FALSE)
#######################