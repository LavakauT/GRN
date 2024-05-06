## ∆ CONSENSUS motif comparison------------
library(dplyr)
library(tibble)
library(universalmotif)
library(stringr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(gridtext)
library(pgirmess)

# library(ggplot2)

dir <- '/Volumes/R425/RNA-seq/arabidopsis/similarity/DD/'
file_name <- 'DD'
filenames <- list.files('/Volumes/R425/RNA-seq/arabidopsis/similarity/DD',
                        pattern="*.txt",
                        full.names=FALSE)
head(filenames)

df <- read.delim(paste0(dir, filenames))
names(df) <- c('motif')

list.motif <- list()
for (i in 1:nrow(df)) {
  m <- create_motif(df[i,], name = df[i,])
  list.motif <- c(list.motif, m)
  list.pwm <- merge_motifs(list.motif, method = "PCC")
  list.pwm@name <- file_name
}


# DAP (Arabidopsis)
dap <- read_meme('/Volumes/R425/RNA-seq/arabidopsis/similarity/ArabidopsisDAPv1.meme')

comparisons <- compare_motifs(c(list.motif, dap),
                              method = "PCC",
                              min.mean.ic = 0,
                              score.strat = "a.mean")

# remove pCRE to pCRE PCC
comparisons2 <- comparisons[-c(1:nrow(df)), ]


for (i in 1:nrow(df)) {
  if(i == 1){
    
    df3 <- data.frame(comparisons2[,i])
    cname <- list.motif[[i]]@name
    colnames(df3) <- cname
    rname <- row.names(df3)
    df3$motif <- rname
    df3 <- df3[order(-df3[cname]), ][1,] # top similar motifs
    df3 <- df3[,c('motif', cname)]
    row.names(df3) <- 1
    
  } else {
    
    df2 <- data.frame(comparisons2[,i])
    cname2 <- list.motif[[i]]@name
    colnames(df2) <- cname2
    rname2 <- row.names(df2)
    df2$motif <- rname2
    df2 <- df2[order(-df2[cname2]), ][1,]
    df2 <- df2[,c('motif', cname2)]
    row.names(df2) <- 1
    
    df3 <- full_join(df3, df2, by = 'motif')
  }
  
}



df3 <- gather(df3, key = 'kmers', value = 'PCC', -motif)
df3 <- df3 %>% filter(PCC != '')
df3$cluster <- file_name


# DAP motif family----
motif.family <- df3


string <- motif.family$motif
# remove dashline and any digit behind it as hamily name
string <- str_replace_all(string, '\\_m1', '')
head(string)

motif.family$motif <- string

motif.family <- motif.family %>%
  group_by(kmers, motif) %>%
  distinct() %>% 
  summarise_all(median)

motif.family <- motif.family %>% 
  column_to_rownames(var = 'kmers')


# DAP only
# motif.family$motif <- str_replace_all(motif.family$motif, '\\_ecoli\\.', '\\_tnt\\.')
tnt <- str_locate_all(motif.family$motif, '\\_tnt\\.')


for (i in 1:length(tnt)) {
  if (i == 1){
    tnt.1 <- data.frame(tnt[[i]])
  } else{
    tnt.2 <- data.frame(tnt[[i]])
    tnt.1 <- rbind(tnt.1, tnt.2)
  }
}

tnt.1$end <- (tnt.1$start)-1
tnt.1$start <- 1

motif.family$motif <- str_sub(motif.family$motif,
                              start = tnt.1$start,
                              end = tnt.1$end)

output <- data.frame(motif.family)
output <- rownames_to_column(output, var = 'k-mers')

write.delim(output, paste0(dir, file_name, '_top_sim_dap.txt'))

# DAP motif family heatmap-------

col_fun = colorRamp2(c(0, 0.5, 1), c('blue', 'white', 'red'))
col_fun2 = colorRamp2(c(0, 0.5, 1), c('white', 'white', 'white'))
motif.family2 <- motif.family %>% filter(motif == 'HSF')
ha = rowAnnotation(PCC = as.numeric(motif.family2$PCC),
                   col = list(PCC = col_fun),
                   annotation_legend_param = list(
                     PCC = list(direction = "horizontal")
                   ))


mfh <- data.frame(motif.family2[,1])
colnames(mfh) <- 'motif'
row.names(mfh) <- rownames(motif.family2)


ht <- Heatmap(mfh,
              name = 'family',
              col = col_fun2,
              column_title = gt_render(
                paste0("<span style = 'font-size:10pt'>**DAP-seq**</span><br>")),
              show_column_names = FALSE,
              row_title = 'Eenriched Kmers from promoter regions',
              row_names_gp = gpar(fontsize = 8),
              row_title_gp = gpar(fontsize = 8),
              show_row_names = TRUE,
              show_heatmap_legend = FALSE,
              row_order = rownames(motif.family2),
              right_annotation = ha,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf(mfh[i, j]), x, y, gp = gpar(fontsize = 10,
                                                              fontface = 'bold.italic'))
                grid.rect(gp = gpar(lwd = 2, fill = "transparent"))
              }
)


tiff(paste0(dir, file_name, '_top_sim_dap.tif'),
     width = 700, height = 700, res = 300)
#, type="cairo"
draw(ht,
     merge_legend = TRUE,
     heatmap_legend_side = "bottom",
     annotation_legend_side = 'bottom')

dev.off()

save.image(file = paste0(dir, 'dap.RData'))




