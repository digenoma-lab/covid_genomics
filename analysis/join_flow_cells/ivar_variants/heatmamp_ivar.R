knitr::opts_chunk$set(echo = TRUE)
library('tidyverse')
library('openxlsx')
library('ggbeeswarm')
library('patchwork')
library('ggplot2')
library('reshape2')
library('ComplexHeatmap')

#Ivar results
setwd('/Users/bbernal/Documents/covid_work/join_flow_cells/ivar_variants/')
ivar_files <- list.files(pattern = "\\.tsv$")
sample <- ivar_files[0:96]
ivar_results <- list()
for(i in sample){
  print(i)
  sample = paste('sample_', 
                 strsplit(strsplit(i, '_')[[1]][3],  '.', fixed = TRUE)[[1]][1], 
                 sep = '')
  variant <- substring(sample, 1, nchar(sample)-3)
  df <- read.table(i, h=TRUE,sep="\t")
  ivar_results[[variant]] <- df$POS}


#generation of binary matrix 
mut_matrix <- t(splitstackshape:::charMat(listOfValues = ivar_results, fill = 0L))
colnames(mut_matrix) <- names(ivar_results)
mut_matrix <- mut_matrix[order(rowSums(mut_matrix),decreasing=T),] #sort by position with more mutations
mut_matrix <- t(mut_matrix) #transpose

mut_matrix <- mut_matrix[-71,]

#counting how many mutations by Ivar has every sample
n_mut_ivar <- data.frame(rowSums(mut_matrix))
n_mut_ivar <- cbind(newColName = rownames(n_mut_ivar), n_mut_ivar)
rownames(n_mut_ivar) <- 1:nrow(n_mut_ivar)
colnames(n_mut_ivar) <- c('sample','num_mut')

#Demographic samples info
demographic_info <- read.table('/Users/bbernal/Documents/covid_work/demographic.csv', h=TRUE,sep=",")
demographic_info = demographic_info[c('sample', 'Comuna', 'Fecha.Análisis', 'Edad')]
demographic_info <- demographic_info %>% 
  mutate(ciudad = case_when(Comuna == 'RANCAGUA'  ~ 'Rancagua' ,
                            .default = 'Santiago'))

demographic_info <- demographic_info[order(match(demographic_info[,1],rownames(mut_matrix))),]

#freya results
freya_results <-read.table('/Users/bbernal/Documents/covid_work/freya_results.txt', h=TRUE,sep=",")
freya_variants <-read.table('/Users/bbernal/Documents/covid_work/freya_variants.tsv', h=TRUE,sep="\t")
freya_variants <- freya_variants[,order(colnames(freya_variants))]
freya_variants <- melt(freya_variants)
#getting variant more prevalent
max <- freya_variants %>% group_by(variable) %>% top_n(1, value) 

#deleting sample 71 
n_mut_ivar<-n_mut_ivar[!grepl("71", 
                              n_mut_ivar$sample),]
demographic_info<-demographic_info[!grepl("71", 
                          demographic_info$sample),]
#adding freya results with sampleid names from ivar (normalizing sample names)
freya_results$sampleId <- n_mut_ivar$sample
#melt to plot 
freya_results <- melt(freya_results, id="sampleId")
freya_results <- filter(freya_results, variable == "pass_filter")

#joining all results
all_results <- cbind(freya_results, mut_ivar = n_mut_ivar$num_mut)
all_results$variant <- max$X
all_results$location <- demographic_info$ciudad
all_results$age <- demographic_info$Edad
all_results$analisis_data <- demographic_info$Fecha.Análisis


#ComplexHeatMap
h_cities= rowAnnotation(location = all_results$location,
                        variant = all_results$variant,
                        col = list(location = c("Rancagua" =  "red", 
                                            "Santiago" = "blue")),
                        gp = gpar(fontsize = 3),
                        border = TRUE)
ht1 <- Heatmap(mut_matrix, 
               name = "iVar mutation",
               rect_gp = gpar(col = "gray70", lwd = 0.1),
               row_names_gp = gpar(fontsize = 7, col ='gray9'),
               column_names_gp = gpar(fontsize = 5, col ='gray9'), # Text size for row names
               col = c("royalblue3", "red2"),
               row_km = 2,
               right_annotation = h_cities)

ht1 


ggplot(data=all_results, aes(x=mut_ivar, y=value, label=sampleId, 
                             colour = variant, shape=location)) + 
  geom_point() +
  #geom_text(hjust=0, vjust=0, size=2) + 
  theme_minimal()


#samtools_depth
setwd('/Users/bbernal/Documents/covid_work/join_flow_cells/samtools_depth/')
samtools_files <- list.files(pattern = "\\.txt$")
sample <- samtools_files[0:96]

depth <- list()
for(i in sample){
  print(i)
  #transform name of sample to sample_NN
  sample = paste('sample_', 
                 strsplit(strsplit(i, '_')[[1]][3],  '.', fixed = TRUE)[[1]][1], 
                 sep = '')
  variant <- substring(sample, 1, nchar(sample)-3)
  df <- read.table(i, h=FALSE,sep="\t")
  colnames(df) <- c('ref', 'pos', 'depth')
  depth[[variant]] <- df}

samtools_results <- lapply(depth, as.data.frame) %>% 
  bind_rows(.id = "sampleId")
samtools_results <- left_join(samtools_results,all_results[c('sampleId', 'location', 'variant')], by="sampleId")

samtools_results<-samtools_results[!grepl("71", 
                                          samtools_results$sample),]

grid1 <- ggplot(samtools_results, aes(x = pos, y=log(depth+1),  label=sampleId, 
                                     colour = variant)) +
  geom_line(size=0.25, alpha=0.7) +
  geom_hline(yintercept = log(100),color="gray20",linetype=1,size=0.5)+
  #geom_text(aes(x=2000, label="100X Coverage", y=log(100)), colour="gray20", vjust = 1.2)+
  facet_grid(cols = vars(location), rows = vars(variant)) +
  #facet_wrap(location~sampleId, scales = 'free_y', nrow=12) +
  #scale_y_continuous(trans='log10',labels = scales::comma) +
  theme_minimal()

grid1


#variats bar plot
p1 <- ggplot(data=freya_variants, aes(x=variable, y=value, fill=X)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(angle=90,hjust=1))
p1