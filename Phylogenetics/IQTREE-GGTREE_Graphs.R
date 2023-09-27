install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("muscle")
BiocManager::install("Biostrings")
library('ggtree')
library('muscle')
library('Biostrings')
library('treeio')
library('tidyverse')

getwd()
my_tree_file_path <- "../Phylogenetics/iqtree2.treefile"
my_tree <- read.newick(my_tree_file_path)

#transforming tree to add variant information
x <- as_tibble(my_tree)

all_results <-read.table('../analysis/all_results.csv', h=TRUE,sep=",")
all_results <- as_tibble(all_results)
all_results
#formatting sample name
x$sampleId <- paste0('sample_',(str_extract(x$label,"(?<=-)\\d+")))
head(x)
#keeping name of reference
x <- x %>% 
  mutate(sampleId = ifelse(label == 'MN908947',label,sampleId))
#joining data
x_new <- treeio::full_join(x,all_results[c('sampleId', 'location', 'variant')], by="sampleId")

#transforming data into a tree object again!
x_new <- as.treedata(x_new)

### MUSCLE ANALYSIS - RUN ONE TIME ###
#extracting fasta files
#folder = "covid_genomics/Phylogenetics/consensus/"
#files = list.files(folder, pattern = ".fa")
#fa_seq = lapply(files,readDNAStringSet)
#fa_seq = do.call(c,fa_seq)
#tipseq_aln <- muscle::muscle(fa_seq) 
#saving file
#writeXStringSet(tipseq_aln, filepath = 'covid_genomics/Phylogenetics/samples_alignment.fasta', format = 'fasta' )

#reading alignment file from MUSCLE
tipseq_aln <- readDNAStringSet("../Phylogenetics/samples_alignment.fasta")
tipseq_aln <- DNAStringSet(tipseq_aln)
#calculating distance between sample
tipseq_dist <- stringDist(tipseq_aln, method = "hamming")
## calculate the percentage of differences
tipseq_d <- as.matrix(tipseq_dist) / width(tipseq_aln[1]) * 100

## convert the matrix to a tidy data frame for facet_plot
dd <- as_tibble(tipseq_d)
dd$seq1 <- rownames(tipseq_d)
td <- gather(dd,seq2, dist, -seq1)
g <- p1$data$variant
names(g) <- p1$data$label
td$clade <- g[td$seq2] 

#this sample is deleted due to low coverage 
td<-td[!grepl("03-03", 
              c(td$seq1)),]
td<-td[!grepl("03-03", 
              c(td$seq2)),]

#adding a new column to set further order of apperance in the graph
td <- td %>%
  mutate(
    alpha = case_when(
      str_detect(clade, "Omicron") ~ "Omicron",
      .default = clade
    )
  )

#setting colors with the same palette 
variant <- setNames(RColorBrewer::brewer.pal(name = "Dark2", n = length(unique(all_results$variant))), unique(all_results$variant))
#sort to make alpha based on appearence 
td<- td[order(td$alpha), decreasing=TRUE]

#plot tree
p1 <- ggtree(x_new, color="grey40", ladderize = TRUE) +
  xlim_tree(0.00325) +
  geom_tiplab(aes(color=variant, label = sampleId), size=3.) +
  scale_color_manual('', values = variant) +
  guides(color = guide_legend(override.aes = list(size = 5, label = "\u25AA")))+
  theme_tree2(legend.position = "bottom")

p1

#adding sequence distnace to phylo tree plot
p2 <- p1 + geom_facet(panel = "Sequence Distance", 
                     data = td, geom = geom_path,  size=0.35,
                     mapping=aes(x = dist, group = seq2, 
                                 color = clade, alpha = alpha),
                     show.legend = F) 

p2
#reshaping proportion of each plot
facet_labeller(p2, c(Tree = "phylogeny")) %>% facet_widths(c(Tree = 1.))
ggsave("PhylogenyIQTREE_SeqDist.png", dpi=300, dev='png', height=25, width=20, units="cm")

p1 <- ggtree(new_tree, color="blue") + geom_tiplab(size=1.7, color='blue') + geom_treescale(fontsize=4, linesize=1, offset=1, color='blue')

ggsave("PhylogenyIQTREE_1.png")

p2 <- ggtree(my_tree, color="blue") + geom_tiplab(size=1.7, color='blue')
p2
ggsave("PhylogenyIQTREE_2.png")

p2 <- ggtree(my_tree, color="blue") + geom_tiplab(size=1.2, color='blue')
ggsave("PhylogenyIQTREE_3.png")

p2 <- ggtree(my_tree, color="blue") + geom_tiplab(size=1.0, color='blue')
ggsave("PhylogenyIQTREE_4.png")

p2 <- ggtree(my_tree, color="blue") + geom_tiplab(size=0.8, color='blue')
ggsave("PhylogenyIQTREE_5.png")

p2 <- ggtree(my_tree, color="blue") + geom_tiplab(size=0.6, color='blue')
ggsave("PhylogenyIQTREE_6.png")

p2 <- ggtree(my_tree, color="blue") + geom_tiplab(size=0.5, color='blue')
ggsave("PhylogenyIQTREE_7.png")


my_tree_file_path <- "iqtree2-Copy.treefile"
my_tree <- read.newick(my_tree_file_path)

p3 <- ggtree(my_tree, color="blue") + geom_tiplab(size=2.0, color='blue') + geom_treescale(fontsize=4, linesize=1, offset=1, color='blue')
ggsave("PhylogenyIQTREE_8.png")

p3 <- ggtree(my_tree, color="blue") + geom_tiplab(size=2.0, color='blue') + theme_tree2()
ggsave("PhylogenyIQTREE_9.png")

p4 <- ggtree(my_tree, color="blue", branch.length='none', layout='circular') + geom_tiplab(size=2.0, color='blue') 
p4
ggsave("PhylogenyIQTREE_Circular.png")
