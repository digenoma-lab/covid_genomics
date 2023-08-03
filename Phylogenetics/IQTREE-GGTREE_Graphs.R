install.packages("BiocManager")
BiocManager::install("ggtree")

library(ggtree)

library(treeio)

install.packages("cowplot")
 
my_tree_file_path <- "iqtree2.treefile"
my_tree <- read.newick(my_tree_file_path)

p1 <- ggtree(my_tree, color="blue") + geom_tiplab(size=1.7, color='blue') + geom_treescale(fontsize=4, linesize=1, offset=1, color='blue')
ggsave("PhylogenyIQTREE_1.png")

p2 <- ggtree(my_tree, color="blue") + geom_tiplab(size=1.7, color='blue')
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
