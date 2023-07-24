knitr::opts_chunk$set(echo = TRUE)
library('tidyverse')
library('openxlsx')
library('ggbeeswarm')
library('patchwork')
library('ggplot2')

## Reads Mapped
samtools_jfc = read.table('join_flow_cells/multiqc-join-flowcell/multiqc_data/multiqc_samtools_stats.txt', header = TRUE)
samtools_pfc = read.table('per_flow_cells/multiqc-per-flowcell-cecret/multiqc_data/multiqc_samtools_stats.txt', header = TRUE)

#selecting data
reads_mapped_jfc <- samtools_jfc %>%
  dplyr::select(Sample, reads_mapped, raw_total_sequences, reads_mapped_percent)
reads_mapped_pfc <- samtools_pfc %>%
  dplyr::select(Sample, reads_mapped, raw_total_sequences, reads_mapped_percent)

reads_mapped_pfc<-reads_mapped_pfc[!grepl("97-97|98-98|99-99|100-100|101-101|102-102|103-103|104-104", 
                                          reads_mapped_pfc$Sample),]

#adding type of data
type = c(rep("join_flow_cell", nrow(reads_mapped_jfc)), rep('per_flow_cell',  nrow(reads_mapped_pfc)))
#join both table
reads_mapped=rbind(reads_mapped_jfc, reads_mapped_pfc)
#join type of sample
reads_mapped=cbind(reads_mapped,type)

#violin plot with dots 
p3 =ggplot(reads_mapped,aes(y=reads_mapped_percent,x=type,color=type,fill=type, label=Sample))+
  geom_violin() +
  #geom_point(position = position_jitter(seed = 2, width = 0.5), size = 1.5) +
  xlab("") +
  ylab("Mapped Reads (%)") +
  theme_minimal() +
  theme(legend.position = "none")
p3


#barplot num of reads
#raw_total_sequences = reads_mapped + reads_unmapped
p1<-ggplot(reads_mapped,aes(x=raw_total_sequences/100000,y=fct_reorder(Sample,raw_total_sequences),fill=type)) + 
  geom_col()+
  facet_wrap(~type,nrow=1,scales = "free_y")+
  geom_vline(data=filter(reads_mapped,type=="join_flow_cell"),aes(xintercept = 50))+
  geom_text(aes(x=50, label="50M target", y=10), colour="black", angle=90, vjust = 1.2,data=filter(reads_mapped,type=="join_flow_cell"))+
  geom_vline(data=filter(reads_mapped,type=="per_flow_cell"),aes(xintercept = 50))+
  geom_text(aes(x=50, label="50M target", y=37), colour="black", angle=90, vjust = 1.2,data=filter(reads_mapped,type=="per_flow_cell"))+
  ylab("Samples") + xlab("Number of reads (M)") +  
  theme_minimal() + theme(axis.text.y = element_text(size=5), legend.position = "none")
p1

#Per sequence quality score
qc_sequence_quality_jfc=read.table('join_flow_cells/multiqc-join-flowcell/multiqc_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt', h=TRUE,sep="\t")
qc_sequence_quality_jfc.T=setNames(data.frame(t(qc_sequence_quality_jfc[,-1])), qc_sequence_quality_jfc[,1])
##add read2
qc_sequence_quality_pfc=read.table('per_flow_cells/multiqc-per-flowcell-cecret/multiqc_data/mqc_fastp-seq-quality-plot_Read_1_After_filtering.txt', header = TRUE, sep = '\t')
qc_sequence_quality_pfc<-qc_sequence_quality_pfc[!grepl("97-97|98-98|99-99|100-100|101-101|102-102|103-103|104-104", 
                                                        qc_sequence_quality_pfc$Sample),]
qc_sequence_quality_pfc.T=setNames(data.frame(t(qc_sequence_quality_pfc[,-1])), qc_sequence_quality_pfc[,1])

#merge of data and transformation of positions values in Row.names column from X1 - X100 to 1-100
seq_qc_merge <- merge(qc_sequence_quality_pfc.T, qc_sequence_quality_jfc.T, by.x=0, by.y=0, all = TRUE, sort = FALSE)
seq_qc_merge$Row.names <- sub('.', '', seq_qc_merge$Row.names)
qc_sequence_quality=seq_qc_merge %>% pivot_longer(-Row.names)

#plot per sequence quality score
t2.rect1 = data.frame (ymin=0, ymax=20, xmin=-Inf, xmax=Inf)
t2.rect2 = data.frame (ymin=20, ymax=28, xmin=-Inf, xmax=Inf)
t1.rect1 = data.frame (ymin=28, ymax=Inf, xmin=-Inf, xmax=Inf)

#Row.names needs to pass as numeric implicit, because it came out of a string and not directly as int
Color=c("#5cb85cff")
p2 = ggplot(data=qc_sequence_quality, aes(x=as.numeric(Row.names), y=value, group=name)) +
  geom_rect(data=t2.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="#d9534fff", color="white",alpha=0.5, inherit.aes = FALSE)+
  geom_rect(data=t2.rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="#f0ad4eff", color="white",alpha=0.5, inherit.aes = FALSE)+
  geom_rect(data=t1.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="#5cb85cff", color="white",alpha=0.5, inherit.aes = FALSE)+
  geom_line(aes(colour = Color),linewidth=0.2) + scale_color_manual(name = "", values=c("Failed" = '#d9534fff', "Pass" = '#5cb85cff', "Warning" = "orange")) +
  #labs(title = "Per Sequence Quality Scores")+
  ylab("Phred score") + xlab("Position (bp)") +
  scale_y_continuous(breaks=seq(0,40,5)) + scale_x_continuous(breaks=seq(0,150,25)) + theme_minimal() + theme(legend.position = "bottom")
p2


##reading tables for coverage
genome_cov_jfc =read.table("join_flow_cells/coverage_jfc.txt",h=T,sep = ",")
genome_cov_pfc =read.table("per_flow_cells/coverage_pfc.txt",h=T,sep = ",")
genome_cov_pfc<-genome_cov_pfc[!grepl("97-97|98-98|99-99|100-100|101-101|102-102|103-103|104-104", 
                                      genome_cov_pfc$sampleId),]
type = c(rep("join_flow_cell", nrow(genome_cov_jfc)), rep('per_flow_cell',  nrow(genome_cov_pfc))) #adding type
genome_cov=rbind(genome_cov_jfc,genome_cov_pfc)
genome_cov=cbind(genome_cov, type)

p4=ggplot(genome_cov,aes(x=CoverageX,y=Value,group=sampleId,color=type)) +
  geom_line(size=1) +
  xlab("Coverage (X)") +
  ylab("Genome Fraction (%)") +
  geom_vline(xintercept = 20,color="grey50",linetype=1,size=1)+
  geom_hline(yintercept = 90,color="orange",linetype=1,size=1)+
  geom_vline(xintercept = 40,color="grey20",linetype=1,size=1)+
  geom_vline(xintercept = 100,color="black",linetype=1,size=1)+
  geom_text(aes(x=20, label="20X coverage", y=30), colour="grey50", angle=90, vjust = 1.2)+
  geom_text(aes(x=40, label="40X coverage", y=30), colour="grey20", angle=90, vjust = 1.2)+
  geom_text(aes(x=100, label="100X coverage", y=30), colour="black", angle=90, vjust = 1.2)+
  geom_text(aes(x=150, label="90% genome fraction", y=90), colour="orange", vjust = 1.2)+
  scale_x_continuous(limits = c(0,200))+
  theme_minimal() + theme(legend.position = "bottom")
p4

#joining plots
mplot=(p1|p2)/(p3|p4) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 15))

pdf("Figure1_DN.pdf",12,11)
mplot
dev.off()
mplot
