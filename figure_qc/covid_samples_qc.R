knitr::opts_chunk$set(echo = TRUE)
library('tidyverse')
library('openxlsx')
library('ggbeeswarm')
library('patchwork')
library('ggplot2')
library('tidyverse')

## Reads Mapped
samtools_jfc = read.table('join_flow_cells/multiqc-join-flowcell/multiqc_data/multiqc_samtools_stats.txt', header = TRUE)
samtools_pfc = read.table('per_flow_cells/multiqc-per-flowcell-cecret/multiqc_data/multiqc_samtools_stats.txt', header = TRUE)

#selecting data
reads_mapped_jfc <- samtools_jfc %>%
  dplyr::select(Sample, reads_mapped, raw_total_sequences, reads_mapped_percent)
reads_mapped_pfc <- samtools_pfc %>%
  dplyr::select(Sample, reads_mapped, raw_total_sequences, reads_mapped_percent)

#adding type of data
type = c(rep("join_flow_cell", nrow(reads_mapped_jfc)), rep('per_flow_cell',  nrow(samtools_pfc)))
#join both table
reads_mapped=rbind(reads_mapped_jfc, reads_mapped_pfc)
#join type of sample
reads_mapped=cbind(reads_mapped,type)

#violin plot with dots 
p3 =ggplot(reads_mapped,aes(y=reads_mapped_percent,x=type,color=type, label=Sample))+
  geom_violin(alpha = 0.75) +
  geom_point(position = position_jitter(seed = 1, width = 0.35), size = 1.5) +
  xlab("") +
  ylab("Mapped Reads (%)") +
  theme_minimal() +
  theme(legend.position = "none")
p3


#barplot num of reads
#raw_total_sequences = reads_mapped + reads_unmapped
p1<-ggplot(reads_mapped,aes(x=raw_total_sequences/100000,y=fct_reorder(Sample,raw_total_sequences),fill=type)) + 
  geom_col()+
  facet_wrap(~type,nrow =1,scales = "free_y")+
  geom_vline(data=filter(reads_mapped,type=="join_flow_cell"),aes(xintercept = 1500))+
  geom_text(aes(x=1500, label="1500M target", y=1), colour="black", angle=90, vjust = 1.2,data=filter(reads_mapped,type=="join_flow_cell"))+
  geom_vline(data=filter(reads_mapped,type=="per_flow_cell"),aes(xintercept = 1500))+
  geom_text(aes(x=1500, label="1500M target", y=1.3), colour="black", angle=90, vjust = 1.2,data=filter(reads_mapped,type=="per_flow_cell"))+
  ylab("Samples") + xlab("Number of reads (M)") +  
  theme_minimal() + theme(axis.text.y = element_text(size=6), legend.position = "none")
p1

#Per sequence quality score
qc_sequence_quality_jfc=read.table('join_flow_cells/multiqc-join-flowcell/multiqc_data/mqc_fastqc_per_base_sequence_quality_plot_1.txt', h=TRUE,sep="\t")
qc_sequence_quality_jfc.T=setNames(data.frame(t(qc_sequence_quality_jfc[,-1])), qc_sequence_quality_jfc[,1])
qc_sequence_quality_pfc=read.table('per_flow_cells/multiqc-per-flowcell-cecret/multiqc_data/mqc_fastp-seq-quality-plot_Read_1_After_filtering.txt', header = TRUE, sep = '\t')
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

##Formating data to make Coverage plot
genome_cov_jfc =read.table("join_flow_cells/multiqc_general_stats.txt",h=T,sep = "\t")
genome_cov_pfc =read.table("per_flow_cells/multiqc_general_stats.txt",h=T,sep = "\t")
genome_cov_jfc <- genome_cov_jfc %>%
  dplyr::select(Sample,QualiMap_mqc.generalstats.qualimap.1_x_pc, QualiMap_mqc.generalstats.qualimap.5_x_pc ,
                QualiMap_mqc.generalstats.qualimap.10_x_pc, QualiMap_mqc.generalstats.qualimap.30_x_pc, QualiMap_mqc.generalstats.qualimap.50_x_pc)

genome_cov_pfc <- genome_cov_pfc %>%
  dplyr::select(Sample,QualiMap_mqc.generalstats.qualimap.1_x_pc, QualiMap_mqc.generalstats.qualimap.5_x_pc ,
                  QualiMap_mqc.generalstats.qualimap.10_x_pc, QualiMap_mqc.generalstats.qualimap.30_x_pc, QualiMap_mqc.generalstats.qualimap.50_x_pc)

#adding columns names and prefix to detect type of sample
colnames(genome_cov_jfc) <- c('Coverage', 1, 5, 10, 30, 50)
genome_cov_jfc$Coverage = paste0('jfc_', genome_cov_jfc$Coverage)
colnames(genome_cov_pfc) <- c('Coverage', 1, 5, 10, 30, 50)
genome_cov_pfc$Coverage = paste0('pfc_', genome_cov_pfc$Coverage)

#joining both datasets and transforming data
genome_cov=rbind(genome_cov_jfc,genome_cov_pfc)
genome_cov.t <- as.data.frame(t(genome_cov))
colnames(genome_cov.t) <- genome_cov.t[1,]
genome_cov.t <- genome_cov.t[-1,]
genome_cov.t <- tibble::rownames_to_column(genome_cov.t, "Coverage") #transform index into a column
genome_cov_filter=genome_cov.t %>% pivot_longer(-Coverage)
genome_cov_filter$type[str_detect(genome_cov_filter$name,"pfc")]<-"per_flow_cell"
genome_cov_filter$type[str_detect(genome_cov_filter$name,"jfc")]<-"join_flow_cell"
genome_cov_filter$value <- as.numeric(genome_cov_filter$value)
genome_cov_filter$Coverage <- as.numeric(genome_cov_filter$Coverage)

genome_cov_filter %>%filter(round(value)==90) %>%group_by(type) %>% summarize(mean=mean(Coverage),sd=sd(Coverage),median=median(Coverage))

#genome coverage plot
p4=ggplot(genome_cov_filter,aes(x=Coverage,y=value,group=name,color=type)) +
  #geom_line(linetype=3) + #dotted line
  geom_line(size=1) +
  xlab("Coverage (X)") +
  ylab("Genome Fraction (%)") +
  #ggtitle("Cumulative Genome Coverage")+
  #xlim(0,100)+
  geom_vline(xintercept = 20,color="grey",linetype=1,size=1)+
  geom_hline(yintercept = 90,color="orange",linetype=1,size=1)+
  geom_vline(xintercept = 40,color="black",linetype=1,size=1)+
  geom_text(aes(x=20, label="20X coverage", y=30), colour="grey", angle=90, vjust = 1.2)+
  geom_text(aes(x=40, label="40X coverage", y=30), colour="black", angle=90, vjust = 1.2)+
  geom_text(aes(x=130, label="90% genome fraction", y=90), colour="orange", vjust = 1.2)+
  theme_minimal() + theme(legend.position = "bottom")
p4

mplot=(p1|p2)/(p3|p4) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 15))

pdf("Figure1_DN.pdf",12,11)
mplot
dev.off()
mplot
