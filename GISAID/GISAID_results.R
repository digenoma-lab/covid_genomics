library('tidyverse')
library('patchwork')
library('ggplot2')
library('reshape2')
library('lubridate')

#reading data
all_results <-read.table('/Users/bbernal/Documents/covid_work/all_results.csv', h=TRUE,sep=",")
metadata = read.table('/Users/bbernal/Documents/covid_work/GISAID/chile_data/GISAID_all_chile.tsv', 
                      h=TRUE,sep="\t", check.names = FALSE)

#formating dates
metadata$Date <-ym(strtrim(metadata$`Collection date`, 7))
variants <- metadata[,c("Lineage","Date")]

#extracting dates of sampling data
samples_dates <- ym(strtrim(all_results$analisis_data, 7))
samples_dates <- na.omit(samples_dates)
samples_dates <- unique(samples_dates)


#generating table for bar plot
hist_variants <- variants %>%
  group_by(Date) %>% 
  count(Lineage)

#formating lineage info to variant names, some links for manual search of VOI and VOC
#https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-classifications.html
#https://gisaid.org/hcov19-variants/

hist_variants <- hist_variants %>%
  mutate(
    name = case_when(
      str_detect(Lineage, "P.2") ~ "Zeta",
      str_detect(Lineage, "B.1.617.1") ~ "Kappa",
      str_detect(Lineage, "B.1.526") ~ "Iota",
      str_detect(Lineage, "B.1.525") ~ "Eta",
      str_detect(Lineage, "B.1.427|B.1.429") ~ "Epsilon",
      str_detect(Lineage, "B.1.1.529|BA.2.86|CH.1.1|BA.2.74|^BA.\\d|XBB.1.\\d") ~ "Omicron",
      str_detect(Lineage, "B.1.617|B.1.617.2|^AY.\\d") ~ "Delta",
      str_detect(Lineage, "B.1.1.7|^Q.\\d") ~ "Alpha",
      str_detect(Lineage, "P.1|^P.1.\\d") ~ "Gamma",
      str_detect(Lineage, "B.1.351|B.1.351.2|B.1.351.3") ~ "Beta",
      str_detect(Lineage, "C.37|C.37.1") ~ "Lambda",
      str_detect(Lineage, "B.1.621|B.1.621.1") ~ "Mu",
      .default = "Other"
    )
  )

#calculating percentage of variants in each month
variant_percentage <- hist_variants %>% 
  group_by(Date) %>%
  count(name) %>%
  mutate(percent = n/sum(n))

# Number of samples 
p1 <- ggplot(hist_variants) + 
  geom_bar(aes(x=Date, y=n, fill=name), stat="identity") +
  annotate("segment", x=samples_dates, xend=samples_dates,
    y=300, yend=2, color="blue4",
    arrow=arrow(length=unit(0.2,"cm"))) +
  scale_x_date(date_labels="%m-%Y",date_breaks  ="1 month") +
  theme_minimal()+
  xlab("Date") + ylab("# samples") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 9))
p1

# variants percentage per month
p2 <- ggplot(variant_percentage) + 
  geom_bar(aes(x=Date, y=percent, fill=name), stat="identity") +
  #  annotate("segment", x=samples_dates, xend=samples_dates,
  #  y=0.25, yend=0.05, color="blue4",
  #  arrow=arrow(length=unit(0.2,"cm"))) +
  scale_x_date(date_labels="%m-%Y",date_breaks  ="1 month") +
  theme_minimal()+
  xlab("Date") + ylab("variants(%)") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, size = 9))

p2
