Covid Figures
================
Alex Di Genova
2024-03-02

# Figure 1

Demographic of sequenced COVID samples. \## Chilean map

Map of Chile hihglithing the Rancagua y Santiago regions
![](MFig_Covid_files/figure-gfm/chimap-1.png)<!-- -->

## Demographic of samples

``` r
df=read.table("../demographics/all_samples.txt", h=T)
dates_formatted <- as.Date(df$Date, format="%m/%d/%Y")
years <- format(dates_formatted, "%Y")
df$year=years
```

we make the plot

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ lubridate 1.9.3     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.2     ✔ tidyr     1.3.1
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggbeeswarm)
count_data <- df %>%
  group_by(year) %>%
  summarise(count = n())

p2=ggplot(df, aes(x = year, y = Age, )) +
     #geom_boxplot(aes(group = year), color="brown", alpha = 0.5) + # Add boxplot
  geom_quasirandom(aes(color = Region, shape = Sex),width = 0.25, height = 0, size=3) +
  # geom_beeswarm(aes(color = Region, shape = Sex),width = 0.25, height = 0, size=3) +
   #geom_jitter(, width = 0.2, height = 0) +
  scale_color_manual(values = c("Metropolitana" = "blue", "O'Higgins" = "red")) +
  labs(title = "Samples",
       x = "Year",
       y = "Age") +
    theme(text = element_text(size = 12)) +
  theme_minimal()
```

    ## Warning in geom_quasirandom(aes(color = Region, shape = Sex), width = 0.25, :
    ## Ignoring unknown parameters: `height`

``` r
p2
```

![](MFig_Covid_files/figure-gfm/demogp-1.png)<!-- -->

## GISAD

## Plotting number of samples and percentage of variants in Chile.

## Number of samples

Here we got the number of samples in Chile from GISAID platform, we
label this with their corresponding name of VOI or VOC. With blue arrows
we point out where our samples are located on the timeline. Around years
2021 and 2022 there is huge increase in the sequencing of COVID samples
in Chile.

![](MFig_Covid_files/figure-gfm/num_samples-1.png)<!-- -->

We compose all plots to generate figure 1

``` r
library(patchwork)

patchwork=p1|(p2/p3)  
patchwork + plot_annotation(tag_levels = 'A')
```

![](MFig_Covid_files/figure-gfm/cmposer-1.png)<!-- -->

``` r
pdf(file="Fig1.pdf",height=7,width=9)
patchwork + plot_annotation(tag_levels = 'A')
dev.off()
```

    ## quartz_off_screen 
    ##                 2

# Figure 2

Philogenetics, variants detected and coverage.

## Phylogenetic tree

We plot the tree including variants type, iVAR snp distribution and
number of SNPs

``` r
library(ggtree)

# Customize further as needed
p4 <- ggtree(x_new, color="grey40", ladderize = TRUE) +
  #xlim_tree(0.00325) +
  geom_tiplab(aes(color=variant, label = sampleId), size=2.) +
  scale_color_manual('', values = variant) +
  scale_fill_manual('', values = variant) +
  guides(color = guide_legend(override.aes = list(size = 5, label = "\u25AA")))+
  theme_tree2(legend.position = "bottom")

# add SNP data by position in the tree
p5= p4 + geom_facet(panel = "SNPs coordinates", data = ivdt_longer, geom = geom_point, 
               mapping=aes(x = pos,group=sampleId, color=variant,fill=variant), shape = '|') +
             geom_facet(panel = "# SNPs", data = ivar_c, geom = geom_col, 
                aes(x = value, color = variant, fill=variant), orientation = 'y', width = .6)
p5
```

![](MFig_Covid_files/figure-gfm/plott-1.png)<!-- -->

``` r
#facet_widths(p8, widths = c(2, 2, 1))
pdf(file="Fig2.pdf",height=7,width=9)
facet_widths(p5, widths = c(2, 2, 1))
```

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <e2>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <96>

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## conversion failure on '▪' in 'mbcsToSbcs': dot substituted for <aa>

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

# Figure 3

Genome coverage

## Samtools depth

Checking depth of the samples and separate them by variant and location.

``` r
suppressMessages(library('reshape2'))
samtools_results <-read.table('samtools_depth.csv', h=TRUE,sep=",")
samtools_results <- melt(samtools_results, id="X")
colnames(samtools_results) <- c('pos', 'sampleId', 'depth')

samtools_results <- left_join(samtools_results,all_results[c('sampleId', 'location', 'variant')], by="sampleId")
samtools_results<-samtools_results[!grepl("71", samtools_results$sample),]

samtools_results <- samtools_results %>% group_by(variant)

samtools_results <- samtools_results %>% group_by(variant) %>%
  mutate(med = mean(depth, na.rm = TRUE))

samtools_depth_range <- samtools_results %>%
  group_by(sampleId, grp = cut(pos, breaks=pretty(pos, n = 60), dig.lab = 5),  variant, med) %>% 
  summarise(count = mean(depth, na.rm = TRUE))
```

    ## `summarise()` has grouped output by 'sampleId', 'grp', 'variant'. You can
    ## override using the `.groups` argument.

### Samtools depth plot

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: Removed 7 rows containing missing values (`geom_line()`).

![](MFig_Covid_files/figure-gfm/depth-1.png)<!-- -->

    ## Warning: Removed 7 rows containing missing values (`geom_line()`).

    ## quartz_off_screen 
    ##                 2

# Figure 4

Heterogeneity analysis

``` r
library(tidyverse)
dh=read.table("summary_het2.txt",h=F)
dh$V1=paste0('sample_',(str_extract(dh$V1,"(?<=-)\\d+")))
dh=left_join(dh,all_results,by=c("V1"="sampleId"))
dh=dh %>% mutate(age2=if_else(age < 30,"young",if_else(age >=30 & age<50,"adult","older"))) %>%  mutate(variant=if_else(variant=="BA.2* [Omicron (BA.2.X)]","Omicron",variant)) %>%remove_missing()
```

    ## Warning: Removed 2 rows containing missing values.

``` r
#count the number of position with ITH
dh$het=dh$V4+dh$V6

#run the lm
a=lm(dh$het~dh$age2+dh$mut_ivar+dh$variant)
summary(a)
```

    ## 
    ## Call:
    ## lm(formula = dh$het ~ dh$age2 + dh$mut_ivar + dh$variant)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -195.647  -67.468   -5.096   55.407  305.431 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        590.365    112.643   5.241 1.17e-06 ***
    ## dh$age2older         1.661     25.769   0.064  0.94876    
    ## dh$age2young       -51.137     25.694  -1.990  0.04981 *  
    ## dh$mut_ivar         -4.945      1.709  -2.894  0.00484 ** 
    ## dh$variantEpsilon  -20.318    145.178  -0.140  0.88903    
    ## dh$variantGamma    -61.479    109.149  -0.563  0.57476    
    ## dh$variantLambda   -61.929    103.711  -0.597  0.55203    
    ## dh$variantOmicron  112.733    113.185   0.996  0.32211    
    ## dh$variantOther   -140.196    105.567  -1.328  0.18777    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 100.3 on 84 degrees of freedom
    ## Multiple R-squared:  0.1605, Adjusted R-squared:  0.08054 
    ## F-statistic: 2.007 on 8 and 84 DF,  p-value: 0.05524

``` r
#ggplot(dh, aes(y=het,x=age2,color=variant,group=age2))  + geom_boxplot() + geom_quasirandom()
```
