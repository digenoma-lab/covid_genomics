analysis_covid
================
2023-08-03

## Reading files

Reading results from iVar and freya from join_flow_cell pipeline and
including demographic information of samples.

    ## Using X as id variables

### Mutation Heatmap

Using ComplexHeatmap library to plot if samples has 1 or any mutations
detected by iVar, using a binary matrix. Variant a location is also
included into the heatmap.

![](analysis_files/figure-gfm/heatmapivar-1.png)<!-- -->

### Number of variants vs heterogeneity

Combining mutations detected by iVar and variations in to positions
detected by Freya. Including variant information and location of
samples.

![](analysis_files/figure-gfm/mutations-1.png)<!-- -->
