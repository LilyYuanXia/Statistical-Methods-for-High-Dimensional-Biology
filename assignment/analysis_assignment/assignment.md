Analysis Assignment
================
Lily Xia

## Question 1: Data Inspection and Basic Manipulation

### Q1.1 Importing the data and getting familiar with it (2 POINT)

`gene_data` is the transcriptomic data, and `meta_data` stores the
samples metadata.

``` r
gene_data <- readRDS("gse60019_expression_matrix.RDS")
meta_data <- readRDS("gse60019_experiment_design.RDS")
```

There are 14479 unique genes from 18 samples in the given expression
matrix. In the sample metadata, there are 4 variables for the 18
samples, including `organism_part`, `cell_type`, `time_point` and
`batch`. These 4 variables are all in factor class. The levels for each
variable are listed below:

``` r
dat <- meta_data[,-1]
data.frame(variable = names(dat),
           class = sapply(dat, class),
           levels = sapply(dat, function(x) paste0(levels(x),  collapse = ", ")),
           row.names = NULL) %>% 
kable()
```

| variable       | class  | levels                                                          |
| :------------- | :----- | :-------------------------------------------------------------- |
| organism\_part | factor | epithelium\_of\_utricle, sensory\_epithelium\_of\_spiral\_organ |
| cell\_type     | factor | surrounding\_cell, sensory\_hair\_cell                          |
| time\_point    | factor | E16, P0, P4, P7                                                 |
| batch          | factor | HWI-EAS00184, HWI-EAS00214, HWI-ST363                           |

### Q1.2 Data manipulation (2 POINTS)

Create an new column in the samples metadata tibble to store the numeric
`age` values which are converted from the `time_point` column. From the
[Mouse Timeline Detailed
website](https://embryology.med.unsw.edu.au/embryology/index.php/Mouse_Timeline_Detailed),
we can see that `E16` is Day 16 in week. Since we assume that the mouse
gestation length is 18, `P0` is Day 18, `P4` is Day 22, `P7` is Day 25.
Here I only show 4 representative rows from the whole datatable due to
the room limit.

``` r
# initial assumption for gestation length
a <- 18
meta_data$age <- ifelse(meta_data$time_point %>% substring(1,1) == "E",
                        as.numeric(meta_data$time_point %>% str_extract("\\d+")),
                        as.numeric(meta_data$time_point %>% str_extract("\\d+")) + a)
meta_data[c(1,5,7,18),] %>% kable()
```

| sample     | organism\_part                         | cell\_type          | time\_point | batch        | age |
| :--------- | :------------------------------------- | :------------------ | :---------- | :----------- | --: |
| GSM1463874 | epithelium\_of\_utricle                | sensory\_hair\_cell | E16         | HWI-EAS00214 |  16 |
| GSM1463876 | sensory\_epithelium\_of\_spiral\_organ | sensory\_hair\_cell | P0          | HWI-ST363    |  18 |
| GSM1463884 | epithelium\_of\_utricle                | sensory\_hair\_cell | P4          | HWI-EAS00214 |  22 |
| GSM1463885 | sensory\_epithelium\_of\_spiral\_organ | surrounding\_cell   | P7          | HWI-EAS00184 |  25 |

### Q1.3 Single gene graphing (3 POINTS)

``` r
vegfa <- gene_data[gene_data$gene == "Vegfa",]
melt_vegfa <- vegfa %>% melt(vars.id = "gene",
                             var = "sample")
full_vegfa <- right_join(melt_vegfa, meta_data, by = "sample")
full_vegfa %>% 
  ggplot(aes(x = age, y = value, color = cell_type)) +
  geom_point() +
  geom_smooth(method='lm', alpha = 0.1) +
  theme_bw() +
  labs(y = "Expression values are given in CPM")
```

![](assignment_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

From the scatterplot of age and expression values in CPM, two linear
regression lines based on different cell types are parallel. According
to this characteristics, there is no sign of interaction between cell
type and age for Vegfa.

## Question 2: Assessing overall data quality

### Q2.1 Overall distributions (4 POINTS)

### Q2.2 How do the samples correlate with one another? (4 POINTS)
