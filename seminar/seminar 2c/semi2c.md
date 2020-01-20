seminar 2c
================
Yuan Xia
20/01/2020

``` r
library(tidyverse)
```

## CLT With Sample Size 5

``` r
set.seed(1)

sampleSize <- 5
numSamples <- 1000

degreeFreedom <- 1

randomChiSqValues <- rchisq(n = numSamples * sampleSize, df = degreeFreedom)

samples <- matrix(randomChiSqValues, nrow = numSamples, ncol = sampleSize)
sampleMeans <- rowMeans(samples)


tibble(x = sampleMeans) %>% 
  ggplot() + 
  geom_line(aes(x = x), stat = "density", color = "blue") +
  geom_point(aes(x = x, y = 0), color = "blue", shape = 1, size = 3)
```

![](semi2c_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
