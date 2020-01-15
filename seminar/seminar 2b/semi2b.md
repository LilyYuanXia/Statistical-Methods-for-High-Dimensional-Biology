Seminar 2b: Graphing using ggplot2
================
Lily Xia
15 January, 2020

``` r
mpg %>% 
ggplot() +
  geom_point(mapping = aes(x = displ, y = hwy, size = class, color = drv))
```

![](semi2b_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
