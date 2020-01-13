---
title: 'Seminar 2a: Introduction to R Markdown'
author: "Yuan Xia"
date: "13 January, 2020"
output: 
  html_document: 
    keep_md: yes
---


```
## ── Attaching packages ──────────── tidyverse 1.2.1 ──
```

```
## ✔ ggplot2 3.2.1     ✔ purrr   0.3.2
## ✔ tibble  2.1.3     ✔ dplyr   0.8.3
## ✔ tidyr   1.0.0     ✔ stringr 1.4.0
## ✔ readr   1.3.1     ✔ forcats 0.4.0
```

```
## ── Conflicts ─────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

# R code
R Code can be directly embedded into an .Rmd file and will be executed during rendering. The result of the code would be printed into resulting report.

```r
summary(cars) %>% knitr::kable()
```

         speed           dist      
---  -------------  ---------------
     Min.   : 4.0   Min.   :  2.00 
     1st Qu.:12.0   1st Qu.: 26.00 
     Median :15.0   Median : 36.00 
     Mean   :15.4   Mean   : 42.98 
     3rd Qu.:19.0   3rd Qu.: 56.00 
     Max.   :25.0   Max.   :120.00 


# Github Link
[Link return to my GitHub](https://github.com/LilyYuanXia)

