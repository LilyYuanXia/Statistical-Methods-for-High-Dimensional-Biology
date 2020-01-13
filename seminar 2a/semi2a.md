---
title: 'Seminar 2a: Introduction to R Markdown'
author: "Yuan Xia"
date: "13 January, 2020"
output: 
  html_document: 
    keep_md: yes
---



### Outlines of Seminar 2a
  * R code chunk
  * Basic features
    * Link
    * Images

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

# Basic features

### Link
[Link return to my GitHub](https://github.com/LilyYuanXia)

### Images
![](https://i1.wp.com/www.chinatownfoundation.org/wp-content/uploads/2015/05/ubc-logo.png?zoom=2&ssl=1)
