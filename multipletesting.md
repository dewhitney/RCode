Multiple testing
========================================================


```r
alpha = 0.05
n = 1:100
fdr = 1 - (1 - alpha)^n
plot(n, fdr)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 

