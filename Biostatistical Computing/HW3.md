## David Whitney's Homework 4
## Biostatistics 562

### Problem 1
For a differential expression simulation study, say 1000 genes are measured on 100 individuals. Let the first 50 individuals constitute *class1* and the second 50 individuals be *class2*. Suppose that expression levels are independent and normally distributed across genes. For individuals in *class2*, let the last 100 genes be distributed N(0.5,1). For all other genes, assume that the distribution of gene expression is standard normal.


```r
set.seed(1990)

person = function(genes = 1000, de = 100, mean) {
    c(rnorm(genes - de), mean + rnorm(de))
}
```


### 1a.
Run one simulation of the experiment and test for differentially expressed genes (genes with mean differences between the two classes) with a false discovery rate cutoff of 10%. How many genes were called significant? How many of these were false discoveries?


```r
experiment = function(N, mean1 = 0, mean2, genes = 1000, de = 100) {
    class1 = replicate(N, person(mean = mean1))
    class2 = replicate(N, person(mean = mean2))
    return(cbind(class1, class2))
}

N = 50
XP = experiment(N = N, mean2 = 0.5)

class = 1:N
p.values = apply(XP, MARGIN = 1, FUN = function(X) {
    t.test(X[class], X[-class])$p.val
})
adjust.vals = p.adjust(p.values, method = "BH")

fdr = 0.1
n.discoveries = sum(adjust.vals < fdr)
n.false = sum(adjust.vals[1:900] < fdr)
```

We find 37 genes that have a significant mean expression between the two classes and 2 of these were false discoveries, having been observed in the first 900 genes.

### 1b.
Repeat this experiment 1000 times. Of the 100 genes we are searching for, what proportion do we find on average (this is known as the true discovery rate TDR)? Of the genes we call significant what proportion on average are false discoveries?  Relate these ideas to the statistical power of our experiment!


```r
do.one = function(N = 50, meandiff = 0.5, fdr = 0.1) {
    XP = experiment(N = N, mean2 = meandiff)
    class = 1:N
    p.values = apply(XP, MARGIN = 1, FUN = function(X) {
        t.test(X[class], X[-class])$p.val
    })
    adjust.vals = p.adjust(p.values, method = "BH")
    
    n.discoveries = sum(adjust.vals < fdr)
    n.true = sum(adjust.vals[901:1000] < fdr)
    n.false = sum(adjust.vals[1:900] < fdr)
    fdp = n.false/n.discoveries
    tdr = n.true/100
    return(c(FDP = fdp, TDR = tdr))
}

set.seed(1776)
B = 1000

# TrueFalse.mean = apply(replicate(B, do.one()), MARGIN=1, FUN=mean)
```

After running 1000 of the above simulation studies, we found that the proportion of discoveries that were false discoveries was 0.0905 and the proportion of true discoveries was 0.3189. Power is the probability of rejecting when there is a true effect, so the true discovery rate we estimate is an estimate of the power.

### 1c.
Why might this lead to concern about the reproducibility in these microarray experiments?</br>
Test?
