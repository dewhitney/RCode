## David Whitney's Homework 4
## Biostatistics 562

### Problem 1
For a differential expression simulation study, say 1000 genes are measured on 100 individuals. Let the first 50 individuals constitute *class1* and the second 50 individuals be *class2*. Suppose that expression levels are independent and normally distributed across genes. For individuals in *class2*, let the last 100 genes be distributed N(0.5,1). For all other genes, assume that the distribution of gene expression is standard normal.




### 1a.
Run one simulation of the experiment and test for differentially expressed genes (genes with mean differences between the two classes) with a false discovery rate cutoff of 10%. How many genes were called significant? How many of these were false discoveries?




We find 37 genes that have a significant mean expression between the two classes and 2 of these were false discoveries, having been observed in the first 900 genes.

### 1b.
Repeat this experiment 1000 times. Of the 100 genes we are searching for, what proportion do we find on average (this is known as the true discovery rate TDR)? Of the genes we call significant what proportion on average are false discoveries?  Relate these ideas to the statistical power of our experiment!




After running 1000 of the above simulation studies, we found that the proportion of discoveries that were false discoveries was 0.0905 and the proportion of true discoveries was 0.3189. Power is the probability of rejecting when there is a true effect, so the true discovery rate we estimate is an estimate of the power.

### 1c.
Why might this lead to concern about the reproducibility in these microarray experiments?</br></br>

If the power of microarray studies is only 0.3189 and the FDR is being controlled at near .0905, then we are unlikely to detect true effects and we will further expect that of the effects that we are detecting, nearly 10% of them will be false discoveries.

### Problem 2
We now explore the scenarios under which using FDR as opposed to FWE gives significant additional power.

### 2a.
We consider the study design from _Problem 1_ but assume now that the mean expression for the last 100 genes of _class2_ follow a N(2,1) distribution. 




When we control the FDR at 10%, we observe a false discovery proportion of 0.0895 and 100% TDR. When we control FWE at .1, we observe a false discovery proportion of 0.0007 and 100% TDR. In this case, the Bonferroni correction outperforms Benjamini-Hochberg by having equal power and a lower Type I error rate.

### 2b.
We consider the study design from _Problem 1_ but assume now that the mean expression for the last 100 genes of _class2_ follow a N(0.2,1) distribution. 




When controlling FDR, the FDP was 0.08793 and TDR 0.00251. For the Bonferroni correction of FWE, the FDP was 0.07617 and the TDR was 0.00181. Both manners of adjusting the p-values had similar powers in this case.
<\br>
<\br>





When the differentially expressed genes followed N(0.01,1), ..., N(1.00,1) distributions, controlling the FDR at .10 yielded an estimated false discovery proportion of 0.08933 and true discovery rate of 0.4309. Controlling the FWE instead yielded FDP of 0.003688 and TDR of 0.2146. In this scenario, the TDR obtained by controlling FDR is about twice the power obtained by controlling the FWE.

### Problem 3



The proportion discoveries was 0.026 when controlling FDR at .10 and 0.025 when controlling the FWE at .10.




The proportion of discoveries was 0.021 for both methods.
