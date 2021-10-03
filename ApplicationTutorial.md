Introduction
============

As described in the [`README.md`](https://github.com/augustofasano/Probit-PFMVB/blob/master/README.md) file, this tutorial contains general guidelines and code to **perform the analyses for the medical applications in Section 4** of the paper.
In particular, the tutorial is divided in two main sections: the first part contains detailed information on how to replicate the **Alzheimer's** application, while the second one shows how to **compute the deviances** reported in Table 2 for the **lesion** dataset. The *parkinson* and *voice* results can be obtained by running the code for the *lesion* case after retrieving the desired datasets from the open source [UCI repository](https://archive.ics.uci.edu/ml/datasets.php).
For implementation purposes, execute the code below considering the same order in which is presented.

The Alzheimer's application
=====================

The code below shows how to **load the data** of the Alzheimer's study presented in Section 4 and explains in detail the `R` code to **implement the different methods for posterior inference** discussed in Section 2. Guidelines to **produce Figures 3 to 6 in the paper** are also provided.

As discussed in Section 4, the focus is to **model presence or absence of Alzheimer’s disease in its early stages as a function of demographic data, genotype and assay results**. The original dataset is available in the `R` library `AppliedPredictiveModeling` and arises from a study of the Washington University ([Craig-Schapiro et al., 2011](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0018850)).

Here, we rely on an interpretable **probit regression** with pairwise interactions, thus obtaining *p* = 9036 predictors collected for 333 individuals. Following [Gelman et al. (2008)](https://projecteuclid.org/euclid.aoas/1231424214) and [Chopin and Ridgway (2017)](https://projecteuclid.org/euclid.ss/1491465628) the original measurements which are used in the training set have been standardized to have mean 0 and standard deviation 0.5.
The same scaling transformation has then been applied to the test set, so that both training and test sets are transformed in the same way.
The vector of the observations falling in the training set, `trainingSet`, is also provided together with the data.
The final standardized dataset is available in the [`data`](https://github.com/augustofasano/Probit-PFMVB/tree/master/data) folder and is called `Alzheimer_Interactions.RData`.

Preliminary operations
======================

Once the files `Alzheimer_Interactions.RData` and `functionsVariational.R` have been downloaded, set the working directory to the folder where they are located. Then, clean the workspace and load `Alzheimer_Interactions.RData` along with the source file `functionsVariational.R` and other useful packages.

``` r
rm(list=ls())
source("functionsVariational.R")
library("TruncatedNormal")
library("truncnorm")
library("transport")

load("Alzheimer_Interactions.RData")
```

The vector `trainingSet` contains the *n* = 300 observations that are going to be used as training set to perform posterior inference, while the remaining 33 units are held out for computing the out-of-sample predictive probabilities.

``` r
seed = 1
set.seed(seed)
n = 300 # training set size

# split training and test sets
X_Test = X[-trainingSet,]
yTest = y[-trainingSet]
X = X[trainingSet,]
y = y[trainingSet]
```

Finally, we conclude the preliminary operations by specifying the **model dimension** *p*, the **prior variance** *ν*<sub>*p*</sub><sup>2</sup> and the **number of i.i.d samples** to generate from the exact and approximate posterior densities. We also precompute **V** **X**<sup>⊺</sup> and (**I**<sub>*n*</sub> + *ν*<sub>*p*</sub><sup>2</sup>**X** **X**<sup>⊺</sup>)<sup> − 1</sup> to be used in the computation of the predictive probabilities.

``` r
p = dim(X)[2] # number of covariates
nu2 = 25 # prior variance
nSample = 2e4 # fix number of samples

# precompute some useful quantities to be used for the predictive probabilities
VXt = t(nu2*X)%*%solve(diag(n)+(nu2*X)%*%t(X))
invIXXt = solve(diag(1,nrow=n,ncol=n)+nu2*(X%*%t(X)))
```

Partially factorized mean-field variational Bayes
=================================================

We start our analysis by implementing our proposed **partially factorized mean-field variational** approximation. Consistent with this goal, we consider the function `getParamsPFM` which provides the optimal *q*<sup>\*</sup><sub>PFM</sub>(**β**, **z**) = *p*(**β** ∣ **z**)*q*<sup>\*</sup><sub>PFM</sub>(**z**) via the **CAVI** in **Algorithm 2**.

The function also outputs the  **approximate posterior means and variances** of the coefficients in **β**. Moreover, we also compute the **approximate predictive probabilities** for the 33 observations in the test set via Monte-Carlo integration as discussed in Section 2.2 of the article.

For comparison, we keep track of the **running time** of the algorithm to converge to the optimal solution and to compute the first two marginal moments, `timeSUN_PFM_algo`, along with the additional time to compute the predictive probabilities, `timeSUN_PFM_inference`.

``` r
tolerance = 1e-3 # tolerance to establish ELBO convergence
# get optimal parameters and moments and predictive probabilities
startTime = Sys.time()
paramsPFM = getParamsPFM(X=X,y=y,nu2=nu2,moments=TRUE,tolerance=tolerance,maxIter=1e4)
timeSUN_PFM_algo = difftime(Sys.time(), startTime, units=("secs"))[[1]]

# predictive probabilities
# obtain a sample from q^*(z) to be used for the predictive probabilities of PFM
startTime = Sys.time()
nSampleZ = 5e3
muTN = paramsPFM$mu
muTN[y==0] = -muTN[y==0]
sampleTruncNorm = matrix(rtruncnorm(n*nSampleZ, a = 0, b = Inf, mean = muTN, sd = sqrt(paramsPFM$sigma2)), nrow = n, ncol = nSampleZ, byrow = F )
sampleTruncNorm[y==0,] = -sampleTruncNorm[y==0,] # need to adjust the sign of the variables for which y_i is 0

predPFM = double(length = length(yTest))
for(i in 1:length(yTest)){
  xNew = matrix(X_Test[i,],ncol = 1)
  Xx = X%*%xNew
  sd = as.double(sqrt(1+nu2*(sum(xNew^2)-nu2*t(Xx)%*%invIXXt%*%Xx)))

  predPFM[i] = mean(pnorm((t(xNew)%*%VXt%*%sampleTruncNorm)/sd))
}
timeSUN_PFM_inference = difftime(Sys.time(), startTime, units=("secs"))[[1]]
```

Mean-field variational Bayes
============================

Here, we implement the code to **obtain the optimal parameters of** *q*<sup>\*</sup><sub>MF</sub>(**β**) as in **Algorithm 1**, using `getParamsMF`.

Also this function outputs the **approximate posterior means and variances** of the coefficients in **β**. Moreover, we also compute the **approximate predictive probabilities** for the 33 units using the exact formula in Section 2.1 of the article.

As done before, the **running times** of the algorithm and of the inference part are monitored.

``` r
startTime = Sys.time()
paramsMF = getParamsMF(X,y,nu2,tolerance,maxIter = 1e4)
timeMF_algo = difftime(Sys.time(), startTime, units=("secs"))[[1]]

# predictive probabilities
startTime = Sys.time()
predMF = double(length = length(yTest))
for(i in 1:length(yTest)){
  xNew = matrix(X_Test[i,],ncol = 1)
  Xx = X%*%xNew
  sd = as.double(sqrt(1+nu2*(sum(xNew^2)-nu2*t(Xx)%*%invIXXt%*%Xx)))

  predMF[i] = as.double(pnorm(t(xNew)%*%paramsMF$meanBeta/sd))
}

timeMF_inference = difftime(Sys.time(), startTime, units=("secs"))[[1]]
```

Sampling from the exact unified skew-normal posterior
========================

In order to have a benchmark to assess the quality of the approximations provided by the above methods, we consider also **i.i.d. samples from the exact unified skew-normal posterior** drawn via the algorithm by [Durante (2019)](https://doi.org/10.1093/biomet/asz034). This can be done using the function `rSUNpost`.

Also in this case, we monitor the **running time** of the sampler and of the associated Monte-Carlo inference strategies.

``` r
# obtain a sample from the posterior distribution
startTime = Sys.time()
betaSUN = rSUNpost(X=X,y=y,nu2=nu2,nSample=nSample)
timeSUN_algo = difftime(Sys.time(), startTime, units=("mins"))[[1]]

# predictive probabilities
startTime = Sys.time()
predSUN = double(length = length(yTest))
for(i in 1:length(yTest)){
  xNew = matrix(X_Test[i,],ncol = 1)

  predSUN[i] = mean(pnorm(t(xNew)%*%betaSUN))
}
timeSUN_inference = difftime(Sys.time(), startTime, units=("secs"))[[1]]
```

Running times and number of iterations
======================================

We can now **check the running times** (in minutes) of the three methods and the **number of iterations** required by the **CAVI** for **PFM-VB** and **MF-VB** to reach the optimal solution.

``` r
timeSUN_algo + timeSUN_inference/60
timeSUN_PFM_algo/60 + timeSUN_PFM_inference/60
timeMF_algo/60 + timeMF_inference/60

paramsPFM$nIter
paramsMF$nIter
```
The **CAVI** for **PFM-VB** converges in 6 iterations, whereas that for  **MF-VB** requires 212 iterations.

**Remark**: Depending on the machine, the running times can slightly differ from those reported in the paper, but they should not significantly deviate from them.

Summaries of the results
========================

We now compute different **summary quantities** which are required to produce the Figures in Section 4 of the paper.

Predictive probabilities
------------------------

We start by saving the **predictive probabilities** previously computed.

``` r
predProb = cbind(predSUN,predMF,predPFM,yTest)
colnames(predProb) = c("SUN", "MF", "PFM", "y")

save(predProb, file="predProb.RData")
```

Wasserstein distances
---------------------

Let us now compute the **Wasserstein distances** between the *p* = 9036 approximate marginal densities provided by the two VB methods and the exact posterior marginals. These distances are computed via Monte-Carlo methods based on 20000 i.i.d. samples from the approximate and exact marginals.


Sampling from the optimal **PFM-VB** approximating density is performed through the function `sampleSUN_PFM`, while samples from the optimal **MF-VB** approximating Gaussian posterior distribution are obtained in a standard way, by exploiting the Choleski decomposition of the optimal covariance matrix **V**. **Note** that `sampleSUN_PFM` produces samples from the joint approximate density and, hence, requires more computational time relative to sampling from the marginals (see the discussion in the article). Although the latter strategy would have been sufficient to compute the **Wasserstein distances**, we still consider `sampleSUN_PFM` to show an implementation of **Algorithm 3** in the article.

The Wasserstein distances are computed with the function `wasserstein1d` from the `R` package `transport`.

``` r
# sample from approximate PFM posterior
set.seed(seed+2)
betaSUN_PFM = sampleSUN_PFM(paramsPFM=paramsPFM,X=X,y=y,nu2=nu2,nSample=nSample)

# sample from approximate MF posterior
set.seed(seed+3)
V = diag(rep(nu2,p),p,p) - t(nu2*X)%*%solve(diag(n)+(nu2*X)%*%t(X))%*%(nu2*X)
L = t(chol(V))
betaMF = L%*%matrix(rnorm(nSample*p),p,nSample)
betaMF = apply(betaMF,2, function(x) paramsMF$meanBeta+x)

wassMF = double(length = p)
for(i in 1:p) {
  wassMF[i] = wasserstein1d(a = betaSUN[i,], b = betaMF[i,], p = 1)
}

wassPFM = double(length = p)
for(i in 1:p) {
  wassPFM[i] = wasserstein1d(a = betaSUN[i,], b = betaSUN_PFM[i,], p = 1)
}
```

Now, we find the indexes of the parameters for which we have the **best and worst approximations**, for both **MF-VB** and **PFM-VB**. Then, we **save the samples obtained for these four parameters** when sampling from the exact posterior and
from the two approximate methods. These quantities are required to produce Figure 4.

``` r
# find best/worst Wasserstein distance for MF and PFM
worstWassMF = which(wassMF==max(wassMF))
bestWassMF = which(wassMF==min(wassMF))
worstWassPFM = which(wassPFM==max(wassPFM))
bestWassPFM = which(wassPFM==min(wassPFM))

# save an array containing the samples for these 4 scenarios
listBestWorst = list()
indexList = c(bestWassMF,worstWassMF,bestWassPFM,worstWassPFM)
count = 1
for(i in indexList){
  listBestWorst[[count]]=cbind(SUN = betaSUN[i,],MF = betaMF[i,],PFM = betaSUN_PFM[i,])
  count = count+1
}
names(listBestWorst) = c("bestMF","worstMF","bestPFM","worstPFM")

save(listBestWorst, file = "BestWorstScenarios.RData")
```

Since the **Wasserstein distances** are computed via Monte-Carlo, there will be some some variability. To quantify such a **Monte-Carlo error**, we also compute the Wasserstein distances  between two different samples of 20000 draws from the exact posterior marginals.

``` r
rm(list=c("betaSUN_PFM","betaMF"))
set.seed(seed+100)
betaSUNcomparison = rSUNpost(X,y,nu2,nSample = nSample)

wassComparison = double(length = p)
for(i in 1:p) {
  wassComparison[i] = wasserstein1d(a = betaSUN[i,], b = betaSUNcomparison[i,], p = 1)
}

# save Wasserstein distances
wass = cbind("MF-VB"=wassMF,"PFM-VB"=wassPFM,"MonteCarloError"=wassComparison)
save(wass, nSample, file="wassersteinDistances.RData") # save also number of samples used to get such distances
```

Moments
-------

Finally, we **compute and save the first two marginal moments** of the exact posterior distribution (via Monte-Carlo) and those provided by the two approximate methods (stored in `paramsPMF` and `paramsMF`, respectively).

``` r
meanSUN = apply(betaSUN,1,mean)
varSUN = apply(betaSUN,1,var)
moments = cbind(meanSUN=meanSUN,varSUN=varSUN,
                meanMF=paramsMF$mean,varMF=paramsMF$diagV,
                meanPFM=paramsPFM$postMoments.meanBeta,varPFM=paramsPFM$postMoments.varBeta)
save(moments, file = "moments.RData")
```

Plotting the results
====================

Let us first clean the workspace and load the required packages.

``` r
rm(list = ls())
library(ggplot2)
library(reshape)
library(RColorBrewer)
```

The code to reproduce each of the Figures 2, S2 and S3 in the paper is reported in the following sub-sections, together with the necessary steps to reproduce Figure 3 and S4.
 The four Figures can be found in the folder [`img`](https://github.com/augustofasano/Probit-PFMVB/tree/master/img).

Figure 2: comparison of moments and predictive probabilities
------------------------------------------------------------

The following code reproduces **Figure 2**.

``` r
load("moments.RData")
load("predProb.RData")

means = as.matrix(moments[,c(3,5)])
colnames(means) = c("MF","PFM")
meanData = melt(means)[,-1]
meanData$group1 = c("Expectations")
meanData$y = c(moments[,1],moments[,1])

variances = as.matrix(sqrt(moments[,c(4,6)]))
colnames(variances) = c("MF","PFM")
varData = melt(variances)[,-1]
varData$group1 = c("Standard Deviations")
varData$y = c(sqrt(moments[,2]),sqrt(moments[,2]))

pred = as.matrix(predProb[,c(2,3)])
colnames(pred) = c("MF","PFM")
predData = melt(pred)[,-1]
predData$group1 = c("Predictive Probabilities")
predData$y = c(predProb[,1],predProb[,1])

data_points = rbind(meanData,varData,predData)
data_points$group1 = factor(data_points$group1, levels=c("Expectations", "Standard Deviations","Predictive Probabilities" ))
F2 = ggplot(data_points, aes(y=value, x=y,color=X2)) +
  geom_point(aes(shape=X2),alpha=0.7)+facet_wrap(group1~., scales="free")+
  geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2)+theme_bw()+
  labs(x="Quantities computed from the exact posterior", y = "Quantities computed from the approximate posterior")+
  theme(axis.title.x = element_text(size=9),axis.title.y =element_text(size=9),plot.margin = margin(0.1, 0.1, 0.05, 0.2, "cm"),legend.position = "none")+
  scale_color_manual(values=c("#BDBDBD","#525252"))+
  scale_shape_manual(values=c(1, 2))

ggsave("F2.png",width=8,height=3)
```

Figures 3 and S4: comparison of moments and predictive probabilities for different *ν*<sub>*p*</sub><sup>2</sup>
------------------------------------------------------------

**Figures 3 and S4** show the same quantities reported in Figure 2, for two different scenarios where *ν*<sub>*p*</sub><sup>2</sup> is allowed to vary with *p*, respectively *ν*<sub>*p*</sub><sup>2</sup>*=2500/p* and *ν*<sub>*p*</sub><sup>2</sup>*=250/p*.

In order to study such scenarios and obtain Figure 3 or S4, it is sufficient to appropriately set *ν*<sub>*p*</sub><sup>2</sup> to the desired value at the beginning and then rerun all the code up to the *Moments* Section, excluding the Section *Wasserstein distances*.

Finally, after loading the required packages, the code used to produce Figure 2 will now output the figure corresponding to the considered scenario.

Figure S2: log-Wasserstein distances comparison
----------------------------------------------

The following code reproduces **Figure S2**. As preliminary step, it also checks the percentage of Wasserstein distances for each of the variational methods which lay within the 2.5%
and 97.5% of the Wasserstein distances between the two different samples of 20000 draws from the same exact posterior marginals quantiles, allowing in principle to assess how the observed variability compares with Monte Carlo error.

``` r
load("wassersteinDistances.RData")
head(wass)

# compare variational Wasserstein distances with Monte Carlo error
up = quantile(wass[,3],probs=0.975)
low = quantile(wass[,3],probs=0.025)

# percentage of PFM wasserstein distances falling into the Monte Carlo error bounds
sum(wass[,2]>low & wass[,2]<up)/dim(wass)[1]*100
# [1] 94.22311

# percentage of MF wasserstein distances falling into the Monte Carlo error bounds
sum(wass[,1]>low & wass[,1]<up)/dim(wass)[1]*100
# [1] 15.94732

# plot histograms Figure S2

errorM = log(wass[,3])
wass = wass[,-3]

up = quantile(errorM,probs=0.975)
low = quantile(errorM,probs=0.025)

data_hist = melt(wass)[,-1]
head(data_hist)
myColors = c(RColorBrewer::brewer.pal(9, "Greys")[3],RColorBrewer::brewer.pal(9, "Greys")[6])

FS2 = ggplot(data_hist, aes(x=log(value))) + geom_histogram(position="identity", alpha=0.5,bins=30,fill=myColors[1],color=myColors[2])+ theme_bw()+facet_grid(.~X2)+labs(x="Logarithm of the Wasserstein distance from the exact posterior", y = "")+theme(axis.title.x = element_text(size=9),axis.title.y =element_text(size=9),plot.margin = margin(0.1, 0.1, 0.1, -0.3, "cm"))+geom_vline(aes(xintercept=c(up)), lty = 2)+geom_vline(aes(xintercept=c(low)), lty = 2)

ggsave("FS2.png",width=7,height=3)
```

Figure S3: Best and worst case scenarios of each approximate method
---------------------------------------------------------------------------------------

The following code reproduces **Figure S3**.

``` r
load("BestWorstScenarios.RData")

bestMF = listBestWorst$bestMF
bestMF = melt(bestMF)[,-1]
bestMF$group1 = c("MF-VB")
bestMF$group2 = c("Best Marginal Approximation")
bestMF$group3 = c(rep(0.5,nSample),rep(0,2*nSample))

worstMF = listBestWorst$worstMF
worstMF = melt(worstMF)[,-1]
worstMF$group1 = c("MF-VB")
worstMF$group2 = c("Worst Marginal Approximation")
worstMF$group3 = c(rep(0.5,nSample),rep(0,2*nSample))

bestPMF = listBestWorst$bestPFM
bestPMF = melt(bestPMF)[,-1]
bestPMF$group1 = c("PFM-VB")
bestPMF$group2 = c("Best Marginal Approximation")
bestPMF$group3 = c(rep(0.5,nSample),rep(0,2*nSample))

worstPMF = listBestWorst$worstPFM
worstPMF = melt(worstPMF)[,-1]
worstPMF$group1 = c("PFM-VB")
worstPMF$group2 = c("Worst Marginal Approximation")
worstPMF$group3 = c(rep(0.5,nSample),rep(0,2*nSample))

dataDensity = rbind(bestMF,worstMF,bestPMF,worstPMF)
dataDensity$X2 = factor(dataDensity$X2,levels=c("SUN","MF","PFM"))

FS3 = ggplot(dataDensity, aes(x=value,fill=X2,color=X2,linetype=X2))+geom_density(alpha=0.1)+facet_grid(group2~group1)+ theme_bw()+scale_fill_manual(values=c("#000000","#FFFFFF","#FFFFFF"))+scale_color_manual(values=c("#D9D9D9","#000000","#000000"))+scale_linetype_manual(values=c("solid","dotted","longdash"))+labs(x="", y = "")+theme(axis.title.x = element_text(size=9),axis.title.y =element_text(size=9),plot.margin = margin(0.1, 0.1, -0.15, -0.3, "cm"),legend.position = "none")+xlim(-30,30)

ggsave("FS3.png",width=7,height=4)
```

Deviances with five-fold cross-validation
------------------------------

In this section, we show how to **compute the deviances** reported in Table 2 for the Alzheimer dataset.

Preliminary operations
------------------------------

Once the files `Alzheimer_Interactions.RData` and `functionsVariational.R` have been downloaded, set the working directory to the folder where they are located. Then, clean the workspace and load `Alzheimer_Interactions.RData` along with the source file `functionsVariational.R` and other useful packages.

After such preliminary opeations, the observations are **allocated  to 5 folds** to be used for cross-validation.

``` r
rm(list = ls())
library(gdata)
library(arm)
source("functionsVariational.R")
library("TruncatedNormal")
library("truncnorm")
library("transport")
library("mvtnorm")
library("sparsevb")
load("Alzheimer_Interactions.RData")

p = dim(X)[2] # number of covariates
nu2 = 25 # prior variance

X_raw = X[,-1]
y_data = y

n_dataset = length(y_data)

# Indicators five fold CV
set.seed(12)
sel_vector = c(rep(1,round(n_dataset/5)),rep(2,round(n_dataset/5)),rep(3,round(n_dataset/5)),rep(4,round(n_dataset/5)),rep(5,n_dataset-round(n_dataset/5)*4))
length(sel_vector)-n_dataset
sel_vector = sample(sel_vector,length(sel_vector))
```

Summary of the results
------------------------------

At this point, we can **loop over the 5 folds**.
At each iteration we take the observations allocated to the corresponding fold as test set and the remaining ones as training set.
The observations in the training set are standardized to have mean 0 and standard deviation 0.5, and the same scaling transformation is applied to the test set.
Then, the training set is used to get the approximate posterior distributions under **PFM-VB**, **MF-VB** and **SVB**. From these, the **predictive probabilities** for the test observations are computed.
Finally, these predictive probabilities are used to **compute the deviances of the test observations**, for each fold.

``` r
Dev_PFM = rep(0,5)
Dev_MF = rep(0,5)
Dev_SVB_gaussian = rep(0,5)


seed = 1
set.seed(seed)

for (f in 1:5){
  sel_set = which(sel_vector!=f)

  # Rescale
  X_data = X_raw

  for (j in 1:dim(X_data)[2]){
    X_data[,j] = (X_data[,j]-mean(X_data[sel_set,j]))/(2*sd(X_data[sel_set,j]))
  }
  X_data = cbind(rep(1,dim(X_data)[1]),X_data)


  # Training data
  y = y_data[sel_set]
  X = X_data[sel_set,]

  # Test data
  y_new = y_data[-sel_set]
  X_new = X_data[-sel_set,]

  n = dim(X)[1] # training set size

  # split training and test sets
  X_Test = X_new
  yTest = y_new

  # precompute some useful quantities to be used for the predictive probabilities
  VXt = t(nu2*X)%*%solve(diag(n)+(nu2*X)%*%t(X))
  invIXXt = solve(diag(1,nrow=n,ncol=n)+nu2*(X%*%t(X)))

  tolerance = 1e-3 # tolerance to establish ELBO convergence

  # get optimal parameters and moments PFM
  startTime = Sys.time()
  paramsPFM = getParamsPFM(X=X,y=y,nu2=nu2,moments=TRUE,tolerance=tolerance,maxIter=1e4)
  timeSUN_PFM_algo = difftime(Sys.time(), startTime, units=("secs"))[[1]]

  # get the predictive probabilities
  startTime = Sys.time()
  nSampleZ = 5e3
  muTN = paramsPFM$mu
  muTN[y==0] = -muTN[y==0]
  sampleTruncNorm = matrix(rtruncnorm(n*nSampleZ, a = 0, b = Inf, mean = muTN, sd = sqrt(paramsPFM$sigma2)), nrow = n, ncol = nSampleZ, byrow = F )
  sampleTruncNorm[y==0,] = -sampleTruncNorm[y==0,] # need to adjust the sign of the variables for which y_i is 0

  predPFM = double(length = length(yTest))
  for(i in 1:length(yTest)){
    xNew = matrix(X_Test[i,],ncol = 1)
    Xx = X%*%xNew
    sd = as.double(sqrt(1+nu2*(sum(xNew^2)-nu2*t(Xx)%*%invIXXt%*%Xx)))

    predPFM[i] = mean(pnorm((t(xNew)%*%VXt%*%sampleTruncNorm)/sd))}

  timeSUN_PFM_inference = difftime(Sys.time(), startTime, units=("secs"))[[1]]

  # get optimal parameters and moments MF
  startTime = Sys.time()
  paramsMF = getParamsMF(X,y,nu2,tolerance,maxIter = 1e4)
  timeMF_algo = difftime(Sys.time(), startTime, units=("secs"))[[1]]

  # get the predictive probabilities
  startTime = Sys.time()
  predMF = double(length = length(yTest))
  for(i in 1:length(yTest)){
    xNew = matrix(X_Test[i,],ncol = 1)
    Xx = X%*%xNew
    sd = as.double(sqrt(1+nu2*(sum(xNew^2)-nu2*t(Xx)%*%invIXXt%*%Xx)))

    predMF[i] = as.double(pnorm(t(xNew)%*%paramsMF$meanBeta/sd))}

  timeMF_inference = difftime(Sys.time(), startTime, units=("secs"))[[1]]


  # get optimal parameters and moments SVB (intercept = TRUE)
  startTime = Sys.time()
  paramsSVB_gaussian = svb.fit(X = X[,-1],
                               Y = y,
                               family = "logistic",
                               slab = "gaussian",
                               sigma = rep(sqrt(nu2), p-1),
                               prior_scale = sqrt(nu2),
                               intercept = TRUE,
                               max_iter = 1e4,
                               tol = tolerance)
  timeSVB_algo = difftime(Sys.time(), startTime, units=("secs"))[[1]]

  # predictive probabilities
  nSample = 5e3 # number of samples to compute predictive probabilities in LOGIT SVB, coherent with PFM
  startTime = Sys.time()
  sampleBetaSVB_gaussian = sampleSVB(paramsSVB_gaussian,nSample)
  sampleBetaSVB_gaussian = rbind(rep(paramsSVB_gaussian$intercept,nSample),sampleBetaSVB_gaussian)
  predSVB_gaussian = double(length = length(yTest))
  for(i in 1:length(yTest)){
    xNew = as.matrix(X_Test[i,])
    linPred = t(xNew)%*%sampleBetaSVB_gaussian

    predSVB_gaussian[i] = mean(plogis(linPred))
  }
  timeSVB_inference = difftime(Sys.time(), startTime, units=("secs"))[[1]]

  Dev_PFM[f] = c(-sum(log((predPFM^yTest)*((1-predPFM)^(1-yTest)))))
  Dev_MF[f] = c(-sum(log((predMF^yTest)*((1-predMF)^(1-yTest)))))
  Dev_SVB_gaussian[f] = c(-sum(log((predSVB_gaussian^yTest)*((1-predSVB_gaussian)^(1-yTest)))))

  print(f)
}
```

Summing the deviances obtained in the 5 folds we get the desired result.

``` r
sum(Dev_PFM)
#[1] 187.5177
sum(Dev_MF)
#[1] 228.7083
sum(Dev_SVB_gaussian)
#[1]  126.2379
```

The lesion application
==================

In this section, we show instead how to **compute the deviances** reported in Table 2 for the lesion dataset.
The procedure is analogous to the one shown above for the Alzheimer dataset.
The original dataset is available at the [UCI repository](https://archive.ics.uci.edu/ml/datasets/Gastrointestinal+Lesions+in+Regular+Colonoscopy).
The dataset has been pre-processed to have one row for each patient, combining in one vector the measurements obtained with the two different lights, as suggested on the dataset website.
Moreover, columns having more than 95% of the observations equal to zero have been removed.
The resulting dataset is available in the [`data`](https://github.com/augustofasano/Probit-PFMVB/tree/master/data) folder and is called `lesion.RData`.

Preliminary operations
------------------------------

Once the files `lesion.RData` and `functionsVariational.R` have been downloaded, set the working directory to the folder where they are located. Then, clean the workspace and load `lesion.RData` along with the source file `functionsVariational.R` and other useful packages.

``` r
rm(list = ls())
library(gdata)
library(arm)
source("functionsVariational.R")
library("TruncatedNormal")
library("truncnorm")
library("transport")
library("mvtnorm")
library("sparsevb")
load("lesion.RData")
```

The matrix of raw covariates `X_raw` and the vector of binary observations `y_data` are loaded.
Now, we **get the number of observations** `n_dataset` and **allocate them to 5 folds** to be used for cross-validation.

``` r
p = dim(X_raw)[2] + 1 # number of covariates (including intecept)
nu2 = 25*100/p # prior variance

n_dataset = length(y_data)

# Indicators five fold CV
set.seed(12)
sel_vector = c(rep(1,round(n_dataset/5)),rep(2,round(n_dataset/5)),rep(3,round(n_dataset/5)),rep(4,round(n_dataset/5)),rep(5,n_dataset-round(n_dataset/5)*4))
length(sel_vector)-n_dataset
sel_vector = sample(sel_vector,length(sel_vector))
```

Summary of the results
------------------------------

At this point, we can **loop over the 5 folds**.
At each iteration we take the observations allocated to the corresponding fold as test set and the remaining ones as training set.
The observations in the training set are standardized to have mean 0 and standard deviation 0.5, and the same scaling transformation is applied to the test set.
Then, the training set is used to get the approximate posterior distributions under both **PFM-VB**, **MF-VB** and **SVB**. From these, the **predictive probabilities** for the test observations are computed, as in the previous Alzheimer's application.
Finally, these predictive probabilities are used to **compute the deviances of the test observations**, for each fold.

``` r
Dev_PFM = rep(0,5)
Dev_MF = rep(0,5)
Dev_SVB_gaussian = rep(0,5)

seed = 1
set.seed(seed)

for (f in 1:5){
  sel_set = which(sel_vector!=f)

  # Rescale
  X_data = X_raw

  for (j in 1:dim(X_data)[2]){
    X_data[,j] = (X_data[,j]-mean(X_data[sel_set,j]))/(2*sd(X_data[sel_set,j]))
  }
  X_data = cbind(rep(1,dim(X_data)[1]),X_data)


  # Training data
  y = y_data[sel_set]
  X = X_data[sel_set,]

  # Test data
  y_new = y_data[-sel_set]
  X_new = X_data[-sel_set,]

  n = dim(X)[1] # training set size

  # split training and test sets
  X_Test = X_new
  yTest = y_new

  # precompute some useful quantities to be used for the predictive probabilities
  VXt = t(nu2*X)%*%solve(diag(n)+(nu2*X)%*%t(X))
  invIXXt = solve(diag(1,nrow=n,ncol=n)+nu2*(X%*%t(X)))

  tolerance = 1e-3 # tolerance to establish ELBO convergence

  # get optimal parameters and moments PFM
  startTime = Sys.time()
  paramsPFM = getParamsPFM(X=X,y=y,nu2=nu2,moments=TRUE,tolerance=tolerance,maxIter=1e4)
  timeSUN_PFM_algo = difftime(Sys.time(), startTime, units=("secs"))[[1]]

  # get the predictive probabilities
  startTime = Sys.time()
  nSampleZ = 5e3
  muTN = paramsPFM$mu
  muTN[y==0] = -muTN[y==0]
  sampleTruncNorm = matrix(rtruncnorm(n*nSampleZ, a = 0, b = Inf, mean = muTN, sd = sqrt(paramsPFM$sigma2)), nrow = n, ncol = nSampleZ, byrow = F )
  sampleTruncNorm[y==0,] = -sampleTruncNorm[y==0,] # need to adjust the sign of the variables for which y_i is 0

  predPFM = double(length = length(yTest))
  for(i in 1:length(yTest)){
    xNew = matrix(X_Test[i,],ncol = 1)
    Xx = X%*%xNew
    sd = as.double(sqrt(1+nu2*(sum(xNew^2)-nu2*t(Xx)%*%invIXXt%*%Xx)))

    predPFM[i] = mean(pnorm((t(xNew)%*%VXt%*%sampleTruncNorm)/sd))}

  timeSUN_PFM_inference = difftime(Sys.time(), startTime, units=("secs"))[[1]]



  # get optimal parameters and moments MF
  startTime = Sys.time()
  paramsMF = getParamsMF(X,y,nu2,tolerance,maxIter = 1e4)
  timeMF_algo = difftime(Sys.time(), startTime, units=("secs"))[[1]]

  # get the predictive probabilities
  startTime = Sys.time()
  predMF = double(length = length(yTest))
  for(i in 1:length(yTest)){
    xNew = matrix(X_Test[i,],ncol = 1)
    Xx = X%*%xNew
    sd = as.double(sqrt(1+nu2*(sum(xNew^2)-nu2*t(Xx)%*%invIXXt%*%Xx)))

    predMF[i] = as.double(pnorm(t(xNew)%*%paramsMF$meanBeta/sd))}

  timeMF_inference = difftime(Sys.time(), startTime, units=("secs"))[[1]]



  # get optimal parameters and moments SVB (intercept = TRUE)
  startTime = Sys.time()
  paramsSVB_gaussian = svb.fit(X = X[,-1],
                               Y = y,
                               family = "logistic",
                               slab = "gaussian",
                               sigma = rep(sqrt(nu2), p-1),
                               prior_scale = sqrt(nu2),
                               intercept = TRUE,
                               max_iter = 1e4,
                               tol = tolerance)
  timeSVB_algo = difftime(Sys.time(), startTime, units=("secs"))[[1]]

  # predictive probabilities
  nSample = 5e3 # number of samples to compute predictive probabilities in LOGIT SVB, coherent with PFM
  startTime = Sys.time()
  sampleBetaSVB_gaussian = sampleSVB(paramsSVB_gaussian,nSample)
  sampleBetaSVB_gaussian = rbind(rep(paramsSVB_gaussian$intercept,nSample),sampleBetaSVB_gaussian)
  predSVB_gaussian = double(length = length(yTest))
  for(i in 1:length(yTest)){
    xNew = as.matrix(X_Test[i,])
    linPred = t(xNew)%*%sampleBetaSVB_gaussian

    predSVB_gaussian[i] = mean(plogis(linPred))
  }
  timeSVB_inference = difftime(Sys.time(), startTime, units=("secs"))[[1]]


  Dev_PFM[f] = c(-sum(log((predPFM^yTest)*((1-predPFM)^(1-yTest)))))
  Dev_MF[f] = c(-sum(log((predMF^yTest)*((1-predMF)^(1-yTest)))))
  Dev_SVB_gaussian[f] = c(-sum(log((predSVB_gaussian^yTest)*((1-predSVB_gaussian)^(1-yTest)))))

  print(f)
}
```

Summing the deviances obtained in the 5 folds we get the desired result.

``` r
sum(Dev_PFM)
#[1] 27.24379
sum(Dev_MF)
#[1] 48.65794
sum(Dev_SVB_gaussian)
# [1] 35.80594
```
