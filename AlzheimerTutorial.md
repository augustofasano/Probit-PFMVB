Introduction
============

As described in the README.md file, this tutorial contains general
guidelines and code to **perform the analyses for the Alzheimer
application in Section 3** of the paper. In particular, you will find
information on how to download the data, detailed <tt>R</tt> code to
**implement the different methods for posterior inference** discussed in
Section 2 and guidelines to produce Figures 1 to 3 in the paper. For
implementation purposes, execute the code below considering the same
order in which is presented.

The Alzheimer dataset
=====================

As discussed in Section 3, the focus is to **model presence or absence
of Alzheimer’s disease in its early stages as a function of demographic
data, genotype and assay results**. The original dataset is available in
the <tt>R</tt> library <tt>AppliedPredictiveModeling</tt> and arises
from a study of the Washington University to determine if biological
measurements from cerebrospinal fluid are useful in modeling and
predicting early stages of Alzheimer’s disease ([Craig-Schapiro et al.,
2011](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0018850)).
Here, we **avoid excessively complex black-box algorithms** and rely on
an **interpretable probit regression**, which improves flexibility by
simply adding pairwise interactions, thus obtaining *p* = 9036
predictors collected for *n* = 333 individuals. Following [Gelman et al.
(2008)](https://projecteuclid.org/euclid.aoas/1231424214) and [Chopin
and Ridgway (2017)](https://projecteuclid.org/euclid.ss/1491465628) the
original measurements have been standardized to have mean 0 and standard
deviation 0.5, before entering such variables and their interactions in
the probit regression.

The final dataset <tt>Alzheimer.RData</tt> is available in the
<tt>data</tt> folder of this repository.

Preliminary operations
======================

Once the files <tt>Alzheimer.RData</tt> and
<tt>functionsVariational.R</tt> have been downloaded, set the working
directory to the folder they have been saved in. Then, clean the
workspace and load the data, source the file
<tt>functionsVariational.R</tt> and load other useful packages.

    rm(list=ls())
    source("functionsVariational.R")
    library("TruncatedNormal")
    library("truncnorm")
    library("transport")

    load("Alzheimer_Interactions.RData")

Now, we take *n* = 300 observations to perform posterior inference,
keeping the remaining 33 observations for out-of-sample analysis.

    seed = 1
    set.seed(seed)
    n = 300 # training set size
    trainingSet = sample(x = dim(X)[1],size = n, replace = FALSE)

    # split training and test sets
    X_Test = X[-trainingSet,]
    yTest = y[-trainingSet]
    X = X[trainingSet,]
    y = y[trainingSet]

Finally, we conclude the preliminary part by getting the model
dimension, specifying the prior variance (in accordance with [Gelman et
al. (2008)](https://projecteuclid.org/euclid.aoas/1231424214)
guidelines) and the number of i.i.d samples to generate from the exact
and approximate posterior densities.

    p = dim(X)[2] # get number of covariates
    nu2 = 25 # prior variance
    nSample = 2e4 #fix number of samples

Exact posterior sampling
========================

After the preliminary operations above have been executed,
**i.i.d. samples from the unified skew-normal posterior** ([Durante,
2019](https://arxiv.org/abs/1802.09565)) can be obtained by simply using
the function <tt>rSUNpost</tt>. In addition, note we also keep track of
the **running time**.

    startTimeSUN = Sys.time()
    betaSUN = rSUNpost(X=X,y=y,nu2=nu2,nSample=nSample)
    endTimeSUN = Sys.time()
    timeSUN = difftime(endTimeSUN, startTimeSUN, units=("mins"))[[1]]

Partially-factorized mean-field variational Bayes
=================================================

This Section contains the code for the following two tasks:

-   **get the optimal parameters of**
    *q*<sub>PFM</sub><sup>\*</sup>(**z**) as in **Algorithm 2**, with
    the function <tt>getParamsPFM</tt>. Then, the optimal joint <span
    class="smallcaps">PFM-VB</span> approximating density immediately
    follows from the equality
    *q*<sub>PFM</sub><sup>\*</sup>(**β**, **z**) = *p*(**β** ∣ **z**)*q*<sub>PFM</sub><sup>\*</sup>(**z**);

-   **sample from the optimal approximating posterior distribution** for
    **β**, *q*<sub>PFM</sub><sup>\*</sup>(**β**): this is performed
    through the function <tt>sampleSUN\_PFM</tt>.

Again, the **running time** is also monitored.

    tolerance = 1e-3 # tolerance to establish ELBO convergence

    # get optimal parameters and moments
    startTimeSUN_PFM = Sys.time()
    paramsPFM = getParamsPFM(X=X,y=y,nu2=nu2,moments=TRUE,tolerance=tolerance,maxIter=1e4)
    endTimeSUN_PFM = Sys.time()
    timeSUN_PFM = difftime(endTimeSUN_PFM, startTimeSUN_PFM, units=("secs"))[[1]]

    # sample from approximate posterior
    set.seed(seed+2)
    betaSUN_PFM = sampleSUN_PFM(paramsPFM=paramsPFM,X=X,y=y,nu2=nu2,nSample=nSample)

Mean-field variational Bayes
============================

This Section contains the code for the following two tasks:

-   **get the optimal parameters of**
    *q*<sub>MF</sub><sup>\*</sup>(**β**) as in **Algorithm 1**, with the
    function <tt>getParamsMF</tt>;

-   **sample from the optimal approximating posterior distribution** for
    **β**, *q*<sub>MF</sub><sup>\*</sup>(**β**): since this is a
    multivariate normal distribution, this is performed in a standard
    way, by exploiting the Choleski decomposition of the optimal
    covariance matrix **V**.

As usual, the **running time** is also monitored.

    startTimeMF = Sys.time()
    paramsMF = getParamsMF(X,y,nu2,tolerance,maxIter = 1e4)
    endTimeMF = Sys.time()
    timeMF = difftime(endTimeMF, startTimeMF, units=("secs"))[[1]]

    # sample from approximate posterior
    set.seed(seed+3)
    L = t(chol(paramsMF$V))
    betaMF = L%*%matrix(rnorm(nSample*p),p,nSample)
    betaMF = apply(betaMF,2, function(x) paramsMF$meanBeta+x)

Running times and number of iterations
======================================

At this point, we can **check the running times** for the three methods
and **number of iterations** of each approximate method with the
following simple code.

**Remark**: depending on the machine you use, the running times can
slightly differ from the ones reported in the paper, but they should not
significantly deviate from them. On the contrary, the number of
iterations for the two methods should match the reported ones, i.e. 212
for <span class="smallcaps">MF-VB</span> and 14 for <span
class="smallcaps">PFM-VB</span>

    timeSUN
    timeMF
    timeSUN_PFM

    paramsMF$nIter
    paramsPFM$nIter

Predictive probabilities
========================

In this section we report the code to compute the **predictive
probabilities** according to the exact unified skew-normal, the <span
class="smallcaps">PFM-VB</span> and the <span
class="smallcaps">MF-VB</span> posterior distributions. See Section 2 of
the paper for further details about how these probabilities can be
computed and the <tt>functionsTutorial.html</tt> for their
implementation.

    # Compute V and VXt outside, using Woodbury's identity
    Omega = diag(rep(nu2,p),p,p)
    invOmega = diag(rep(1/nu2,p),p,p)
    VXt = t(nu2*X)%*%solve(diag(n)+(nu2*X)%*%t(X))
    V = Omega - VXt%*%(nu2*X)

    nTest = length(yTest)
    predProb = matrix(nrow = nTest, ncol = 3)
    colnames(predProb) = c("SUN", "MF", "PFM")
    for(i in 1:nTest){
      xNew = matrix(X_Test[i,],ncol = 1)
      predProb[i,] = c(predictMC(xNew,betaSUN),
                       predictMF(xNew,paramsMF),
                       predictPFM(xNew,paramsPFM,X,y,1e4))
    }

    predProb = cbind(predProb, y=yTest)

    save(predProb, file="predProb.RData")

Summaries of the results
========================

In this Section, we compute the different **summary quantities** of our
results which have been used to produce the figures in Section 3 of the
paper.

Wasserstein distances
---------------------

We start by computing the **Wasserstein distances** among the *p* = 9036
posterior empirical marginal distributions induced by the 20000
i.i.d. samples from the exact posterior and the associated ones obtained
with the 20000 i.i.d. samples under <span class="smallcaps">MF-VB</span>
and <span class="smallcaps">PFM-VB</span> approximations.

We compute such quantities with the function <tt>wasserstein1d</tt>
available in the <tt>transport</tt> package, which has been loaded at
the beginning.

    wassMF = double(length = p)
    for(i in 1:p) {
      wassMF[i] = wasserstein1d(a = betaSUN[i,], b = betaMF[i,], p = 1)
    }

    wassPFM = double(length = p)
    for(i in 1:p) {
      wassPFM[i] = wasserstein1d(a = betaSUN[i,], b = betaSUN_PFM[i,], p = 1)
    }

Now, we find the indexes of the parameters for which we have the **best
and worst approximations**, for both the methods <span
class="smallcaps">MF-VB</span> and <span
class="smallcaps">PFM-VB</span>. Then, we **save the samples obtained
for these four parameters** when sampling from the exact posterior and
from the two approximate methods.

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

Since these distances are computed among couples of empirical
distributions obtained with 20000 samples, and are affected by the
**variability due to the Monte-Carlo sampling**, we compute also the
Wasserstein distances of the marginals one would obtain when **comparing
two samples from the exact unified skew-normal distribution**.

Hence, we generate other 20000 i.i.d. samples from the true posterior
and compute the Wasserstein distances among the marginals obtained in
the two samples from the exact posterior, in order to have a benchmark
and **assess what amount in the variability of the Wasserstein distances
is due to Monte-Carlo error**.

All the Wasserstein distances are then saved.

    rm(list=c("betaSUN_PFM","betaMF")) # free memory
    set.seed(seed+100)
    betaSUNcomparison = rSUNpost(X,y,nu2,nSample = nSample)

    wassComparison = double(length = p)
    for(i in 1:p) {
      wassComparison[i] = wasserstein1d(a = betaSUN[i,], b = betaSUNcomparison[i,], p = 1)
    }

    # save Wasserstein distances
    wass = cbind("MF-VB"=wassMF,"PFM-VB"=wassPFM,"MonteCarloError"=wassComparison)
    save(wass, nSample, file="wassersteinDistances.RData") # save also number of samples used to get such distances

Moments
-------

Finally, we **compute and save the first two marginal moments** for the
exact posterior distribution and the two approximate methods.

The moments of the exact posterior distribution are simply obtained
empirically from the i.i.d. samples, while the ones of the approximate
methods are directly available in closed-form and are accessible among
the quantities stored in the <tt>paramsMF</tt> and <tt>paramsPFM</tt>
variables.

    meanSUN = apply(betaSUN,1,mean)
    varSUN = apply(betaSUN,1,var)
    moments = cbind(meanSUN=meanSUN,varSUN=varSUN,
                    meanMF=paramsMF$mean,varMF=diag(paramsMF$V),
                    meanPFM=paramsPFM$postMoments.meanBeta,varPFM=paramsPFM$postMoments.varBeta)
    save(moments, file = "moments.RData")

Plotting the results
====================

At this point, **we have computed all the quantities we need**, so we
clean the workspace and load the needed packages.

The code to reproduce each of the figures in the paper is reported in
the following sub-sections.

    rm(list = ls())
    library(ggplot2)
    library(reshape)
    library(RColorBrewer)

Figure 1: log-Wasserstein distances comparison
----------------------------------------------

The following code reproduces Figure 1.

    load("wassersteinDistances.RData")
    head(wass)
    errorM = log(wass[,3])
    wass = wass[,-3]

    up = quantile(errorM,probs=0.975)
    low = quantile(errorM,probs=0.025)

    data_hist = melt(wass)[,-1]
    head(data_hist)
    myColors = c(RColorBrewer::brewer.pal(9, "Greys")[3],RColorBrewer::brewer.pal(9, "Greys")[6])

    F1 = ggplot(data_hist, aes(x=log(value))) + geom_histogram(position="identity", alpha=0.5,bins=30,fill=myColors[1],color=myColors[2])+ theme_bw()+facet_grid(.~X2)+labs(x="Logarithm of the Wasserstein distance from the exact posterior", y = "")+theme(axis.title.x = element_text(size=9),axis.title.y =element_text(size=9),plot.margin = margin(0.1, 0.1, 0.1, -0.3, "cm"))+geom_vline(aes(xintercept=c(up)), lty = 2)+geom_vline(aes(xintercept=c(low)), lty = 2)

    ggsave("F1.png",width=7,height=3)

![](https://raw.githubusercontent.com/augustofasano/Probit-PFMVB/master/img/F1.png)

Figure 2: density comparison in the best and worst scenarios of each approximate method
---------------------------------------------------------------------------------------

The following code reproduces Figure 2.

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
    F2 = ggplot(dataDensity, aes(x=value,fill=X2,color=X2,linetype=X2))+
      geom_density(alpha=0.1)+facet_grid(group2~group1)+ theme_bw()+scale_fill_manual(values=c("#000000","#FFFFFF","#FFFFFF"))+scale_color_manual(values=c("#D9D9D9","#000000","#000000"))+scale_linetype_manual(values=c("solid","dotted","longdash"))+labs(x="", y = "")+theme(axis.title.x = element_text(size=9),axis.title.y =element_text(size=9),plot.margin = margin(0.1, 0.1, -0.15, -0.3, "cm"),legend.position = "none")+xlim(-29,29)

    ggsave("F2.png",width=7,height=4)

![](https://raw.githubusercontent.com/augustofasano/Probit-PFMVB/master/img/F2.png)

Figure 3: comparison of moments and predictive probabilities
------------------------------------------------------------

The following code reproduces Figure 3.

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
    F3 = ggplot(data_points, aes(y=value, x=y,color=X2)) +
      geom_point(aes(shape=X2),alpha=0.7)+facet_wrap(group1~., scales="free")+  geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2)+ theme_bw()+labs(x="Quantities computed from the exact posterior", y = "Quantities computed from the approximate posterior")+theme(axis.title.x = element_text(size=9),axis.title.y =element_text(size=9),plot.margin = margin(0.1, 0.1, 0.05, 0.2, "cm"),legend.position = "none")+scale_color_manual(values=c("#BDBDBD","#525252"))+
      scale_shape_manual(values=c(1, 2))

    ggsave("F3.png",width=8,height=3)

![](https://raw.githubusercontent.com/augustofasano/Probit-PFMVB/master/img/F3.png)
