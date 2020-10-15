Introduction
============
As described in the [`README.md`](https://github.com/augustofasano/Probit-PFMVB/blob/master/README.md) file, this tutorial provides details on the **functions required to implement the methods presented in Section 2 of the paper**. The `R` source file can be found in [`functionsVariational.R`](https://github.com/augustofasano/Probit-PFMVB/blob/master/functionsVariational.R).

Implemented functions
---------------------
The **list of implemented functions** is reported below. Each of them is then analyzed in detail in the following.

-   `getParamsPFM`: returns the parameters of the optimal PFM approximation (**Algorithm 2** in the paper)
-   `sampleSUN_PFM`: samples from the optimal PFM approximating density (**Algorithm 3** in the paper)
-   `getParamsMF`: returns the parameters of the optimal MF approximation (**Algorithm 1** in the paper)
-   `rSUNpost`: samples from the exact **SUN** posterior distribution ([Durante, 2019](https://doi.org/10.1093/biomet/asz034))

### `getParamsPFM`

This function implements the **CAVI** to obtain the optimal PFM approximating density. See **Algorithm 2** in the paper.

**Input**:

-   `X`: *n* × *p* matrix of explanatory variables
-   `y`: binary vector of response variables
-   `nu2`: prior variance for *β*<sub>*j*</sub>’s coefficients (*ν*<sub>*p*</sub><sup>2</sup> in the paper)
-   `moments`: logical, do you want to obtain marginal moments in the output?
-   `tolerance`: absolute change in the ELBO\[*q*<sub>PFM</sub>(**z**)\] used to establish convergence
-   `maxIter`: maximum number of allowed iterations before stopping

**Output**: A list containing

-   `mu`: optimal value **μ**<sup>\*</sup>
-   `sigma2`: optimal value **σ**<sup> \* 2</sup>
-   `nIter`: number of iterations required by the **CAVI**, either because it converged or because the maximum number of iterations `maxIter` was reached
-   (optional, if `moments` is set to TRUE) `postMoments`: list containing the posterior mean (`postMoments.meanBeta`) and marginal posterior variances (`postMoments.varBeta`) of **β**

**Remark**: The ELBO\[*q*<sub>PFM</sub>(**z**)\] = ELBO\[*q*<sub>PFM</sub>(**β**,**z**)\] expression is reported below.
For ease of reading, call
**Ω**<sub>*z*</sub> = **I**<sub>*n*</sub> + *ν*<sub>*p*</sub><sup>2</sup>**X** **X**<sup>⊺</sup>.
Notice that
*p*(**z**, **y**) = *p*(**z**)p(**y**|**z**) = *ϕ*<sub>*n*</sub>( **0**; **Ω**<sub>*z*</sub>)∏<sub>*i* = 1, …, *n*</sub> **1**\[*z*<sub>*i*</sub> > 0\], thus

ELBO\[*q*<sub>PFM</sub>(**z**)\] = C -
0.5∑<sub>*i* = 1, …, *n*</sub>\[(**Ω**<sub>*z*</sub><sup> − 1</sup>)<sub>*ii*</sub>*E*\[*z*<sub>*i*</sub><sup>2</sup>\] - *σ*<sub>*i*</sub><sup> \* − 2</sup>*E*\[*z*<sub>*i*</sub><sup>2</sup>\] + 2*z̄*<sub>*i*</sub><sup>(*t*)</sup>*μ*<sub>*i*</sub><sup>(*t*)</sup>/*σ*<sub>*i*</sub><sup> \* 2</sup> - (*μ*<sub>*i*</sub><sup>(*t*)</sup>/*σ*<sub>*i*</sub><sup>\*</sup>)<sup>2</sup> - 2log (*Φ*((2*y*<sub>*i*</sub> − 1)*μ*<sub>*i*</sub><sup>(*t*)</sup>/*σ*<sub>*i*</sub><sup>\*</sup>))\] - ∑<sub>*i>j*</sub>(**Ω**<sub>*z*</sub><sup> − 1</sup>)<sub>*ij*</sub>*z̄*<sub>*i*</sub><sup>(*t*)</sup>*z̄*<sub>*j*</sub><sup>(*t*)</sup>,

where *C* is a constant term and all the expectations are with respect
to *q*<sub>*PFM*</sub><sup>(*t*)</sup>(**z**), so that

*z̄*<sub>*i*</sub><sup>(*t*)</sup> = *E*\[*z*<sub>*i*</sub>\] = *μ*<sub>*i*</sub><sup>(*t*)</sup> +
(2*y*<sub>*i*</sub> − 1)*σ*<sub>*i*</sub><sup>\*</sup> *ϕ*(*μ*<sub>*i*</sub><sup>(*t*)</sup>/*σ*<sub>*i*</sub><sup>\*</sup>) \[*Φ*((2*y*<sub>*i*</sub> − 1)*μ*<sub>*i*</sub><sup>(*t*)</sup>/*σ*<sub>*i*</sub><sup>\*</sup>) \]<sup> − 1</sup>,

*E*\[*z*<sub>*i*</sub><sup>2</sup>\] = *μ*<sub>*i*</sub><sup>(*t*) 2</sup> + *σ*<sub>*i*</sub><sup>\* 2</sup> + (2*y*<sub>*i*</sub> − 1)*μ*<sub>*i*</sub><sup>(*t*)</sup>*σ*<sub>*i*</sub><sup>\*</sup>*ϕ*(*μ*<sub>*i*</sub><sup>(*t*)</sup>/*σ*<sub>*i*</sub><sup>\*</sup>) \[*Φ*((2*y*<sub>*i*</sub> − 1)*μ*<sub>*i*</sub><sup>(*t*)</sup>/*σ*<sub>*i*</sub><sup>\*</sup>) \]<sup> − 1</sup>.



``` r
getParamsPFM = function(X,y,nu2,moments = TRUE,tolerance = 1e-2, maxIter = 1e4) {
  ######################################################
  # PRECOMPUTATION
  ######################################################
  # define model dimensions
  n = dim(X)[1]
  p = dim(X)[2]

  # compute H = X%*%V%*%t(X) and Omega_z directly or with Woodbury
  if(p<=n) {
    # define prior covariance matrix and its inverse
    Omega = diag(rep(nu2,p),p,p)
    invOmega = diag(rep(1/nu2,p),p,p)
    V = solve(t(X)%*%X+invOmega)
    H = X%*%V%*%t(X)
    invOmZ = diag(1,nrow=n,ncol=n) - H # needed for ELBO
  } else{
    XXt = X%*%t(X)
    invOmZ = solve(diag(1,nrow=n,ncol=n)+nu2*XXt) # needed for ELBO
    H = nu2*XXt%*%invOmZ
  }

  # compute optimal sigma2
  h = diag(diag(H))
  sigma2 = matrix(1/(1-diag(H)), ncol = 1)
  sigma = sqrt(sigma2)

  # compute matrix to write the CAVI update in a vectorized form
  A = diag(as.double(sigma2), nrow = n, ncol = n)%*%(H - h)

  # other useful quantities needed for ELBO
  diagInvOmZ = diag(invOmZ)
  coeffMean_Z2 = diagInvOmZ-1/sigma2

  # initialization of variables
  mean_Z2 = matrix(0,n,1)
  mu = matrix(0,n,1)

  musiRatio = as.double(mu/sigma)
  phiPhiRatio = exp(dnorm(musiRatio,log = T)-pnorm((2*y-1)*musiRatio,log.p = T))
  meanZ = mu + (2*y-1)*sigma*phiPhiRatio

  elbo = -Inf
  diff = 1
  nIter=0

  ######################################################
  # CAVI ALGORITHM
  ######################################################

  while(diff > tolerance & nIter < maxIter) {
    elboOld = elbo
    sumLogPhi = 0

    for(i in 1:n) {
      mu[i] = A[i,]%*%meanZ

      # compute first (needed for algorithm) and second (needed for ELBO) moments
      musiRatio = mu[i]/sigma[i]
      phiPhiRatio = exp(dnorm(musiRatio, log = T) - pnorm((2*y[i]-1)*musiRatio, log.p = T))
      meanZ[i] = mu[i] + (2*y[i]-1)*sigma[i]*phiPhiRatio
      mean_Z2[i] = mu[i]^2+sigma2[i]+(2*y[i]-1)*mu[i]*sigma[i]*phiPhiRatio # needed for ELBO
      sumLogPhi = sumLogPhi + pnorm((2*y[i]-1)*musiRatio, log.p = T)
    }

    # computation of ELBO (up to an additive constant not depending on mu)
    elbo = -(t(meanZ)%*%invOmZ%*%meanZ -
               sum((meanZ^2)*diagInvOmZ) +
               sum(mean_Z2*coeffMean_Z2))/2 -
      sum(meanZ*mu/sigma2) + sum((mu^2)/sigma2)/2 + sumLogPhi

    diff = abs(elbo-elboOld)
    nIter = nIter+1

    if(nIter%%100==0) {print(paste0("iter: ", nIter, ", ELBO: ", elbo))}
  }

  # get the optimal parameters of the normals before truncation, now that convergence has been reached
  mu = A%*%meanZ

  results = list(mu = mu, sigma2 = sigma2, nIter = nIter)

  ######################################################
  # (OPTIONAL) CLOSED-FORM MOMENTS' COMPUTATION
  ######################################################

  if(moments == TRUE) {
    # compute V and V%*%t(X), directly or with Woodbury
    if(p<=n) {
      diagV = diag(V) # V already computed
      VXt = V%*%t(X)
    } else{ # use Woodbury
      VXt = t(nu2*X)%*%solve(diag(n)+(nu2*X)%*%t(X))
      #V = Omega - VXt%*%(nu2*X)
      diagV = nu2*(1-colSums(t(VXt) * X))
    }

    musiRatio = mu/sigma
    phiPhiRatio = exp(dnorm(musiRatio, log = T) - pnorm((2*y-1)*musiRatio, log.p = T))

    meanZ = mu + (2*y-1)*sigma*phiPhiRatio
    postVarZ = as.double(sigma2*(1-(2*y-1)*musiRatio*phiPhiRatio - phiPhiRatio^2))

    W = apply(VXt,1,function(x) sum(x*x*postVarZ))

    meanBeta = VXt%*%meanZ
    varBeta = diagV + W

    moments_PFM = list(meanBeta=meanBeta,varBeta=matrix(varBeta,ncol = 1))

    results = c(results,postMoments=moments_PFM)
  }

  return(results)
}
```

### `sampleSUN_PFM`

This function **samples from the optimal unified skew-normal PFM approximating density** for **β**. See **Algorithm 3** in the paper. Here, we implement an efficient version of such a routine which samples from the **joint** approximate density, thus allowing Monte Carlo inference on any complex functional of **β**. **Note** that, as discussed in the article, if the interest is only on the marginals of **β**, then sampling from each *β*<sub>*j*</sub> can be done at a much higher speed.

**Input**:

-   `paramsPFM`: output of the function `getParamsPFM`
-   `X`: *n* × *p* matrix of explanatory variables
-   `y`: binary vector of response variables
-   `nu2`: prior variance for *β*<sub>*j*</sub>’s coefficients (*ν*<sub>*p*</sub><sup>2</sup> in the paper)
-   `nSample`: number of i.i.d. samples from *q*<sup>\*</sup><sub>PFM</sub>(**β**) to generate

**Output**: A *p* × `nSample` matrix, where each column is a sample from *q*<sup>\*</sup><sub>PFM</sub>(**β**).

``` r
sampleSUN_PFM = function(paramsPFM, X, y, nu2, nSample) {
  # get model dimensions
  n = dim(X)[1]
  p = dim(X)[2]

  # get parameters useful for sampling
  muTN = paramsPFM$mu
  muTN[y==0] = -muTN[y==0] # generate all the truncated as left truncated, with support (0,Inf)

  Omega = diag(rep(nu2,p),p,p)
  invOmega = diag(rep(1/nu2,p),p,p)

  if(p<=n) { # compute V directly or with Woodbury
    V = solve(t(X)%*%X+invOmega)
    VXt = V%*%t(X)
  } else{
    VXt = t(nu2*X)%*%solve(diag(n)+(nu2*X)%*%t(X))
    V = Omega - VXt%*%(nu2*X)
  }
  V = 0.5*(V+t(V))

  # sample the truncated normal component
  sampleTruncNorm = matrix(rtruncnorm(n*nSample, a = 0, b = Inf, mean = muTN, sd = sqrt(paramsPFM$sigma2)), nrow = n, ncol = nSample, byrow = F ) # rtruncnorm(10, a = 0, b = Inf, mean = c(-5,5), sd = c(1,1))  to understand the function
  sampleTruncNorm[y==0,] = -sampleTruncNorm[y==0,] # need to adjust the sign of the variables for which y_i is 0

  # linearly trasform the truncated normal samples
  B = VXt%*%sampleTruncNorm

  # sample the multivariate normal component
  L = t(chol(V))
  sampleMultNorm = L%*%matrix(rnorm(nSample*p),p,nSample)

  # obtain a sample from the approximate posterior distribution by combining the multivariate normal a truncated normal components
  betaSUN_PFM = B + sampleMultNorm
}
```

### `getParamsMF`

This function implements the **CAVI** to obtain the optimal MF approximating density. See **Algorithm 1** in the paper.

**Input**:

-   `X`: *n* × *p* matrix of explanatory variables
-   `y`: binary vector of response variables
-   `nu2`: prior variance for *β*<sub>*j*</sub>’s coefficients (*ν*<sub>*p*</sub><sup>2</sup> in the paper)
-   `tolerance`: absolute change in the ELBO\[*q*<sub>MF</sub>(**β**, **z**)\] used to establish convergence
-   `maxIter`: maximum number of allowed iterations before stopping

**Output**: A list containing

-   `meanBeta`: optimal mean parameter **β**<sup>\*</sup> for the mean-field Gaussian approximation
-   `diagV`: optimal marginal posterior variances for the mean-field Gaussian approximation, i.e. the diagonal elements of the  matrix **V**
-   `nIter`: number of iterations required by the **CAVI**, either because it converged or because the maximum number of iterations `maxIter` was reached

``` r
getParamsMF = function(X,y,nu2, tolerance = 1e-2, maxIter = 1e4){
  ######################################################
  # PRECOMPUTATION
  ######################################################
  # define model dimensions
  n = dim(X)[1]
  p = dim(X)[2]

  # compute V and other useful quantities
  if(p<=n) {
    Omega = diag(rep(nu2,p),p,p)
    invOmega = diag(rep(1/nu2,p),p,p)
    V = solve(t(X)%*%X+invOmega)
    diagV = diag(V)
    VXt = V%*%t(X)
    H = X%*%V%*%t(X) # needed for ELBO
  } else{
    VXt = t(nu2*X)%*%solve(diag(1, nrow = n, ncol = n)+(nu2*X)%*%t(X))
    # V = Omega - VXt%*%(nu2*X)
    diagV = nu2*(1-colSums(t(VXt) * X))
    XXt = X%*%t(X)
    H = nu2*XXt%*%solve(diag(1,nrow=n,ncol=n)+nu2*XXt) # needed for ELBO
  }

  # other useful quantites
  XVVXt = t(VXt)%*%VXt # needed for ELBO
  signH = H
  signH[y==0,] = -H[y==0,]

  # initialization of variables
  meanZ = matrix(0,n,1)
  diff = 1
  elbo = -Inf
  nIter=0

  ######################################################
  # CAVI ALGORITHM
  ######################################################

  while(diff > tolerance & nIter < maxIter) {
    elboOld = elbo

    # update parameters
    mu = H%*%meanZ
    meanZ = mu +(2*y-1)*exp(dnorm(mu, log = T) - pnorm((2*y-1)*mu, log.p = T))

    # compute ELBO
    elbo = (t(meanZ)%*%XVVXt%*%meanZ)/nu2 + sum(pnorm(signH%*%meanZ, log.p = T))

    # compute change in ELBO
    diff = abs(elbo-elboOld)

    nIter = nIter+1
    if(nIter %% 100 == 0) {print(paste0("iter: ", nIter, ", ELBO: ", elbo))}
  }
  meanBeta = VXt%*%meanZ

  return(list(meanBeta = meanBeta, diagV = diagV, nIter = nIter))
}
```

### `rSUNpost`

This function **samples from the exact unified skew-normal posterior distribution** *p*(**β** ∣ **y**), obtained by [Durante
(2019)](https://doi.org/10.1093/biomet/asz034), under multivariate normal prior for **β** having the form *p*(**β**) = *ϕ*(**β**; *ν*<sub>*p*</sub><sup>2</sup>**I**<sub>*p*</sub>).

**Input**:

-   `X`: *n* × *p* matrix of explanatory variables
-   `y`: binary vector of response variables
-   `nu2`: prior variance for *β*<sub>*j*</sub>’s coefficients (*ν*<sub>*p*</sub><sup>2</sup> in the paper)
-   `nSample`: number of i.i.d. samples from *p*(**β** ∣ **y**) to generate

**Output**: A *p* × `nSample` matrix, where each column is a sample from *p*(**β** ∣ **y**)

``` r
rSUNpost = function(X,y,nu2,nSample) {
  # define model dimensions
  n = dim(X)[1]
  p = dim(X)[2]

  # get parameters useful for sampling
  Omega = diag(rep(nu2,p),p,p)
  invOmega = diag(rep(1/nu2,p),p,p)

  signX = X
  signX[y==0,] = -X[y==0,]
  Gamma_post_unnormalized = diag(1,n,n)+(nu2*signX)%*%t(signX)
  inv_Gamma_post_unnormalized = solve(Gamma_post_unnormalized)
  s = diag(sqrt(Gamma_post_unnormalized[cbind(1:n,1:n)]),n,n)
  s_1 = diag(1/s[cbind(1:n,1:n)],n,n)
  gamma_post = matrix(0,n,1) # because prior mean is set to 0
  Gamma_post = s_1%*%Gamma_post_unnormalized%*%s_1

  V = Omega-t(nu2*signX)%*%inv_Gamma_post_unnormalized%*%(signX*nu2)
  V = 0.5*(V+t(V))
  L = t(chol(V))

  # compute multiplicative coefficients for the truncated multivariate normal component
  coefTruncNorm = t(nu2*signX)%*%inv_Gamma_post_unnormalized%*%s

  # sample the multivariate normal component
  sampleMultNorm = matrix(rnorm(nSample*p),p,nSample)

  sampleTruncNorm = t(rtmvnorm(n = nSample, mu = rep(0,n), sigma = Gamma_post, lb = -gamma_post, u = rep(Inf,n)))

  # combine the multivariate normal and truncated normal components
  sampleSUN = L%*%sampleMultNorm+coefTruncNorm%*%sampleTruncNorm
}
```
