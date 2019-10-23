Introduction
============

As described in the <tt>README.md</tt> file, this tutorial explains
**how to use the functions implementing the methods presented in Section
2 in the paper**, reporting also the **code in full detail**. Applying
such functions to the Alzheimer dataset contained in the file
<tt>Alzheimer.RData</tt> gives the results presented in Section 3: this
is explained in detail in the tutorial <tt>Alzheimer.html</tt>. If you
want to source directly all the functions, you can download and source
the file <tt>functionsVariational.R</tt>, which contains all the
<tt>R</tt> code reported in this tutorial in one unique file.

Implemented functions
---------------------

The **summary list of implemented functions** is reported below. Each of
them is then analyzed in detail in the following.

-   <tt>getParamsPFM</tt>: returns the **parameters of the optimal PFM
    approximating density** (Algorithm 2);
-   <tt>predictPFM</tt>: returns the **predictive probability of success
    for a new vector of explanatory variables**, according to the
    **optimal PFM approximation** (Proposition 2);
-   <tt>sampleSUN\_PFM</tt>: **samples from the optimal PFM
    approximating density** (Algorithm 3);
-   <tt>getParamsMF</tt>: returns the **parameters of the optimal MF
    approximating density** (Algorithm 1), [Consonni and
    Marin (2007)](https://www.sciencedirect.com/science/article/pii/S0167947306003951)
    ;
-   <tt>predictMF</tt>: returns the **predictive probability of success
    for a new vector of explanatory variables**, according to the
    **optimal MF approximation**;
-   <tt>rSUNpost</tt>: **samples from the exact posterior distribution**
    ([Durante, 2019](https://arxiv.org/abs/1802.09565));
-   <tt>predictMC</tt>: **returns the predictive probability of success
    for a new vector of explanatory variables**, computed via
    **Monte-Carlo integration**.

### <tt>getParamsPFM</tt>

This function implements **Algorithm 2** computing the parameters
**μ**<sup>\*</sup> and **σ**<sup> \* 2</sup> of the optimal PFM
approximating distribution of **z**,
*q*<sub>PFM</sub><sup>\*</sup>(**z**) = ∏<sub>*i* = 1, …, *n*</sub> *q*<sub>PFM</sub><sup>\*</sup>(*z*<sub>*i*</sub>),
with
*q*<sub>PFM</sub><sup>\*</sup>(*z*<sub>*i*</sub>) = *ϕ*(*z*<sub>*i*</sub> − *μ*<sub>*i*</sub><sup>\*</sup>; *σ*<sub>*i*</sub><sup> \* 2</sup>)**1**\[(2*y*<sub>*i*</sub> − 1)*z*<sub>*i*</sub> &gt; 0\]
reported in Theorem 2. Once they are available, the optimal
approximating joint posterior distribution for (**β**, **z**) is simply
*q*<sub>PFM</sub><sup>\*</sup>(**β**, **z**) = *p*(**β** ∣ **z**) ⋅ *q*<sub>PFM</sub><sup>\*</sup>(**z**).

**Input**:

-   <tt>X</tt>: *n* × *p* matrix of explanatory variables;
-   <tt>y</tt>: binary vector of response variables;
-   <tt>nu2</tt>: prior variance for *β*<sub>*i*</sub>’s coefficients
    (*ν*<sup>2</sup> in the paper);
-   <tt>moments</tt>: logical, do you want to obtain marginal moments in
    the output?
-   <tt>tolerance</tt>: absolute change in the
    ELBO\[*q*<sub>PFM</sub>(**z**)\] used to establish convergence;
-   <tt>maxIter</tt>: maximum number of allowed iterations before
    stopping.

**Output**:

List containing:

-   <tt>mu</tt>: optimal value **μ**<sup>\*</sup>;
-   <tt>sigma2</tt>: optimal value **σ**<sup> \* 2</sup>;
-   <tt>nIter</tt>: number of iteration before the algorithm stopped,
    either because it converged or because the maximum number of
    iterations <tt>maxIter</tt> was reached;
-   (optional, if <tt>moments</tt> is set to TRUE) <tt>postMoments</tt>:
    list containing the posterior mean (<tt>postMoments.meanBeta</tt>)
    and marginal posterior variances (<tt>postMoments.varBeta</tt>) of
    **β**.

<!-- -->

    getParamsPFM = function(X,y,nu2,moments = TRUE,tolerance = 1e-2, maxIter = 1e4) { # stopping is determined by change in the ELBO
      ######################################################
      # PRECOMPUTATION
      ######################################################
      # define model dimensions
      n = dim(X)[1]
      p = dim(X)[2]
      
      # define prior covariance matrix and its inverse
      Omega = diag(rep(nu2,p),p,p)
      invOmega = diag(rep(1/nu2,p),p,p)
      
      # compute S = X%*%V%*%t(X) and  directly or with Woodbury
      if(p<=n) {
        V = solve(t(X)%*%X+invOmega)
        S = X%*%V%*%t(X)
        invOmZ = diag(1,nrow=n,ncol=n) - S # needed for ELBO
      } else{
        XXt = X%*%t(X)
        invOmZ = solve(diag(1,nrow=n,ncol=n)+nu2*XXt) # needed for ELBO
        S = nu2*XXt%*%invOmZ
      }

      # compute optimal sigma2
      h = diag(diag(S))
      sigma2 = matrix(1/(1-diag(S)), ncol = 1)
      sigma = sqrt(sigma2)
      
      # compute matrix to write the CAVI update in a vectorized form
      A = diag(as.double(sigma2), nrow = n, ncol = n)%*%(S - h)
      
      # other useful quantities needed for ELBO
      diagInvOmZ = diag(invOmZ)
      coeffMean_Z2 = diagInvOmZ-1/sigma2
      
      # initialization of variables
      meanZ = matrix(0,n,1)
      mean_Z2 = matrix(0,n,1)
      mu = matrix(0,n,1)
      elbo = -1
      diff = 1
      nIter=0
      
      ######################################################
      # CAVI ALGORITHM
      ######################################################
      
      while(diff > tolerance & nIter < maxIter) {
        elboOld = elbo
        
        for(i in 1:n) {
          mu[i] = A[i,]%*%meanZ
          
          # compute first (needed for algorithm) and second (needed for ELBO) moments
          musiRatio = mu[i]/sigma[i]
          phiPhiRatio = dnorm(musiRatio)/pnorm((2*y[i]-1)*musiRatio)
          meanZ[i] = mu[i] + (2*y[i]-1)*sigma[i]*phiPhiRatio
          mean_Z2[i] = mu[i]^2+sigma2[i]+(2*y[i]-1)*mu[i]*sigma[i]*phiPhiRatio # needed for ELBO
        }
        
        # computation of ELBO (up to an additive constant not depending on z)
        elbo = -(t(meanZ)%*%invOmZ%*%meanZ -
                   sum((meanZ^2)*diagInvOmZ) +
                   sum(mean_Z2*coeffMean_Z2))/2 +
              t(meanZ/sigma2)%*%mu
        
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
          # V already computed
          VXt = V%*%t(X)
        } else{ # use Woodbury
          VXt = t(nu2*X)%*%solve(diag(n)+(nu2*X)%*%t(X))
          V = Omega - VXt%*%(nu2*X)
        }
        
        musiRatio = mu/sigma
        phiPhiRatio = dnorm(musiRatio)/pnorm((2*y-1)*musiRatio)
        
        meanZ = mu + (2*y-1)*sigma*phiPhiRatio
        postSdZ = as.double(sqrt(sigma2*(1-(2*y-1)*musiRatio*phiPhiRatio - phiPhiRatio^2)))
        
        C = t(t(VXt)*postSdZ) # same as VXt%*%diag(postSdZ) but much faster
        
        W = apply(C,1,function(x) sum(x*x))
        
        meanBeta = VXt%*%meanZ
        varBeta = diag(V) + W
        
        moments_PFM = list(meanBeta=meanBeta,varBeta=matrix(varBeta,ncol = 1))
        
        results = c(results,postMoments=moments_PFM)
      }
      
      return(results)
    }

### <tt>predictPFM</tt>

This function returns the **approximate predictive probability**
pr<sub>PFM</sub>(*y*<sub>NEW</sub> = 1 ∣ **y**) = ∫*Φ*(**x**<sub>NEW</sub><sup>⊺</sup>**β**)*q*<sub>PFM</sub><sup>\*</sup>(**β**)*d***β** = 𝔼<sub>*q*<sub>PFM</sub><sup>\*</sup>(**z**)</sub>{*Φ*\[**x**<sub>NEW</sub><sup>⊺</sup> **V** **X**<sup>⊺</sup> **z**(1 + **x**<sub>NEW</sub><sup>⊺</sup> **V** **x**<sub>NEW</sub>)<sup> − 1/2</sup>\]},
see **Proposition 2** in the paper.

Such a predictive probability is then computed via Monte-Carlo
integration as

`nSample`<sup> − 1</sup>∑<sub>*r* = 1, …, `nSample`</sub> *Φ*\[**x**<sub>NEW</sub><sup>⊺</sup> **V** **X**<sup>⊺</sup> **z**<sup>(*r*)</sup>(1 + **x**<sub>NEW</sub><sup>⊺</sup> **V** **x**<sub>NEW</sub>)<sup> − 1/2</sup>\],

where **z**<sup>(*r*)</sup> are i.i.d. samples from
*q*<sub>PFM</sub><sup>\*</sup>(**z**), which are straightforward to
sample as one only needs to sample from univariate truncated normals,
thanks to the factorization of *q*<sub>PFM</sub><sup>\*</sup>(**z**).

**Input**:

-   <tt>xNew</tt>: *p* × 1 matrix (column vector) of the new
    observation’s covariates;
-   <tt>paramsPFM</tt>: output of the function <tt>getParamsPFM</tt>;
-   <tt>X</tt>: *n* × *p* matrix of explanatory variables;
-   <tt>y</tt>: binary vector of response variables;
-   <tt>nSample</tt>: number of i.i.d. samples from
    *q*<sub>PFM</sub><sup>\*</sup>(**z**) to be used for the computation
    of the approximate predictive probability.

**Output**: approximate predictive probability
pr<sub>PFM</sub>(*y*<sub>NEW</sub> = 1 ∣ **y**).

**Remark**: in order to properly work, **this function needs the
matrix**
**V** = (*ν*<sup> − 2</sup>**I**<sub>*p*</sub> + **X**<sup>⊺</sup>**X**)<sup> − 1</sup>
**to be already computed** and stored in the global environment. Since
it is a *p* × *p* matrix, recomputing it for each new observation or
even passing it as a function parameter could be very inefficient.
Computation of *V* can be performed efficiently exploiting **Woodbury’s
identity**.

    predictPFM = function(xNew,paramsPFM,X,y,nSample) {
      # we take xNew as a COLUMN vector
      # get number of observations and useful quantities for sampling
      n = dim(X)[1]
      signY = 2*y-1
      muTN = signY*paramsPFM$mu # generate all the truncated as left truncated, with support (0,Inf)
      
      # sample the truncated normals
      sampleTruncNorm = matrix(rtruncnorm(n*nSample, a = 0, b = Inf, mean = muTN, sd = sqrt(paramsPFM$sigma2)), nrow = n, ncol = nSample, byrow = F ) # rtruncnorm(10, a = 0, b = Inf, mean = c(-5,5), sd = c(1,1))  to understand the function
      sampleTruncNorm = apply(sampleTruncNorm,2, function(x) x*signY) # need to adjust the sign of the variables for which y_i is 0
      
      # obtain predictive probability Monte-Carlo estimate
      sd = as.double(sqrt(1+t(xNew)%*%V%*%xNew))
      predProb = mean(pnorm((t(xNew)%*%VXt%*%sampleTruncNorm)/sd))
    }

### <tt>sampleSUN\_PFM</tt>

This function **samples from the approximating optimal
partially-factorized mean-field distribution**
*q*<sub>PFM</sub><sup>\*</sup>(**β**), after the optimal parameters
**μ**<sup>\*</sup> and **σ**<sup> \* 2</sup> have been obtained using
the function <tt>getParamsPFM</tt>. The sampling is performed exploiting
the factorization
*q*<sub>PFM</sub><sup>\*</sup>(**β**, **z**) = *p*(**β** ∣ **z**) ⋅ *q*<sub>PFM</sub><sup>\*</sup>(**z**),
keeping only the samples for **β**. It implements **Algorithm 3**, with
some slight tuning modifications.

**Input**:

-   <tt>paramsPFM</tt>: output of the function <tt>getParamsPFM</tt>;
-   <tt>X</tt>: *n* × *p* matrix of explanatory variables;
-   <tt>y</tt>: binary vector of response variables;
-   <tt>nu2</tt>: prior variance for *β*<sub>*i*</sub>’s coefficients
    (*ν*<sup>2</sup> in the paper);
-   <tt>nSample</tt>: number of i.i.d. samples from
    *q*<sub>PFM</sub><sup>\*</sup>(**β**) to generate.

**Output**: *p* × `nSample` matrix, where each column is a sample from
*q*<sub>PFM</sub><sup>\*</sup>(**β**).

    sampleSUN_PFM = function(paramsPFM, X, y, nu2, nSample = 1e4) {
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
      sampleTruncNormm[y==0,] = -sampleTruncNormm[y==0,] # need to adjust the sign of the variables for which y_i is 0
      
      # linearly trasform the truncated normal samples
      B = VXt%*%sampleTruncNorm
      
      # sample the multivariate normal component
      L = t(chol(V))
      sampleMultNorm = L%*%matrix(rnorm(nSample*p),p,nSample)
      
      # obtain a sample from the approximate posterior distribution by combining the multivariate normal a truncated normal components
      betaSUN_PFM = B + sampleMultNorm
    }

### <tt>getParamsMF</tt>

This function implements **Algorithm 1** computing the parameters
**β**<sup>\*</sup> and **V** of the optimal MF approximating
distribution of **β**, i.e.
*q*<sub>MF</sub><sup>\*</sup>(**β**) = *ϕ*(**β** − **β**<sup>\*</sup>; **V**).

**Input**:

-   <tt>X</tt>: *n* × *p* matrix of explanatory variables;
-   <tt>y</tt>: binary vector of response variables;
-   <tt>nu2</tt>: prior variance for *β*<sub>*i*</sub>’s coefficients
    (*ν*<sup>2</sup> in the paper);
-   <tt>tolerance</tt>: absolute change in the
    ELBO\[*q*<sub>MF</sub>(**β**, **z**)\] used to establish
    convergence;
-   <tt>maxIter</tt>: maximum number of allowed iterations before
    stopping;

**Output**:

List containing:

-   <tt>meanBeta</tt>: optimal mean parameter **β**<sup>\*</sup> for the
    mean-field gaussian approximation;
-   <tt>V</tt>: optimal covariance matrix **V** for the mean-field
    Gaussian approximation;
-   <tt>nIter</tt>: number of iteration before the algorithm stopped,
    either because it converged or because the maximum number of
    iterations <tt>maxIter</tt> was reached;

<!-- -->

    getParamsMF = function(X,y,nu2, tolerance = 1e-2, maxIter = 1e4){
      ######################################################
      # PRECOMPUTATION
      ######################################################
      # define model dimensions
      n = dim(X)[1]
      p = dim(X)[2]
      
      # define prior covariance matrix and its inverse
      Omega = diag(rep(nu2,p),p,p)
      invOmega = diag(rep(1/nu2,p),p,p)
      
      # compute V (optimal covariance matrix for beta) directly or with Woodbury
      if(p<=n) {
        V = solve(t(X)%*%X+invOmega)
        VXt = V%*%t(X)
      } else{
        VXt = t(nu2*X)%*%solve(diag(1, nrow = n, ncol = n)+(nu2*X)%*%t(X))
        V = Omega - VXt%*%(nu2*X)
      }
      V = 0.5*(V+t(V))
      
      # other useful quantities needed for ELBO
      signX = X
      signX[y==0,] = -X[y==0,]
      
      # initialization of variables
      meanZ = matrix(0,n,1)
      meanBeta = matrix(0,p,1)
      diff = 1
      elbo = -1
      nIter=0
      
      ######################################################
      # CAVI ALGORITHM
      ######################################################
      
      while(diff > tolerance & nIter < maxIter) {
        elboOld = elbo
        meanZOld = meanZ
        meanBetaOld = meanBeta
        for(i in 1:n) {
          mu = X[i,]%*%meanBeta
          if(y[i]>0) {
            meanZ[i,] = mu + dnorm(mu)/pnorm(mu)
          } else {
            meanZ[i,] = mu - dnorm(mu)/(1-pnorm(mu))
          }
        }
        
        meanBeta = VXt%*%meanZ
        
        # compute ELBO
        elbo = sum(meanBeta^2)/nu2 + sum(log(pnorm(signX%*%meanBeta)))
        
        # compute change in ELBO
        diff = abs(elbo-elboOld)

        nIter = nIter+1
        if(nIter %% 100 == 0) {print(paste0("iter: ", nIter, ", ELBO: ", elbo))}
      }
      return(list(meanBeta = meanBeta, V = V, nIter = nIter))
    }

### <tt>predictMF</tt>

This function returns the **predictive probability of success for a new
vector of explanatory variables, according to the optimal MF
approximation**, i.e. it outputs
pr<sub>MF</sub>(*y*<sub>NEW</sub> = 1 ∣ **y**) = *Φ*\[**x**<sub>NEW</sub><sup>⊺</sup>(1 + **x**<sub>NEW</sub>**V** **x**<sub>NEW</sub>)<sup> − 1/2</sup>\].

**Input**:

-   <tt>xNew</tt>: *p* × 1 matrix (column vector) of the new
    observation’s covariates;
-   <tt>paramsMF</tt>: output of the function <tt>getParamsMF</tt>.

**Output**: mean-field predictive probability
pr<sub>MF</sub>(*y*<sub>NEW</sub> = 1 ∣ **y**).

    predictMF = function(xNew,paramsMF){
      pnorm(t(xNew)%*%paramsMF$meanBeta/as.double(sqrt(1+t(xNew)%*%paramsMF$V%*%xNew)))
    }

### <tt>rSUNpost</tt>

This function **samples from the true unified skew-normal posterior
distribution** *p*(**β** ∣ **y**), obtained by [Durante
(2019)](https://arxiv.org/abs/1802.09565), under multivariate normal
prior for **β** having the form
*p*(**β**) = *ϕ*(**β**; *ν*<sup>2</sup>**I**<sub>*p*</sub>).

**Input**:

-   <tt>X</tt>: *n* × *p* matrix of explanatory variables;
-   <tt>y</tt>: binary vector of response variables;
-   <tt>nu2</tt>: prior variance for *β*<sub>*i*</sub>’s coefficients
    (*ν*<sup>2</sup> in the paper);
-   <tt>nSample</tt>: number of i.i.d. samples from *p*(**β** ∣ **y**)
    to generate.

**Output**: *p* × `nSample` matrix, where each column is a sample from
*p*(**β** ∣ **y**).

    rSUNpost = function(X,y,nu2,nSample) {
      # define model dimensions
      n = dim(X)[1]
      p = dim(X)[2]
      
      # get parameters useful for sampling
      Omega = diag(rep(nu2,p),p,p)
      invOmega = diag(rep(1/nu2,p),p,p)
      nu = sqrt(nu2)
      
      D = X
      D[y==0,] = -X[y==0,]
      Gamma_post_unnormalized = (nu2*D)%*%t(D)+diag(1,n,n)
      inv_Gamma_post_unnormalized = solve(Gamma_post_unnormalized)
      s = diag(sqrt(Gamma_post_unnormalized[cbind(1:n,1:n)]),n,n)
      s_1 = diag(1/s[cbind(1:n,1:n)],n,n)
      gamma_post = matrix(0,n,1) # because prior mean is set to 0
      Gamma_post = s_1%*%Gamma_post_unnormalized%*%s_1
      
      Var_V0 = diag(1,p,p)-t(nu*D)%*%inv_Gamma_post_unnormalized%*%(D*nu)
      Var_V0 = 0.5*(Var_V0+t(Var_V0))
      L = t(chol(Var_V0))
      
      # compute multiplicative coefficients for the multivariate normal component
      coefMultNorm = nu*L
      
      # compute multiplicative coefficients for the truncated multivariate normal component
      coefTruncNorm = t(nu2*D)%*%inv_Gamma_post_unnormalized%*%s
      
      # sample the multivariate normal component
      sampleMultNorm = matrix(rnorm(nSample*p),p,nSample)
      
      # sample the truncated multivariate normal component
      if(n == 1) {
        sampleTruncNorm = matrix(rtruncnorm(n = Nsample, a = -gamma_post, b = Inf, mean = 0, sd = 1), nrow = 1, ncol = nSample)
      } else{
        sampleTruncNorm = mvrandn(l = -gamma_post, u = rep(Inf,n), Sig = Gamma_post, n = nSample)
      }
      
      # combine the multivariate normal and truncated normal components
      sampleSUN = coefMultNorm%*%sampleMultNorm+coefTruncNorm%*%sampleTruncNorm
    }

### <tt>predictMC</tt>

This function returns the **predictive probability of success for a new
vector of explanatory variables** **x**<sub>NEW</sub>.
pr(*y*<sub>NEW</sub> = 1 ∣ **y**) = ∫*Φ*(**x**<sub>NEW</sub><sup>⊺</sup>**β**)*p*(**β** ∣ **y**)*d***β**,
computed via Monte-Carlo integration as

`nSample`<sup> − 1</sup>∑<sub>*r* = 1, …, `nSample`</sub> *Φ*(**x**<sub>NEW</sub><sup>⊺</sup>**β**<sup>(*r*)</sup>)

using the i.i.d. samples
**β**<sup>(1)</sup>, …, **β**<sup>`(nSample)`</sup> from the posterior
distribution *p*(**β** ∣ **y**) obtained with the function
<tt>rSUNpost</tt>.

**Input**:

-   <tt>xNew</tt>: *p* × 1 matrix (column vector) of the new
    observation’s covariates;
-   <tt>betaSUN</tt>: output of the function <tt>rSUNpost</tt>.

**Output**: predictive probability pr(*y*<sub>NEW</sub> = 1 ∣ **y**).

    predictMC = function(xNew,betaSUN) {
      PredProb = mean(pnorm(t(xNew)%*%betaSUN))
    }
