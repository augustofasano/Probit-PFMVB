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

predictMF = function(xNew,paramsMF){
  pnorm(t(xNew)%*%paramsMF$meanBeta/as.double(sqrt(1+t(xNew)%*%paramsMF$V%*%xNew)))
}

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

predictMC = function(xNew,betaSUN) {
  PredProb = mean(pnorm(t(xNew)%*%betaSUN))
}