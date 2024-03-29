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
    V = solve(crossprod(X)+invOmega)
    H = X%*%V%*%t(X)
    invOmZ = diag(1,nrow=n,ncol=n) - H # needed for ELBO
  } else{
    XXt = tcrossprod(X)
    invOmZ = solve(diag(1,nrow=n,ncol=n)+nu2*XXt) # needed for ELBO
    H = nu2*XXt%*%invOmZ
  }
  
  # compute optimal sigma2
  # h = diag(diag(H))
  sigma2 = as.double(1/(1-diag(H)), ncol = 1)
  sigma = sqrt(sigma2)
  
  # compute matrix to write the CAVI update in a vectorized form
  A = sigma2*H
  A[cbind(1:n,1:n)] = 0
  
  # other useful quantities needed for ELBO
  diagInvOmZ = diag(invOmZ)
  
  # initialization of variables
  mu = matrix(0,n,1)
  
  # inizialize coherently the vector of expectations meanZ
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
      sumLogPhi = sumLogPhi + pnorm((2*y[i]-1)*musiRatio, log.p = T)
    }
    
    # computation of ELBO (up to an additive constant not depending on mu)
    elbo = -(t(meanZ)%*%invOmZ%*%meanZ -
               sum((meanZ^2)*diagInvOmZ))/2 -
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
      VXt = t(nu2*X)%*%invOmZ
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
  
  # compute V and other useful quantities
  if(p<=n) {
    Omega = diag(rep(nu2,p),p,p)
    invOmega = diag(rep(1/nu2,p),p,p)
    V = solve(crossprod(X)+invOmega)
    diagV = diag(V)
    VXt = V%*%t(X)
    H = X%*%VXt # needed for ELBO
  } else{
    XXt = tcrossprod(X)
    invOmZ = solve(diag(1,nrow=n,ncol=n)+nu2*XXt)
    VXt = t(nu2*X)%*%invOmZ
    diagV = nu2*(1-colSums(t(VXt) * X))
    H = nu2*XXt%*%invOmZ # needed for ELBO
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
    elbo = -0.5*(t(meanZ)%*%XVVXt%*%meanZ)/nu2 + sum(pnorm(signH%*%meanZ, log.p = T))
    
    # compute change in ELBO
    diff = abs(elbo-elboOld)
    
    nIter = nIter+1
    if(nIter %% 100 == 0) {print(paste0("iter: ", nIter, ", ELBO: ", elbo))}
  }
  meanBeta = VXt%*%meanZ
  
  return(list(meanBeta = meanBeta, diagV = diagV, nIter = nIter))
}

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
  
  # sample the truncated multivariate normal component
  sampleTruncNorm = t(rtmvnorm(n = nSample, mu = rep(0,n), sigma = Gamma_post, lb = -gamma_post, u = rep(Inf,n)))
  
  # combine the multivariate normal and truncated normal components
  sampleSUN = L%*%sampleMultNorm+coefTruncNorm%*%sampleTruncNorm
}

sampleSVB = function(paramsSVB,nSample){
  # samples from the sparse-vb spike-and-slab approximation
  p = length(paramsSVB$mu)
  betaSample = matrix(0, nrow = p, ncol = nSample)
  for (i in 1:nSample) {
    inclusionIndicator = which(runif(n = p) <= paramsSVB$gamma)
    if(length(inclusionIndicator)>0){
      betaSample[inclusionIndicator,i] = rnorm(n = length(inclusionIndicator),
                                               mean = paramsSVB$mu[inclusionIndicator],
                                               sd = paramsSVB$sigma[inclusionIndicator])
    }
  }
  betaSample
}