suoget_coef_no_q = function(yN,yM,yS,XN,XM,XS, N_n,N_m,N_s,tau, alpha0, eta0, offset, maxiter,verbose){
  
    #Client-side computation protect y
    tol.coef = 0.1
    etaN=eta0[1:N_n]
    etaM=eta0[(N_n+1):(N_n+N_m)]
    etaS=eta0[(N_n+N_m+1):(N_n+N_m+N_s)]
    offsetN=offset[1:N_n]
    offsetM=offset[(N_n+1):(N_n+N_m)]
    offsetS=offset[(N_n+N_m+1):(N_n+N_m+N_s)]
    muN = etaN
    muM = etaM
    muS = etaS
    mu.etaN = rep(1,N_n)
    mu.etaM = rep(1,N_m)
    mu.etaS = rep(1,N_s)
    YN = etaN - offsetN + (yN - muN)/mu.etaN
    YM = etaM - offsetM + (yM - muM)/mu.etaM
    YS = etaS - offsetS + (yS - muS)/mu.etaS
    Y=c(YN,YM,YS)
    mu.eta=c(mu.etaN,mu.etaM,mu.etaS)
    mu=c(muN,muM,muS)
    var_mu=rep(1,(N_n+N_m+N_s))
    sqrtW = mu.eta / sqrt(var_mu)
    W = sqrtW^2
    gc()
    
 for(i in 1:maxiter){
    invsigma<-invsig(W,tau) #client1
    #invsigma<-invsigma$invsig
    #invsig client 1 sends to other clients
    invsigmaN<-invsigma[(1:N_n),(1:N_n)]  
    invsigmaM<-invsigma[(N_n+1):(N_n+N_m),(N_n+1):(N_n+N_m)]
    invsigmaS<-invsigma[(N_n+N_m+1):(N_n+N_m+N_s),(N_n+N_m+1):(N_n+N_m+N_s)]
    
    #Client-side computation
    sigmaiXM<-wensigmaiX(invsigmaM, XM)  
    sigmaiXN<-wensigmaiX(invsigmaN, XN)
    sigmaiXS<-wensigmaiX(invsigmaS, XS)
    XsigmaiXM<-getxsigmaix(XM, sigmaiXM)
    XsigmaiXN<-getxsigmaix(XN, sigmaiXN)
    XsigmaiXS<-getxsigmaix(XS, sigmaiXS)
    
    XsigmaiX<-XsigmaiXN+XsigmaiXM+XsigmaiXS  #clients to server
    sigmaiX<-rbind(sigmaiXN,sigmaiXM) #clients to server
    sigmaiX<-rbind(sigmaiX,sigmaiXS)
    
    #client 1 computation,send to server,then server computes.
    sigmaiY<-wensigmaiY(invsigma, Y)  
    re.coef=server_q(invsigma,sigmaiX, XsigmaiX,sigmaiY,Y,W,tau)
    alpha = re.coef$alpha
    eta = re.coef$eta + offset
    suo=re.coef
    gc()
    
   if(verbose) {
      cat("Tau:\n")
      print(tau)
      cat("Fixed-effect coefficients:\n")
      print(alpha)
   }
   #Client-side computation
    etaN=eta0[1:N_n]
    etaM=eta0[(N_n+1):(N_n+N_m)]
    etaS=eta0[(N_n+N_m+1):(N_n+N_m+N_s)]
    offsetN=offset[1:N_n]
    offsetM=offset[(N_n+1):(N_n+N_m)]
    offsetS=offset[(N_n+N_m+1):(N_n+N_m+N_s)]
    muN = etaN
    muM = etaM
    muS = etaS
    mu.etaN = rep(1,N_n)
    mu.etaM = rep(1,N_m)
    mu.etaS = rep(1,N_s)
   YN = etaN - offsetN + (yN - muN)/mu.etaN
   YM = etaM - offsetM + (yM - muM)/mu.etaM
   YS = etaS - offsetS + (yS - muS)/mu.etaS
   Y=c(YN,YM,YS)
   mu.eta=c(mu.etaN,mu.etaM,mu.etaS)
   #Y = eta - offset + (y - mu)/mu.eta
   var_mu=rep(1,(N_n+N_m+N_s))
   sqrtW = mu.eta / sqrt(var_mu)
   W = sqrtW^2
   gc()
   
   if( max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol.coef))< tol.coef){
     break
   }
     alpha0 = alpha
   }
   
   re = list(suo=suo, Y=Y, alpha=alpha, eta=eta, W=W, cov=re.coef$cov, Sigma_iY = re.coef$sigmaiY, Sigma_iX = re.coef$sigmaiX, mu=mu)
}



suoglmmkin.ai_PCG_Rcpp_Quan_no = function(N_n,N_m,N_s,yN,yM,yS,XN,XM,XS,eta0,alpha0, isCovariateOffset,tau=c(0,0), fixtau = c(0,0), tauInit = c(0,0),offset,family, maxiter, tol, verbose) {
 
 t_begin = proc.time()
 n=N_n+N_m+N_s
 if(is.null(offset)){
    offset = rep(0, n)
 }
 offsetN=offset[1:N_n]
 offsetM=offset[(N_n+1):(N_n+N_m)]
 offsetS=offset[(N_n+N_m+1):n]
 etaN=eta0[1:N_n]
 etaM=eta0[(N_n+1):(N_n+N_m)]
 etaS=eta0[(N_n+N_m+1):n]
 muN = etaN
 muM = etaM
 muS = etaS
 mu.etaN = rep(1,N_n)
 mu.etaM = rep(1,N_m)
 mu.etaS = rep(1,N_s)
 YN = etaN - offsetN + (yN - muN)/mu.etaN
 YM = etaM - offsetM + (yM - muM)/mu.etaM
 YS = etaS - offsetS + (yS - muS)/mu.etaS
 Y=c(YN,YM,YS)
 mu.eta=c(mu.etaN,mu.etaM,mu.etaS)
 mu=c(muN,muM,muS)
 var_mu=rep(1,(N_n+N_m+N_s))
 sqrtW = mu.eta / sqrt(var_mu)
 W = sqrtW^2
 
 alpha = alpha0 
 alpha0=alpha
 eta0 = eta0
 

if(verbose) {
    cat("Fixed-effect coefficients:\n")
    print(alpha)
  }



if(family %in% c("poisson", "binomial")) {
tau[1] = 1
fixtau[1] = 1
}

q = 1

if(sum(tauInit[fixtau == 0]) == 0){
    #tau[fixtau == 0] = var(Y)/(q+1)
    tau[1] = 1
    tau[2] = 0
    if (abs(var(Y)) < 0.1){
    stop("WARNING: variance of the phenotype is much smaller than 1. Please consider invNormalize=T\n")
    }
}else{
    tau[fixtau == 0] = tauInit[fixtau == 0]
  }

tau0=tau
cat("initial tau is ", tau,"\n")

re.coef<-suoget_coef_no_q(yN,yM,yS,XN,XM,XS, N_n,N_m,N_s,tau,alpha0, eta0, offset, maxiter,verbose)
re<-re.coef$suo
tau[2] = max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$trace[2])/n)
tau[1] = max(0, tau0[1] + tau0[1]^2 * (re$YPA0PY - re$trace[1])/n)
gc()

if (verbose) {
    cat("Variance component estimates:\n")
    print(tau)
}

 for (i in seq_len(maxiter)) {
    
    W = sqrtW^2
   
    if(verbose) 
    cat("\nIteration ", i, tau, ":\n")
    alpha0 = alpha
    tau0 = tau
    cat("tau0_v1: ", tau0, "\n")
    eta0 = eta0
    
    re.coef<-suoget_coef_no_q(yN,yM,yS,XN,XM,XS, N_n,N_m,N_s,tau, alpha0, eta0, offset, maxiter,verbose)
    fit = suofitglmm_q(tau, re.coef$suo$score,re.coef$suo$AI, tol = tol)
    gc()
    
  if(tau[1]!=0){
    tau = as.numeric(fit$tau)
    cov = re.coef$cov
    alpha = re.coef$alpha
    eta = re.coef$eta
    cat("cov: ", cov, "\n")

    if(verbose) {
      cat("Variance component estimates:\n")
      print(tau)
      cat("Fixed-effect coefficients:\n")
      print(alpha)
    }
    Y = re.coef$Y
    mu = re.coef$mu


    if(tau[1]<=0 | tau[2] <= 0) break
    
    if(max(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) break
         
    if(max(tau) > tol^(-2)) {
        warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
      	i = maxiter
      	break
      }
}else{#if(tau[1]!=0){
     break
    }
  }

  if(verbose) cat("\nFinal " ,tau, ":\n")
  
  t_end_null = proc.time()
  cat("t_end_null - t_begin, fitting the NULL model without LOCO took\n")
  print(t_end_null - t_begin)
  gc()
  
  re.coef = suoget_coef_no_q(yN,yM,yS,XN,XM,XS, N_n,N_m,N_s,tau, alpha, eta, offset, maxiter,verbose)
  cov = re.coef$cov
  alpha = re.coef$alpha
  eta = re.coef$eta
  Y = re.coef$Y
  mu = re.coef$mu
  wen=re.coef
  
  muN=mu[1:N_n]
  muM=mu[(N_n+1):(N_n+N_m)]
  muS=mu[(N_n+N_m+1):n]
  resN = yN - muN
  resM = yM - muM
  resS = yS - muS
  res=c(resN,resM,resS)
  mu2 = rep((1/(tau[1])),length(res))
  converged = ifelse(i < maxiter, TRUE, FALSE)
  gc()
  
  if(!isCovariateOffset){
  obj.noK = NULL
  glmmResult = list(wen=wen,theta=tau, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged, obj.noK=obj.noK, traitType="Quan", isCovariateOffset = isCovariateOffset)
  }else{
  glmmResult = list(wen=wen, theta=tau, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged,  traitType="Quan", isCovariateOffset = isCovariateOffset)
}
 return(glmmResult)
}
