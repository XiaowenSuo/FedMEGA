mu_eta<- function(eta){
  p <-exp(eta)/(1+exp(eta))
  return(p*(1-p))
}

suoget_coef_no = function(yN,yM,yS,XN,XM,XS, N_n,N_m,N_s,tau, alpha0, eta0, offset, maxiter,verbose){
    
    #Client-side computation protect y
    tol.coef = 0.1
    etaN=eta0[1:N_n]
    etaM=eta0[(N_n+1):(N_n+N_m)]
    etaS=eta0[(N_n+N_m+1):(N_n+N_m+N_s)]
    muN = plogis(etaN)
    muM = plogis(etaM)
    muS = plogis(etaS)
    mu.etaN = mu_eta(etaN)
    mu.etaM = mu_eta(etaM)
    mu.etaS = mu_eta(etaS)
    offsetN=offset[1:N_n]
    offsetM=offset[(N_n+1):(N_n+N_m)]
    offsetS=offset[(N_n+N_m+1):(N_n+N_m+N_s)]
    YN = etaN - offsetN + (yN - muN)/mu.etaN
    YM = etaM - offsetM + (yM - muM)/mu.etaM
    YS = etaS - offsetS + (yS - muS)/mu.etaS
    
    Y=c(YN,YM,YS)
    mu.eta=c(mu.etaN,mu.etaM,mu.etaS)
    mu=c(muN,muM,muS)
    sqrtW = mu.eta / sqrt((mu * (1 - mu)))
    W = sqrtW^2


 for(i in 1:maxiter){
    invsigma<-invsig(W,tau) #sever 
    #invsigma<-invsigma$invsig
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
    
    #Server-side computation
    sigmaiY<-wensigmaiY(invsigma, Y)  
    re.coef=server(invsigma,sigmaiX, XsigmaiX,sigmaiY,Y,W,tau)
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
   etaN=eta[1:N_n]
   etaM=eta[(N_n+1):(N_n+N_m)]
   etaS=eta[(N_n+N_m+1):(N_n+N_m+N_s)]
   muN = plogis(etaN)
   muM = plogis(etaM)
   muS = plogis(etaS)
   mu=c(muN,muM,muS)
   mu.etaN = mu_eta(etaN)
   mu.etaM = mu_eta(etaM)
   mu.etaS = mu_eta(etaS)
   YN = etaN - offsetN + (yN - muN)/mu.etaN
   YM = etaM - offsetM + (yM - muM)/mu.etaM
   YS = etaS - offsetS + (yS - muS)/mu.etaS
   Y=c(YN,YM,YS)
   mu.eta=c(mu.etaN,mu.etaM,mu.etaS)
   #Y = eta - offset + (y - mu)/mu.eta
   sqrtW = mu.eta / sqrt((mu * (1 - mu)))
   W = sqrtW^2

   if( max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol.coef))< tol.coef){
     break
   }
     alpha0 = alpha
   }

   re = list(suo=suo, Y=Y, alpha=alpha, eta=eta, W=W, cov=re.coef$cov, Sigma_iY = re.coef$sigmaiY, Sigma_iX = re.coef$sigmaiX, mu=mu)
}

Fed_ScoreTest_NULL_Model = function(mu2, XN,XM,XS,N_n,N_m,N_s){
  V = as.vector(mu2)
  #server sends Vj to client j
  VN=V[1:N_n]
  VM=V[(N_n+1):(N_n+N_m)]
  VS=V[(N_n+N_m+1):(N_n+N_m+N_s)]
  #clients compute XVj and send XVj to client j=1
  XVn=t(XN*VN)
  XVm=t(XM*VM)
  XVs=t(XS*VS)
  #client j=1 receive XVj and aggregate XV
  XV =cbind(XVn,XVm,XVs)
  #clients compute XVXj and send XVXj to client j=1
  XVXn=t(XN) %*% (t(XVn))
  XVXm=t(XM) %*% (t(XVm))
  XVXs=t(XS) %*% (t(XVs))
  #client j=1 receive XVXj and aggregate XVX
  XVX = XVXn+XVXm+XVXs 
  #client j=1 compute XVX_inv and send it to other clients
  XVX_inv = solve(XVX)  
  #clients compute XXVX_invj and send XXVX_invj to client j=1
  XXVX_invn= as.matrix(XN) %*% XVX_inv
  XXVX_invm= as.matrix(XM) %*% XVX_inv
  XXVX_invs= as.matrix(XS) %*% XVX_inv
  #client j=1 receive XXVX_invj and aggregate XXVX_inv
  XXVX_inv = rbind(XXVX_invn,XXVX_invm,XXVX_invs)
  #client j=1 compute XXVX_invXV
  XXVX_invXV=G_G(XXVX_inv,XV) 
  Z=XXVX_invXV
  re = list(XV = XV, XVX = XVX, XXVX_inv = XXVX_inv, V = V,Z=Z)
  return(re) 
}

suoglmmkin.ai_PCG_Rcpp_Binary_no = function(N_n,N_m,N_s,yN,yM,yS,XN,XM,XS,fit0, alpha0,isCovariateOffset,tau=c(0,0), fixtau = c(0,0), tauInit = c(0,0),offset,family, maxiter, tol, verbose) {
 
 n=N_n+N_m+N_s
 if(is.null(offset)){
    offset = rep(0, n)
 }
 offsetN=offset[1:N_n]
 offsetM=offset[(N_n+1):(N_n+N_m)]
 offsetS=offset[(N_n+N_m+1):n]
 etaN=fit0[1:N_n]
 etaM=fit0[(N_n+1):(N_n+N_m)]
 etaS=fit0[(N_n+N_m+1):n]
 muN = plogis(etaN)
 muM = plogis(etaM)
 muS = plogis(etaS)
 mu.etaN = mu_eta(etaN)
 mu.etaM = mu_eta(etaM)
 mu.etaS = mu_eta(etaS)
 YN = etaN - offsetN + (yN - muN)/mu.etaN
 YM = etaM - offsetM + (yM - muM)/mu.etaM
 YS = etaS - offsetS + (yS - muS)/mu.etaS
 Y=c(YN,YM,YS)
 mu.eta=c(mu.etaN,mu.etaM,mu.etaS)
 mu=c(muN,muM,muS)
 eta=c(etaN,etaM,etaS)
 eta0=eta
 alpha0=alpha0

 if(family %in% c("poisson", "binomial")) {
  tau[1] = 1
  fixtau[1] = 1
 }

 q = 1

 if(tauInit[fixtau == 0] == 0){
  tau[fixtau == 0] = 0.1
 }else{
  tau[fixtau == 0] = tauInit[fixtau == 0]
 }
 cat("inital tau is ", tau,"\n")
 tau0=tau


 re.coef<-suoget_coef_no(yN,yM,yS,XN,XM,XS, N_n,N_m,N_s,tau, alpha0, eta0, offset, maxiter,verbose)
 re<-re.coef$suo
 
 tau[2] = max(0, tau0[2] + tau0[2]^2 * (re$YPAPY - re$trace)/n)
 gc()
 
if (verbose) {
    cat("Variance component estimates:\n")
    print(tau)
}

 for (i in seq_len(maxiter)) {

    if(verbose) 
    cat("\nIteration ", i, tau, ":\n")
    alpha0 = re.coef$alpha
    tau0 = tau
    eta0=eta
    cat("tau0_v1: ", tau0, "\n")
    t_begin_Get_Coef = proc.time()
    re.coef<-suoget_coef_no(yN,yM,yS,XN,XM,XS, N_n,N_m,N_s,tau, alpha0, eta0, offset, maxiter,verbose)
    t_end_Get_Coef =  proc.time()
    cat("Updating fix effect coeffcieints took\n")
    print(t_end_Get_Coef - t_begin_Get_Coef)
    fit = suofitglmm(tau, re.coef$suo$score,re.coef$suo$AI, tol = tol)
    t_end_fitglmmaiRPCG= proc.time()
    cat("Updating variance component estimate took\n")
    print(t_end_fitglmmaiRPCG - t_end_Get_Coef)
    gc()
    
    tau = as.numeric(fit$tau)
    cov = re.coef$cov
    alpha = re.coef$alpha
    eta = re.coef$eta
    Y = re.coef$Y
    mu = re.coef$mu

    print(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol))
    cat("tau: ", tau, "\n")
    cat("tau0: ", tau0, "\n")


    if(tau[2] == 0){
    break
    } 
    
      
    if(max(abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol){
      break
    }  
      

    if(max(tau) > tol^(-2)) {
        warning("Large variance estimate observed in the iterations, model not converged...", call. = FALSE)
      	i = maxiter
      	break
      }
}

  if(verbose) cat("\nFinal " ,tau, ":\n")
  
  
  re.coef = suoget_coef_no(yN,yM,yS,XN,XM,XS, N_n,N_m,N_s,tau, alpha0, eta0, offset, maxiter,verbose)
  cov = re.coef$cov
  alpha = re.coef$alpha
  eta = re.coef$eta
  Y = re.coef$Y
  mu = re.coef$mu
  
  wen=re.coef
  
  converged = ifelse(i < maxiter, TRUE, FALSE)
  gc()
  
  #client
  muN=mu[1:N_n]
  muM=mu[(N_n+1):(N_n+N_m)]
  muS=mu[(N_n+N_m+1):n]
  resN = yN - muN
  resM = yM - muM
  resS = yS - muS
  #server
  res=c(resN,resM,resS)
  mu2 = mu * (1-mu)

  if(!isCovariateOffset){
  obj.noK = Fed_ScoreTest_NULL_Model(mu2,XN,XM,XS,N_n,N_m,N_s)
  glmmResult = list(wen=wen,theta=tau, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged,  traitType="binary", isCovariateOffset = isCovariateOffset,obj.noK= obj.noK )
  }else{
  glmmResult = list(wen=wen, theta=tau, coefficients=alpha, linear.predictors=eta, fitted.values=mu, Y=Y, residuals=res, cov=cov, converged=converged,  traitType="binary", isCovariateOffset = isCovariateOffset,obj.noK= obj.noK)
}
 return(glmmResult)
}
