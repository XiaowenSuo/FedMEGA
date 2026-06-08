#####Binary  sample number 15000 SNP 200000
library(pryr)
initial_mem <- mem_used()
t_begin_suoglmmnull = proc.time()
library(data.table)
library(SAIGE)
library(Rcpp)
library(psych)
phenoFile='inputM15w20w.txt'
dataM= data.table:::fread(phenoFile, header = T, stringsAsFactors = FALSE, data.table=F)
phenoCol='pheno'
covarColList=c('Gender','PC1','PC2','PC3','PC4','PC5')
M_all<-list()
M_all$X<-dataM[,covarColList]
M_all$y<-dataM[,phenoCol]
head(M_all$X)
phenoFile='inputN15w20w.txt'
dataN= data.table:::fread(phenoFile, header = T, stringsAsFactors = FALSE, data.table=F)
phenoCol='pheno'
covarColList=c('Gender','PC1','PC2','PC3','PC4','PC5')
N_all<-list()
N_all$X<-dataN[,covarColList]
N_all$y<-dataN[,phenoCol]
phenoFile='inputS15w20w.txt'
dataS = data.table:::fread(phenoFile, header = T, stringsAsFactors = FALSE, data.table=F)
phenoCol='pheno'
covarColList=c('Gender','PC1','PC2','PC3','PC4','PC5')
S_all<-list()
S_all$X<-dataS[,covarColList]
S_all$y<-dataS[,phenoCol]

###########################step0
library(pryr)
library(data.table)

initial_mem <- mem_used()
t_begin_suoglmmnull = proc.time()


read_pheno <- function(file){
  dt = fread(file, header = T, stringsAsFactors = FALSE, data.table=F)
  covs = c('Gender','PC1','PC2','PC3','PC4','PC5')
  X = cbind(1, dt[,covs]) 
  y = as.vector(dt[,"pheno"])
  list(X=as.matrix(X), y=y)
}

M_all = read_pheno('inputM15w20w.txt')
N_all = read_pheno('inputN15w20w.txt')
S_all = read_pheno('inputS15w20w.txt')

XN = N_all$X
yN = N_all$y
XM = M_all$X
yM = M_all$y
XS = S_all$X
yS = S_all$y


logistic <- function(x) {
  plogis(x)
}

compute_gradients_and_hessians <- function(X_local, y_local, alpha) {
  alpha <- as.numeric(alpha)

  linpred <- as.vector(X_local %*% alpha)
  p <- pmin(pmax(plogis(linpred), 1e-8), 1 - 1e-8)

  grad <- as.numeric(crossprod(X_local, p - y_local))

  W <- p * (1 - p)
  H <- crossprod(X_local, X_local * W)

  list(grad = grad, Hessian = H)
}

update_alpha <- function(alpha, grad, H, lam = 1e-6) {
  alpha <- as.numeric(alpha)
  grad <- as.numeric(grad)

  np <- length(alpha)
  H_reg <- H + diag(lam, np)

  step <- tryCatch(
    solve(H_reg, grad),
    error = function(e) qr.solve(H_reg, grad)
  )

  as.numeric(alpha - step)
}

distributed_logistic_regression <- function(XN, yN, XM, yM, XS, yS,
                                            alpha_init,
                                            threshold = 1e-6,
                                            max_iter = 500) {
  alpha <- as.numeric(alpha_init)

  for (iter in seq_len(max_iter)) {
    N_local <- compute_gradients_and_hessians(XN, yN, alpha)
    M_local <- compute_gradients_and_hessians(XM, yM, alpha)
    S_local <- compute_gradients_and_hessians(XS, yS, alpha)

    global_grad <- N_local$grad + M_local$grad + S_local$grad
    global_Hessian <- N_local$Hessian + M_local$Hessian + S_local$Hessian

    alpha_new <- update_alpha(alpha, global_grad, global_Hessian)

    if (max(abs(alpha_new - alpha)) < threshold) {
      alpha <- alpha_new
      break
    }

    alpha <- alpha_new

    if (iter == max_iter) {
      warning("Maximum number of iterations reached without convergence.")
    }
  }

  alpha
}


m<-ncol(XN)
alpha_init <- rep(0, m)
final_alpha <- distributed_logistic_regression(XN, yN,XM, yM,XS, yS, alpha_init, threshold = 1e-6)
print(final_alpha)

cost_time = proc.time()-t_begin_suoglmmnull
final_mem = mem_used()


#########################step1

mu_eta<- function(eta){
  p <-exp(eta)/(1+exp(eta))
  return(p*(1-p))
}


sourceCpp(file.path('code', 'cpp', 'inv_sig.cpp'))
sourceCpp(file.path('code', 'cpp', 'server.cpp'))
sourceCpp(file.path('code', 'cpp', 'suofitglmm.cpp'))
sourceCpp(file.path('code', 'cpp', 'wensigmaiX.cpp'))
sourceCpp(file.path('code', 'cpp', 'wenxsigmaix.cpp'))
sourceCpp(file.path('code', 'cpp', 'wensigmaiy.cpp'))
sourceCpp(file.path('code', 'cpp', 'G.cpp'))
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


#Enter the sample numbers for each site
N_n=length(N_all$y)
N_m=length(M_all$y)
N_s=length(S_all$y)

#XN<-cbind(1,N_all$X)
XN=as.matrix(XN)
#XM<-cbind(1,M_all$X)
XM=as.matrix(XM)
#XS<-cbind(1,S_all$X)
XS=as.matrix(XS)
yN=N_all$y
yM=M_all$y
yS=S_all$y


alpha0=final_alpha
etaN <- drop(XN %*% alpha0)
etaM <- drop(XM %*% alpha0)
etaS <- drop(XS %*% alpha0)
eta0 <- c(etaN,etaM,etaS)

fit0=eta0
offset=NULL
family='binomial'
verbose=TRUE
isCovariateOffset=FALSE
maxiter=20
tol = 0.02
tau=c(0,0)
fixtau = c(0,0)
tauInit = c(0,0)
result=suoglmmkin.ai_PCG_Rcpp_Binary_no(N_n,N_m,N_s,yN,yM,yS,XN,XM,XS,fit0, alpha0,isCovariateOffset,tau=c(0,0), fixtau = c(0,0), tauInit = c(0,0),offset,family, maxiter, tol, verbose)
t_end_suoglmmnull =  proc.time()
cat("null model took\n")
print(t_end_suoglmmnull - t_begin_suoglmmnull)
gc()

##########################step 2

t_begin_scoretest = proc.time()
start_mem <- mem_used()
library(pryr)
library(BEDMatrix)
library(data.table)
library(missMethods)
library(stats)
library(Rcpp)

genoM<-BEDMatrix('sampleM15w20w.bed')
genoN<-BEDMatrix('sampleN15w20w.bed')
genoS<-BEDMatrix('sampleS15w20w.bed')
new_colnames <- sub("([0-9]+)_.*", "\\1", colnames(genoM))
colnames(genoM) <- new_colnames
new_colnames <- sub("([0-9]+)_.*", "\\1", colnames(genoN))
colnames(genoN) <- new_colnames
new_colnames <- sub("([0-9]+)_.*", "\\1", colnames(genoS))
colnames(genoS) <- new_colnames
genoS<-genoS[,]
genoM<-genoM[,]
genoN<-genoN[,]

##score=G*res   
res<-result$residuals
resN<-res[1:N_n]
resN<-as.matrix(resN,ncol=1)
resM<-res[(N_n+1):(N_m+N_n)]
resM<-as.matrix(resM,ncol=1)
resS<-res[(N_m+N_n+1):(N_m+N_n+N_s)]
resS<-as.matrix(resS,ncol=1)
dim(resS)
sourceCpp(file.path('code', 'cpp', 'suovar.cpp'))  
SCORE=suoscore(genoN,genoM,genoS,resN,resM,resS)
head(SCORE)

sourceCpp(file.path('code', 'cpp', 'FEDGG.cpp')) 
Z<-result$obj.noK$Z$XXVX_invXV
V<-result$obj.noK$V


XXVX_invXVN <- Z[, 1:N_n, drop = FALSE]
XXVX_invXVM <- Z[, (N_n + 1):(N_n + N_m), drop = FALSE]
XXVX_invXVS <- Z[, (N_n + N_m + 1):(N_n + N_m + N_s), drop = FALSE]

P_null <- result$wen$suo$P
W_null <- as.vector(V)

compute_g_tilde <- function(gN, gM, gS,
                            XXVX_invXVN, XXVX_invXVM, XXVX_invXVS,
                            N_n, N_m, N_s) {
  gN <- as.matrix(gN)
  gM <- as.matrix(gM)
  gS <- as.matrix(gS)

  ZgN <- XXVX_invXVN %*% gN
  ZgM <- XXVX_invXVM %*% gM
  ZgS <- XXVX_invXVS %*% gS
  Zg <- ZgN + ZgM + ZgS

  Zg_N <- Zg[1:N_n, , drop = FALSE]
  Zg_M <- Zg[(N_n + 1):(N_n + N_m), , drop = FALSE]
  Zg_S <- Zg[(N_n + N_m + 1):(N_n + N_m + N_s), , drop = FALSE]

  g_tilde_N <- gN - Zg_N
  g_tilde_M <- gM - Zg_M
  g_tilde_S <- gS - Zg_S

  rbind(g_tilde_N, g_tilde_M, g_tilde_S)
}

estimate_grammar_gamma_ratio_binary <- function(genoN, genoM, genoS,
                                                P_null, W_null,
                                                XXVX_invXVN, XXVX_invXVM, XXVX_invXVS,
                                                N_n, N_m, N_s,
                                                n_ratio_snp = 30,
                                                seed = 2025) {
  n_snp <- ncol(genoM)
  n_ratio_snp <- min(n_ratio_snp, n_snp)

  set.seed(seed)
  ratio_snp <- sample(seq_len(n_snp), n_ratio_snp, replace = FALSE)

  ratio_values <- numeric(n_ratio_snp)

  for (k in seq_along(ratio_snp)) {
    j <- ratio_snp[k]

    g_tilde <- compute_g_tilde(
      gN = genoN[, j, drop = FALSE],
      gM = genoM[, j, drop = FALSE],
      gS = genoS[, j, drop = FALSE],
      XXVX_invXVN = XXVX_invXVN,
      XXVX_invXVM = XXVX_invXVM,
      XXVX_invXVS = XXVX_invXVS,
      N_n = N_n,
      N_m = N_m,
      N_s = N_s
    )

    exact_var <- as.numeric(t(g_tilde) %*% P_null %*% g_tilde)
    approx_var <- as.numeric(t(g_tilde) %*% (W_null * g_tilde))

    ratio_values[k] <- exact_var / approx_var
  }

  ratio_values <- ratio_values[
    is.finite(ratio_values) &
      ratio_values > 0
  ]

  if (length(ratio_values) == 0) {
    stop("No valid GRAMMAR-Gamma variance ratios were obtained.")
  }

  mean(ratio_values)
}

varRatio_null <- estimate_grammar_gamma_ratio_binary(
  genoN = genoN,
  genoM = genoM,
  genoS = genoS,
  P_null = P_null,
  W_null = W_null,
  XXVX_invXVN = XXVX_invXVN,
  XXVX_invXVM = XXVX_invXVM,
  XXVX_invXVS = XXVX_invXVS,
  N_n = N_n,
  N_m = N_m,
  N_s = N_s,
  n_ratio_snp = 30,
  seed = 2025
)

cat("Estimated binary GRAMMAR-Gamma variance ratio:", varRatio_null, "\n")

S <- ncol(genoM)

var <- parallel_fedgg(
  S = S,
  genoM2 = genoM,
  genoN2 = genoN,
  genoS2 = genoS,
  XXVX_invXVM = XXVX_invXVM,
  XXVX_invXVN = XXVX_invXVN,
  XXVX_invXVS = XXVX_invXVS,
  V = V,
  varRatio_null = varRatio_null,
  N_m = N_m,
  N_n = N_n,
  N_s = N_s
)

head(var)







snp<-fread('snp200000.txt',header=FALSE) ##SNP name document
#snp<-snp$x
final <- data.frame(
  SNP=snp,
  T = SCORE,  
  var =var
)
final$pval = pchisq((final$T)^2/final$var, lower.tail = FALSE, df=1)
head(final)
t_end_scoretest =  proc.time()
cat("score test took\n")
print(t_end_scoretest - t_begin_scoretest)
gc()
end_mem <- mem_used()
memory_used <- end_mem - start_mem
memory_used_GB <- memory_used/ (1024^3)
print(paste("score test Memory used:", memory_used_GB, "GB"))
memory_used2 <- end_mem - initial_mem
memory_used_GB <- memory_used2 / (1024^3)
print(paste("Total Memory used:", memory_used_GB, "GB"))

write.table(final, file = "resultbinary15w20w.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
oksnp<-subset(final,final$pval<5e-8)
dim(oksnp)
print(oksnp)


