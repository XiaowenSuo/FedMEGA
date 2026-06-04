library(pryr)
initial_mem <- mem_used()
t_begin_suoglmmnull = proc.time()
library(data.table)
library(SAIGE)
library(Rcpp)
library(psych)
phenoFile='inputM_q15w20w.txt'
dataM= data.table:::fread(phenoFile, header = T, stringsAsFactors = FALSE, data.table=F)
phenoCol='pheno'
covarColList=c('Gender','PC1','PC2','PC3','PC4','PC5')
M_all<-list()
M_all$X<-dataM[,covarColList]
M_all$y<-dataM[,phenoCol]
head(M_all$X)
phenoFile='inputN_q15w20w.txt'
dataN= data.table:::fread(phenoFile, header = T, stringsAsFactors = FALSE, data.table=F)
phenoCol='pheno'
covarColList=c('Gender','PC1','PC2','PC3','PC4','PC5')
N_all<-list()
N_all$X<-dataN[,covarColList]
N_all$y<-dataN[,phenoCol]
phenoFile='inputS_q15w20w.txt'
dataS = data.table:::fread(phenoFile, header = T, stringsAsFactors = FALSE, data.table=F)
phenoCol='pheno'
covarColList=c('Gender','PC1','PC2','PC3','PC4','PC5')
S_all<-list()
S_all$X<-dataS[,covarColList]
S_all$y<-dataS[,phenoCol]

gradient_linear <- function(X, y, beta) {
  n <- nrow(X)
  grad <- numeric(length(beta))
  for (i in 1:n) {
    grad <- grad + 2 * (X[i,] %*% beta - y[i]) * X[i,]
  }
  return(grad)
}

hessian_linear <- function(X) {
  n <- nrow(X)
  p <- ncol(X)
  H <- matrix(0, p, p)
  for (i in 1:n) {
    H <- H + t(X[i, , drop = FALSE]) %*% X[i, ,drop = FALSE ]
  }
  return(H)
}

compute_gradients_and_hessians_linear <- function(X_local, y_local, beta) {
  grad <- gradient_linear(X_local, y_local, beta)
  H <- hessian_linear(X_local)
  return(list(grad = grad, Hessian = H))
}

update_beta <- function(beta, grad, H) {
  beta_new <- beta - solve(H) %*% grad
  return(beta_new)
}

distributed_linear_regression <- function(XN, yN, XM, yM, XS, yS, beta_init, threshold = 1e-6) {
  beta <- beta_init
  converged <- FALSE
  while (!converged) {
    N_local <- compute_gradients_and_hessians_linear(XN, yN, beta)
    M_local <- compute_gradients_and_hessians_linear(XM, yM, beta)
    S_local <- compute_gradients_and_hessians_linear(XS, yS, beta)
    
    global_grad <- Reduce('+', lapply(list(N_local, M_local, S_local), function(result) result$grad))
    global_Hessian <- Reduce('+', lapply(list(N_local, M_local, S_local), function(result) result$Hessian))
    beta_new <- update_beta(beta, global_grad, global_Hessian)
   
    
    if (max(abs(beta_new - beta)) < threshold) {
      converged <- TRUE
    }
    beta <- beta_new
  }
  return(beta)
}

N_n=length(N_all$y)
N_m=length(M_all$y)
N_s=length(S_all$y)

XN<-cbind(1,N_all$X)
XN=as.matrix(XN)
XM<-cbind(1,M_all$X)
XM=as.matrix(XM)
XS<-cbind(1,S_all$X)
XS=as.matrix(XS)
yN=N_all$y
yM=M_all$y
yS=S_all$y
m<-ncol(XN)
alpha_init <- rep(0, m)
final_alpha <- distributed_linear_regression(XN, yN,XM, yM,XS, yS, alpha_init, threshold = 1e-6)
print(final_alpha)

final_etaN<-XN %*% final_alpha
final_etaM<-XM %*% final_alpha
final_etaS<-XS %*% final_alpha
final_eta<-rbind(final_etaN,final_etaM,final_etaS)



##################

sourceCpp(file.path('code', 'cpp', 'inv_sig.cpp'))
sourceCpp(file.path('code', 'cpp', 'server_q.cpp'))
sourceCpp(file.path('code', 'cpp', 'suofitglmm_q.cpp'))
sourceCpp(file.path('code', 'cpp', 'wensigmaiX.cpp'))
sourceCpp(file.path('code', 'cpp', 'wenxsigmaix.cpp'))
sourceCpp(file.path('code', 'cpp', 'wensigmaiy.cpp'))
sourceCpp(file.path('code', 'cpp', 'G.cpp'))
suoget_coef_no_q = function(yN,yM,yS,XN,XM,XS, N_n,N_m,N_s,tau, alpha0, eta0, offset, maxiter,verbose){
  
    #Client-side computation锛歱rotect y
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



alpha0=final_alpha
eta0=final_eta

offset=NULL
family='poisson'
verbose=TRUE
isCovariateOffset=FALSE
maxiter =20
tol = 0.02
tau=c(0,0)
fixtau = c(0,0)
tauInit = c(0,0)
t_begin_suoglmmnull = proc.time()
result=suoglmmkin.ai_PCG_Rcpp_Quan_no(N_n,N_m,N_s,yN,yM,yS,XN,XM,XS,eta0,alpha0, isCovariateOffset,tau=c(0,0), fixtau = c(0,0), tauInit = c(0,0),offset,family, maxiter, tol, verbose)
t_end_suoglmmnull =  proc.time()
cat("null model took\n")
print(t_end_suoglmmnull - t_begin_suoglmmnull)
gc()


##############2


t_begin_test = proc.time()
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

sigma<-result$theta[1]
#sigma_e2<-result$theta[2]
n=N_n+N_m+N_s
V <- sigma* diag(n)
V_inv <- solve(V)

V_invN<-V_inv[(1:N_n),(1:N_n)]
V_invM<-V_inv[(N_n+1):(N_m+N_n),(N_n+1):(N_m+N_n)]  
V_invS<-V_inv[(N_m+N_n+1):(N_m+N_n+N_s),(N_m+N_n+1):(N_m+N_n+N_s)]

y0<-(mean(yN)*N_n+N_m*mean(yM)+mean(yS)*N_s)/n
yN0<-yN-y0
yM0<-yM-y0
yS0<-yS-y0

mean_M <- colMeans(genoM)  
mean_N <- colMeans(genoN)
mean_S <- colMeans(genoS)
global_mean <-(mean_N *N_n+N_m*mean_M+mean_S*N_s)/n
genoM0 <- sweep(genoM, 2, global_mean)
genoN0 <- sweep(genoN, 2, global_mean)
genoS0<- sweep(genoS, 2, global_mean)

aN<-t(genoN0)%*%V_invN%*%yN0
aM<-t(genoM0)%*%V_invM%*%yM0
aS<-t(genoS0)%*%V_invS%*%yS0
a<-aN+aM+aS
cN <- colSums(genoN0^2)  
cM <- colSums(genoM0^2)
cS <- colSums(genoS0^2)
c<-cN+cM+cS

set.seed(2025)
n_snp <- ncol(genoN0)
random_integers <- sample(seq_len(n_snp), 100, replace = FALSE)
print(random_integers)   
genoN1<-genoN0[,random_integers]
genoM1<-genoM0[,random_integers]
genoS1<-genoS0[,random_integers]
bN <- colSums(genoN1 * (V_invN %*% genoN1))  
bM<- colSums(genoM1 * (V_invM %*% genoM1))
bS <- colSums(genoS1 * (V_invS %*% genoS1))
b1<-bN+bM+bS
cN1 <- colSums(genoN1^2)  
cM1 <- colSums(genoM1^2)
cS1 <- colSums(genoS1^2)
c1<-cN1+cM1+cS1
gamma=mean(b1/c1)
print(gamma)

var<-1/(c*gamma)
se<-sqrt(var)
beta<-var*a
head(beta)
z=beta/se
snp<-fread('snp200000.txt',header=FALSE)
snp<-snp$V1
final <- data.frame(
  SNP=snp,
  beta=beta,
  var =var
)
final$pval = 2 * pnorm(abs(z), lower.tail = FALSE)
head(final)
oksnp1<-subset(final,final$pval<2.5e-7)
dim(oksnp1)
causalsnp<-fread("causalsnp15w.txt",header=F)
sum(oksnp1$SNP %in% causalsnp$V1)

t_end_test =  proc.time()
cat(" test took\n")
print(t_end_test - t_begin_test)
gc()
end_mem <- mem_used()
memory_used <- end_mem - start_mem
memory_used_GB <- memory_used/ (1024^3)
print(paste("test Memory used:", memory_used_GB, "GB"))
memory_used2 <- end_mem - initial_mem
memory_used_GB <- memory_used2 / (1024^3)
print(paste("Total Memory used:", memory_used_GB, "GB"))


write.table(final, file = "resultquan15w20w.txt", row.names = FALSE, col.names=TRUE, quote = FALSE)



