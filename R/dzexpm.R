dzexpm <-
function(y_ful,x_ful,n_ful,n,u1,u2,theta,p,iter,loopiter){

######################################################
###########      function declearation 
######################################################
 g2 = function(mu)
 {
   temp   = mu^2/2
   s_mu   = sign(mu)
   s_mu2  = sign(mu*mu)
   s_neg  = -(s_mu - rep(1, length(mu)))/2
   s_pos  =  (s_mu + rep(1, length(mu)))/2
   s_zero = rep(1, length(mu)) - s_mu2
   temp   = (s_neg -0.5*s_zero)*(0.5 + pgamma(temp, shape=1.5, scale=1)/2) + (s_pos +   0.5*s_zero)*(0.5-pgamma (temp, shape=1.5, scale=1)/2)
   temp   = temp - (2*mu)*dnorm(mu, mean = 0, sd = 1) + (mu^2)*pnorm(-mu, mean = 0, sd = 1)
   temp   = temp/(1+(mu^2))
   if(temp < 0.0000001)
   {
     temp = 0.0000001
   }
   return(temp)  
 }
#######################################################
 s_1 = function(s0, mu, p)
 {
   temp = (s0^2)/(1/p-1) 
   temp = temp*(1/g2(mu) -1)
   if( temp < 0)
   {
     temp = 0
   }
   return(sqrt(temp))  
 }
####################################################
 diagonal = function(s0, s1, eps) #tested
 {
   d     = diag((sign(eps)*(s0-s1)/2 + (s0+s1)/2))
   n_neg = sum(-1*sign(eps) +1)/2
   return(list(d=d, n_neg=n_neg ))
 }
####################################################
 update_beta = function(y, X, Sigma_inv)
  {
  temp    = t(X)%*%Sigma_inv
  Sigma_b = solve(temp%*%X)
  mu_b    = Sigma_b%*%temp%*%y
  beta    = t(chol(Sigma_b))%*%rnorm(n= length(X[1,]), mean=0, sd=1) + mu_b
  return(c(beta))
  }
####################################################
 update_mu = function(mu, st_mu, s0, p, eps, cova_inv, C, n_neg, mu_0, sigma2_mu, ic)
 {
   mu_new = rnorm(1, mean=mu, sd=st_mu) #normal proposal
   s1_new = s_1(s0=s0,mu=mu_new, p=p)
   if(s1_new==0)
   {
     s1_new=0.0000001
   }
   else if(500000 < s1_new)
   {
     s1_new = 500000
   }
   out                 = diagonal(s0=s0, s1=s1_new, eps=c(eps))
   d_new               = out$d
   temp                = diag(sqrt(diag(d_new)))
   cova_inv_new        = solve((temp%*% C)%*%temp)  
   log_det             = n_neg*(log(s1/s1_new))
   temp1               = (mu*mu - mu_new*mu_new)/(2*sigma2_mu)
   temp2               = mu_0*(mu_new - mu)/sigma2_mu
   likeratio           = 0.5*(sum(eps*((cova_inv - cova_inv_new)%*%eps)) + log_det) 
   ratio               = temp1 + temp2 + likeratio  
   if (0 >= ratio)
   {
     u = - rexp(1, 1)
     if( u < ratio)
     {
       mu       = mu_new
       ic       = ic + 1
       d        = d_new
       s1       = s1_new
       cova_inv = cova_inv_new
     }
   }  
   else
   {
     mu       = mu_new
     ic       = ic + 1
     d        = d_new
     s1       = s1_new
    cova_inv = cova_inv_new
   }
   return(list(mu=mu, ic=ic, d=d, s1=s1, cova_inv = cova_inv))
 }
#####################################################################
 update_s0 = function(s0, st_s0, mu, p, eps, cova_inv, C, n_neg, alpha_s0, beta_s0, ic)
 {
   s0_new  = rnorm(n=1, mean=log(s0), sd=st_s0) #log normal proposal
   s0_new  = exp(s0_new)
   s1_new  = s_1(s0_new, mu, p)
   if(s1_new==0)
   {
     s1_new = 0.0000001
   }
   else if(500000 < s1_new)
   {
     s1_new = 500000
   }
   out                 = diagonal(s0, s1_new, c(eps))
   d_new               = out$d
   temp                = diag(sqrt(diag(d_new)))
   cova_inv_new        = solve((temp%*% C)%*%temp)
   log_det             = n_neg*(log(s1/s1_new))
   temp1               = (s0 - s0_new)/beta_s0
   temp2               = log(s0_new/s0)*alpha_s0       #gamma prior
   likeratio           = 0.5*(sum(eps*((cova_inv - cova_inv_new)%*%eps)) + log_det) 
   ratio               = temp1 + temp2 + likeratio  
   if (0 >= ratio)
   {
     u = - rexp(1, 1)
     if( u < ratio)
     {
       s0       = s0_new
       ic       = ic + 1
       d        = d_new
       s1       = s1_new
       cova_inv = cova_inv_new
     }
   }  
   else
   {
     s0       = s0_new
     ic       = ic + 1
     d        = d_new
     s1       = s1_new
     cova_inv = cova_inv_new
   }
   return(list(s0=s0, ic=ic, d=d, s1=s1, cova_inv = cova_inv))
 }
##################################################################
 update_tau=function(eps_e, beta_tau, alpha_tau)
 {  
    alpha = (alpha_tau + length(eps_e)/2)
    beta  = beta_tau   + sum(eps_e*eps_e)/2
    tau = 1/rgamma(n=1, shape= alpha, rate = beta) 
    return(tau)
 }
##################################################################
 update_eps_e = function (tau, Sigma_inv, Z)
 {
   n           = length(Z)
   Sigma_eps_e = solve(Sigma_inv + (1/tau)*diag(1,n))
   mu_eps_e    = Sigma_eps_e %*%(Sigma_inv%*%Z)
   eps_e       = t(chol(Sigma_eps_e))%*%rnorm(n) + mu_eps_e
   return(eps_e)
 }
##################################################################
 rho = function(u1, u2, theta)
 {
   n = length(u1)
   corr = matrix(rep(1, (n*n)), ncol=n)
   for(i in 1:(n-1))
   {
     itemp = i + 1
     for(j in itemp:n)
     {
       temp      = sqrt( (u1[i]-u1[j])* (u1[i]-u1[j])+(u2[i]-u2[j])* (u2[i]-u2[j]))
       corr[i,j] = exp(-theta*temp)
       corr[j,i] = corr[i,j]
     }
   }
   return(corr)
 }
###################################################################
 dic        = function(n,eps,sig_inv,det_dic)
 {      dic1   =-0.5*n*log(2*pi)
        dic2   = det_dic
        dic3   =-0.5*(t(eps)%*%sig_inv)%*%(eps)
        dictemp=dic1+dic2+dic3
        dic    =-2*dictemp
        return(dic)
 }
######################################################
#######        default set for the parameters
######################################################
 st_mu     = 0.4     #gives % acceptance ratio for the fixed parameters 
 st_s0     = 0.001   #gives %acceptance ratio for the fixed parameters
  
 alpha_s0  = 1 #shape, makes the prior mean 1 and variance inf
 beta_s0   = 1 #scale
 
 alpha_tau = 1
 beta_tau  = 1
 
 mu_0      = 0
 sigma2_mu = 1
 rmu=0
 rs0=0
######################################################
#######        default set for the results' indecies
######################################################
 
 #iter = 100
 #loopiter = 1000
 ic_vec = bias = mpe = SD= 0  #### index for the results
 ic_mu = ic_s0 = 0
 
 s0_dic=0
 s1_dic=0

 #p=0.1   ## assumed in our model
   mu          =  0
   tau         =  0.5       
   s0          =  s0_true = 2     ## we mean the sigma_0 squared, so be careful
   s1          =  s_1(s0, mu, p)  ## we mean the sigma_1 squared, this is assumed in our  model
   beta        =  matrix(c(1,1),nrow=2)
#################################################
###          index and structure clearation
#################################################
##data input
  set.seed(101570)
  n_pred = n_ful-n
  ind    = sample(1:n_ful, n)
  n_pre  = rep(0,n_pred)
  j = 1
  for(i in 1:n_ful)
  {
     if(prod(rep(i,n)-ind)!=0)
    { n_pre[j] = i; j = j + 1 }
  }
   
  y  = y_ful[ind]
  x  = cbind(matrix(1,nrow=n,ncol=1),x_ful[ind])
  u3 = rbind(matrix(u1[ind],ncol=1),matrix(u1[n_pre],ncol=1))
  u4 = rbind(matrix(u2[ind],ncol=1),matrix(u2[n_pre],ncol=1))

  C       =  rho(u3,u4, theta)
  C11     =  C[1:n,1:n]
  C11_inv =  solve(C11)
  C22     =  C[(n+1):n_ful,(n+1):n_ful]
  C21     =  C[(n+1):n_ful,1:n]
  Ctemp   =  C21%*%C11_inv
  Cpred   =  C22 - Ctemp%*%t(C21)
  C_pred_sqrt   = t(chol(Cpred))

   x_pred       = cbind(matrix(1,nrow=n_pred,ncol=1),x_ful[n_pre])
   y_pred_true  = y_ful[n_pre]     ## holdout points
   eps          = y- x%*%beta
   out          = diagonal(s0, s1, c(eps))
   d            = out$d
   d_sqrt       = diag(sqrt(diag(d))) 
   d_sqrt_inv   = diag(1/sqrt(diag(d)))

   ic_mat= c(rep(0, n_pred))
   cov   = c(rep(0, n_pred))
#################################################
###          gear up for MCMC               #####
#################################################
   y_pred_mat   = matrix(rep(0, iter*n_pred), nrow=n_pred)
   eps_e        = rep(0, n)

   beta_mat     = matrix(rep(0,iter*2),ncol=iter)
   eps_e_mat    = matrix(rep(0,iter*n),ncol=iter)
   mu_mat       = c(0,iter)
   s0_mat       = c(0,iter)
   tau_mat      = c(0,iter)
   Dthetabar_mat= c(0,iter)
   Dbar_mat     = c(0,iter) 
   mu_t            = matrix(0,nrow=loopiter,ncol=iter)
   tau_t           = matrix(0,nrow=loopiter,ncol=iter)
   s0_t            = matrix(0,nrow=loopiter,ncol=iter)
   beta_t1         = matrix(0,nrow=loopiter,ncol=iter)
   beta_t2         = matrix(0,nrow=loopiter,ncol=iter)
#################################################
###          iteration for MCMC             #####
#################################################
  for(i in 1:iter){
     for(j in 1:loopiter){ 
       cova_inv    = (d_sqrt_inv%*%C11_inv)%*%d_sqrt_inv
       y_tilde     = y - eps_e
       beta        = update_beta(y=y_tilde, X=x, Sigma_inv=cova_inv)
       beta_t1[j,i]= beta[1]
       beta_t2[j,i]= beta[2]
        
       eps         = y -x%*%beta - eps_e
       out         = diagonal(s0=s0, s1=s1, eps=c(eps))
       d           = out$d
       n_neg       = out$n_neg
       d_sqrt      = diag(sqrt(diag(d)))
       d_sqrt_inv  = diag(1/sqrt(diag(d)))
       cova_inv    = (d_sqrt_inv%*%C11_inv)%*%d_sqrt_inv

       out         = update_mu(mu=mu, st_mu=st_mu, s0=s0, p=p, eps=eps, cova_inv=cova_inv,   C=C11,  n_neg=n_neg, mu_0=mu_0, sigma2_mu=sigma2_mu, ic= ic_mu)
       if(ic_mu != out$ic){
         mu        = out$mu
         ic_mu     = out$ic
         s1        = out$s1
         d         = out$d
         cova_inv  = out$cova_inv
         d_sqrt    = diag(sqrt(diag(d)))
         d_sqrt_inv = diag(1/sqrt(diag(d)))
       }
       rm(out)
       mu_t[j,i]  = mu

       tau        = update_tau(eps_e=eps_e, beta_tau= beta_tau, alpha_tau = alpha_tau) 
       tau_t[j,i] = tau
       Z          = y - x%*%beta
       eps_e      = update_eps_e(tau=tau, Z=Z, Sigma_inv =cova_inv)  
       eps        = y - x%*%beta - eps_e
       out        = diagonal(s0=s0, s1=s1, eps=c(eps))
       d          = out$d     
       cova_inv   = (d_sqrt_inv%*%C11_inv)%*%d_sqrt_inv
       out        = update_s0(s0=s0, st_s0=st_s0, mu=mu, p=p, eps=eps, cova_inv=cova_inv,   C=C11,  n_neg=n_neg, alpha_s0=alpha_s0, beta_s0=beta_s0, ic=ic_s0)
       if(ic_s0 != out$ic){
         s0       = out$s0
         ic_s0    = out$ic
         s1       = out$s1
         d        = out$d
         cova_inv = out$cova_inv
         d_sqrt   = diag(sqrt(diag(d)))
         d_sqrt_inv = diag(1/sqrt(diag(d)))
       }
       rm(out)
       s0_t[j,i]    = s0      
    }  ##loopiter ends
##############################################################
############ results for the iter i 
##############################################################
       eps_dic       = y -x%*%beta - eps_e
       out_dic       = diagonal(s0=s0, s1=s1, eps=c(eps_dic))
       d_dic         = out_dic$d
       n_neg_dic     = out_dic$n_neg
       d_sqrt_inv_dic= diag(1/sqrt(diag(d_dic)))
       log_det_dic_1 = (n_neg_dic)*(log(s1))
       log_det_dic_2 = (n-n_neg_dic)*(log(s0))
       log_det_temp  = log_det_dic_1+log_det_dic_2
       det_dic       = (-0.5)*log_det_temp+0.5*log(det(C11_inv)) 
       rm(out_dic)
       cova_inv_dic      = (d_sqrt_inv_dic%*%C11_inv)%*%d_sqrt_inv_dic
       Dbar_mat[i]       = dic(n=n,eps=eps_dic,sig_inv=cova_inv_dic,det_dic=det_dic)
###############################################################
###############     parameters' results    
###############################################################
     beta_mat[,i]      = beta 
     mu_mat[i]         = mu
     s0_mat[i]         = s0 
     eps_e_mat[,i]     = eps_e
     tau_mat[i]        = tau
  
     eps               = diag(d_sqrt_inv)*(y - x%*%beta - eps_e) - mu
     eps_pred          = C_pred_sqrt%*%rnorm(n=n_pred, mean=0, sd=1) + Ctemp%*%eps + mu
     out               = diagonal(s0, s1, c(eps_pred))
     d_pred            = out$d
     eps_pred          = sqrt(diag(d_pred))*eps_pred
     y_pred_mat[,i]    = x_pred %*%beta + eps_pred + sqrt(tau)*rnorm(n_pred)
 } ##iter ends
#####################################
############       Coverage 
#####################################   
  for(i in 1: n_pred){ 
      if(y_pred_true[i] < quantile(y_pred_mat[i,],prob=c(0.975)) &&  quantile(y_pred_mat[i,],prob=c(0.025))  <y_pred_true[i] )
      {
      ic_mat[i] =ic_mat[i]+1
      }
     cov[i]=sum(100*ic_mat[i])
   }
  ic_vec    = 100*sum(ic_mat)/n_pred
#####################################
############       bias,mpe,SD 
#####################################
   temp1     = y_pred_mat - kronecker(matrix(rep(1, iter), nrow=1), y_pred_true)
   temp2     = temp1*temp1
   bias      = median(apply(temp1,1,median))
   mpe       = median(apply(temp2,1,median))
   SD        = median(apply(y_pred_mat,1,sd))
   #print(c(ic_vec, bias, mpe,SD))
#####################################
############       DIC 
#####################################  
   beta_dic1      = sum(beta_mat[1,])/iter
   beta_dic2      = sum(beta_mat[2,])/iter
   beta_dic       = matrix(c(beta_dic1,beta_dic2),nrow=2)
   mu_dic         = sum(mu_mat)/iter
   eps_dic_1      = y -x%*%beta_dic - rowMeans(eps_e_mat)
   s0_dic         = sum(s0_mat)/iter
   s1_dic         = s_1(s0_dic, mu_dic, p)
       mout_dic         = diagonal(s0=s0_dic, s1=s1_dic, eps=c(eps_dic_1))
       md_dic           = mout_dic$d
       mn_neg_dic       = mout_dic$n_neg
       md_sqrt_inv_dic  = diag(1/sqrt(diag(md_dic)))
       mlog_det_dic_1    = (mn_neg_dic)*(log(s1_dic))
       mlog_det_dic_2    = (n-mn_neg_dic)*(log(s0_dic))
       mlog_det_temp     = mlog_det_dic_1+mlog_det_dic_2
       mdet_dic          = (-0.5)*(mlog_det_temp)+0.5*log(det(C11_inv)) 
 
       rm(mout_dic)
       cova_inv_dic_1  = (md_sqrt_inv_dic%*%C11_inv)%*%(md_sqrt_inv_dic)
   Dbar         =mean(Dbar_mat)
   Dtheta_bar   =dic(n=n,eps=eps_dic_1,sig_inv=cova_inv_dic_1,det_dic=mdet_dic)   
   DIC          =2*Dbar-Dtheta_bar
#####################################
############       Acceptance rates
##################################### 
# rate_mu=ic_mu/(loopiter*iter)
# rate_s0=ic_s0/(loopiter*iter)
#####################################
############       results return 
#####################################  
data.frame(DIC=DIC,coverage=ic_vec,bias.med=mean(bias),mpe.med=mean(mpe),SD.med=mean(SD))
}
