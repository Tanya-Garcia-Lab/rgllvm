
#################################
## Main function to run method ##
#################################

rgllvm <- function(y,x,z,
                   n,m,p,
                   beta0=rep(0,p)){

  ###############################################
  ## Set parameters for where to store results ##
  ###############################################
  p1 <- p
  maxm <- max(m)
  lb <- p1

  ##############
  ## Set tolerance for
  eps=0.001
  tol=1e-6

  ################################
  ## set terms in common module ##
  ################################
  set_dim <- TRUE
  do.allocations(set_dim,n,maxm,lb)

  ####################################
  ## set n1 and m1 in common module ##
  ####################################
  store.n1m1(n,m)

  #####################################
  ## get vv_combinations and indices ##
  #####################################
  getvv.terms(n,m,maxm)

  ## Put data in format so that glmer can work
  data.set <- make.data.set(n,m,p,y,x)$data.set.out

  ## Set beta parameter
  beta.new <- beta0
  beta.init <- beta.new

  ## Solve semipar estimator in f90
  data.f90 <- make.data.f90(y,x)

  est.f90 <- estimation.logistic(p1,beta.init,data.f90$ynew,
                                   data.f90$znew)

  if(est.f90$eflag==1){
      beta.est <- est.f90$betaest
      beta.var <- est.f90$var
  } else {
    beta.est <- NULL
    beta.var <- NULL
  }

  ## get mle estimate
  mle.out <- get.mle(data.set,nAGQ=20)
  beta.mle <- mle.out$betas
  beta.mle.var <- as.matrix(mle.out$betas.var)
  sigma2.mle <- mle.out$sigma2
  sigma2.mle.sd <- mle.out$sigma2.sd

  ## get glmmPQL estimate
  pql.out <- get.pql(data.set)
  beta.pql <- pql.out$betas
  beta.pql.var <- as.matrix(pql.out$betas.var)


  ################################
  ## set terms in common module ##
  ################################
  set_dim <- FALSE
  do.allocations(set_dim,n,maxm,lb)



  return(list(beta.est=beta.est,beta.var=beta.var,
       beta.mle=beta.mle,beta.mle.var=beta.mle.var,
       beta.pql=beta.pql,beta.pql.var=beta.pql.var,
       sigma2.mle=sigma2.mle,sigma2.mle.sd=sigma2.mle.sd))
}

###############
## Libraries ##
###############
#library(lme4)  ## for normal-based MLE
#library(rootSolve)  ## for solving root of an equation
#library(xtable) ## for LaTeX Table
#library(MASS) ## for glmmPQL

############################
## functions used for F90 ##
############################

do.allocations <- function(set_dim,nset,maxm0,lb0){
  ##storage for f90
  storage.mode(set_dim) <- "logical"
  storage.mode(nset) <- "integer"
  storage.mode(maxm0) <- "integer"
  storage.mode(lb0) <- "integer"

  out <- .Fortran("do_allocations",set_dim,nset,maxm0,lb0,
                  PACKAGE="rgllvm")
}


## function to store n1 and m1
store.n1m1 <- function(n10,m10){
  ## storage for f90
  storage.mode(n10) <- "integer"
  storage.mode(m10) <- "integer"

  out <- .Fortran("store_n1m1",n10,m10,
                  PACKAGE="rgllvm")
}



##############################
## Function to set vv_terms ##
##############################
getvv.terms <- function(n,m,maxm){
  ## storage for f90
  storage.mode(n) <- "integer"
  storage.mode(m) <- "integer"
  storage.mode(maxm) <- "integer"

  out <- .Fortran("getvv_terms",n,m,maxm,
                  PACKAGE="rgllvm")
}


#########################################
## function to do estimation procedure ##
#########################################
estimation.logistic <- function(lb,betat,ynew,z){
  ##storage for f90
  storage.mode(betat) <- "double"
  storage.mode(ynew) <- "double"
  storage.mode(z) <- "double"


  eflag <- 0
  betaest <- array(0,dim=lb,dimnames=list(paste("lb",1:lb,sep="")))
  var <- array(0,dim=c(lb,lb),dimnames=list(paste("lb",1:lb,sep=""),
      paste("lb",1:lb,sep="")))
  var_beta <-  array(0,dim=lb,dimnames=list(paste("lb",1:lb,sep="")))

  storage.mode(eflag) <- "integer"
  storage.mode(betaest) <- "double"
  storage.mode(var) <- "double"
  storage.mode(var_beta) <- "double"

  out <- .Fortran("estimation_logistic",betat,ynew,z,eflag=eflag,
      betaest=betaest,var=var,var_beta=var_beta,
      PACKAGE="rgllvm")
  eflag <- out$eflag
  betaest <- out$betaest
  var <- out$var
  var_beta <- out$var_beta
  list(eflag=eflag,betaest=betaest,var=var,var_beta=var_beta)
}



##############################################
## function to make data compatible for f90 ##
##############################################
make.data.f90 <- function(y,z){
  ynew <- y
  ynew[is.na(ynew)] <- 0

  znew <- z
  znew[is.na(znew)] <- 0

  list(ynew=ynew,znew=znew)
}



###############################################
## Function to  make data suitable for glme4 ##
###############################################

make.data.set <- function(n,m,p,y,z){
  ## vector indicating which family each observation belongs
  family <- rep(1:n,times=m)

  ## make  y into a vector of values
  new.y <- as.vector(t(y))
  new.y <- as.vector(na.omit(new.y))

  ## make z into a matrix of values
  z.values <- NULL
  for (s in 1:p){
    new.z <- z[,,s]
    new.z <- as.vector(t(new.z))
    new.z <- as.vector(na.omit(new.z))
    z.values <- cbind(z.values,new.z)
  }

  data.set.out <- data.frame(cbind(family,new.y,z.values))
  colnames(data.set.out) <- c("family","Y",paste("Z",1:p,sep=""))
  list(data.set.out=data.set.out)
}


######################################
## Function to get normal-based MLE ##
######################################

#' @import lme4
#' @import stats
get.mle <- function(data.set,nAGQ=20){

  ##Y <- data.set$Y
  ##family <- data.set$family
  ##rm.columns <- as.vector(na.omit(match(colnames(data.set),c("family","Y"))))
  ##other.data <- data.frame(data.set[,-rm.columns])
  ##colnames(other.data) <- paste("Z",1:ncol(other.data),sep="")
  ##data.use <- data.frame(Y,family,other.data)
  ##fm <- glmer(Y ~ -1+ . +(1|family),
  ##            data=data.use,
  ##            family=binomial,nAGQ=20)


  znam <- paste("Z",1:(ncol(data.set)-2),sep="")
  znam <- c(znam,"(1|family)")
  fmla <- stat::as.formula(paste("Y~-1 +",paste(znam,collapse="+")))
##  fmla <- as.formula(paste("Y~",paste(znam,collapse="+")))

  fm <- lme4::glmer(fmla,data=data.set,family=binomial,nAGQ=nAGQ,
     	control=lme4::glmerControl(optimizer="bobyqa"))



  ## Extract needed information for fixed effects
  Vcov <- stats::vcov(fm, useScale = FALSE)
  v.delta <- Vcov
  betas <- lme4::fixef(fm)
  se <- sqrt(diag(Vcov))


##  zval <- betas / se
##  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
##  out.fixed <- cbind(betas, se, zval, pval)

  ## Extract needed information for random effects (returns Std. deviation, not std. error!)
  ##s2 <- VarCorr(fm)$`data.set$family`
  ##ntmp <- nrow(ranef(fm)$`data.set$family`)
  ##CI <- (ntmp-1) * s2/  qchisq(c(0.975, 0.025), df = ntmp - 1)
  ##out.rand <- c(s2,CI)
  ##names(out.rand) <- c("variance","CI-lo","CI-hi")

  ## Easier way to extract information for random effects
  ##remat <- summary(fm)@REmat
  ##s2 <- as.numeric(remat[,"Variance"])
  ##s2.sd <- as.numeric(remat[,"Std.Dev."])
  s2 <-0
  s2.sd <-0

  list(betas=betas,betas.var=v.delta,sigma2=s2,sigma2.sd=s2.sd)
  ##list(out.fixed=out.fixed,out.rand=out.rand)
}



##########################
## Function for glmmPQL ##
##########################

#' @import stats
#' @import MASS
#' @importFrom lme4 fixef
get.pql <- function(data.set){

  znam <- paste("Z",1:(ncol(data.set)-2),sep="")
  fmla <- stats::as.formula(paste("Y~-1 +",paste(znam,collapse="+")))

  fm <- MASS::glmmPQL(fixed=fmla,random=~1|family,
     		data=data.set,family=binomial,niter=100,verbose=FALSE)

  ## Extract needed information for fixed effects
  Vcov <- stats::vcov(fm, useScale = FALSE)
  v.delta <- Vcov
  betas <- lme4::fixef(fm)
  se <- sqrt(diag(Vcov))


##  zval <- betas / se
##  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
##  out.fixed <- cbind(betas, se, zval, pval)

  ## Extract needed information for random effects (returns Std. deviation, not std. error!)
  ##s2 <- VarCorr(fm)$`data.set$family`
  ##ntmp <- nrow(ranef(fm)$`data.set$family`)
  ##CI <- (ntmp-1) * s2/  qchisq(c(0.975, 0.025), df = ntmp - 1)
  ##out.rand <- c(s2,CI)
  ##names(out.rand) <- c("variance","CI-lo","CI-hi")

  ## Easier way to extract information for random effects
  ##remat <- summary(fm)@REmat
  ##s2 <- as.numeric(remat[,"Variance"])
  ##s2.sd <- as.numeric(remat[,"Std.Dev."])
  s2 <-0
  s2.sd <-0

  list(betas=betas,betas.var=v.delta,sigma2=s2,sigma2.sd=s2.sd)
  ##list(out.fixed=out.fixed,out.rand=out.rand)
}


#############################################
## Functions to get our proposed estimator ##
#############################################

## function to see which vv.combinations to include
select.vv <- function(vv.combinations,wi){
  tmp <- apply(vv.combinations,1,sum)
  index <- which(tmp <=wi & tmp >=(wi-1))
  return(as.matrix(vv.combinations[index,]))
}

## Compute E(V_i|W_i,Z_i;beta)

expect.v <- function(beta,y,z,m,n,p,tol=1e-6){
  ev <- array(NA,dim=c(n,max(m)),
              dimnames=list(1:n,paste("m",1:max(m),sep="")))

  for(i in 1:n){
    Ai.inv <-	get.Ainv(m[i])
    wi <- sum(y[i,1:m[i]])

    m.tmp <- m[i]-1
    mat.tmp <- matrix(rep(c(0,1),each=m.tmp),ncol=m.tmp,nrow=2,byrow=TRUE)
    vv.combinations <- as.matrix(do.call(`expand.grid`,as.data.frame(mat.tmp)))

    z.tmp <- t(z[i,1:m[i],])
    theta <- t(z.tmp) %*% beta

    if(wi==0){
      ev[i,2:m[i]] <- 0
    } else if(wi==m[i]){
      ev[i,2:m[i]] <-1
    } else {
      numerator <- 0
      denominator <- 0
      vv.new <- select.vv(vv.combinations,wi)

      for(u in 1:nrow(vv.new)){
        V <- c(wi,vv.new[u,])
        main.tmp <- exp( as.vector(theta) %*% Ai.inv %*% V)
        numerator <- numerator + vv.new[u,] * main.tmp
        denominator <- denominator + main.tmp
      }
      ##if(abs(denominator) <tol){
      ##  ev[i,2:m[i]] <- 0
      ##} else {
      ev[i,2:m[i]] <- numerator/denominator
      ##}
    }
  }
  ev[,1] <- 0

  return(ev)
}

## function to compute A_i^{-1}
get.Ainv <- function(mi){
  Ai.inv <- diag(1,mi)
  Ai.inv[1,-1] <- -1
  return(Ai.inv)
}

## Create terms in score vector
seff.terms <- function(beta,y,z,m,n,p){
  ## get E(V_i|W_i,Z_i;beta)
  ev <- expect.v(beta,y,z,m,n,p)

  ## place to store terms in the score vector
  score.out <- array(0,dim=c(p,n),
                     dimnames=list(paste("p",1:p,sep=""),1:n))

  for(i in 1:n){
    ##print(i)
    ##m.tmp <- m[i]-1
    ##z1.tmp <- matrix(rep(z[i,1,],m.tmp),ncol=p,nrow=m.tmp,byrow=TRUE)
    ##colnames(z1.tmp) <- paste("p",1:p,sep="")

    ##z.subtract <- (z[i,2:m[i],]-z1.tmp)
    ##tmpz <-  z.subtract
    ##tmpv <- (y[i,2:m[i]]-ev[i,2:m[i]])
    ##tmp.prod <- sweep(t(tmpz),MARGIN=2,tmpv,`*`)
    ##if(nrow(tmp.prod)==1){
    ##	score[,i] <- tmp.prod
    ##} else {
    ##	score[,i] <- apply(tmp.prod,1,sum)
    ##}


    ## Attempt 2
    Ai.inv <-	get.Ainv(m[i])
    z.tmp <- t(z[i,1:m[i],])

    Vi <- y[i,]
    Vi[1] <- 0
    score.out[,i] <- z.tmp %*% Ai.inv %*% (Vi[1:m[i]]-ev[i,1:m[i]])


  }
  list(score.out=score.out)
}



## Create score vector
make.seff <- function(y,z,m,n,p){
  function(beta){
    get.terms <- seff.terms(beta,y=y,z=z,m=m,n=n,p=p)
    score <- get.terms$score.out
    out <- apply(score,1,sum)
    return(out)
  }
}

## Compute variance of beta using semipar approach
get.var <- function(beta,y,z,m,n,p){
  get.terms <- seff.terms(beta,y,z,m,n,p)
  score <- get.terms$score.out

  out <- 0

  for(i in 1:n){
    tmp <- as.matrix(score[,i])
    out <- out + tmp %*% t(tmp)
  }

  ## Compute B^{-1}
  out <- solve(out)

  return(out)
}


## function to convert alphac to beta
get.beta <- function(alphac){
  beta <- alphac[-1]-alphac[1]
  return(beta)
}

## for testing
myfunc1 <- function(theta){
  num <- 1+exp(theta+1)
  den <- 1+exp(theta)
  return(log(num/den))
}

myfunc2 <- function(theta){
  num <- exp(theta)
  den <- 1+ exp(theta)
  return(num/den)
}



#####################
## ANALZZE RESTULS ##
#####################

## out : results from simulations()
#' @import stats
ana.results <- function(nsimu,p,beta.est, beta.var,
                        beta.mle, beta.mle.var,
			beta.pql, beta.pql.var,
		        sigma2.mle,sigma2.mle.sd,beta0,
			real.data=FALSE){

 if(real.data==FALSE){
  betat <- beta0

  haus.out <- array(0,dim=c(nsimu,8),
			dimnames=list(1:nsimu,c("MLEstatistic","k","pvalue","reject",
					"PQLstatistic","k","pvalue","reject")))

  for(s in 1:nsimu){
    ## get Hausman test statistic
    haus.out[s,1:4] <- hausman.test(beta.est[s,],beta.mle[s,],
  	       		beta.var[s,,],beta.mle.var[s,,])
    haus.out[s,5:8] <- hausman.test(beta.est[s,],beta.pql[s,],
  	       		beta.var[s,,],beta.pql.var[s,,])
  }

  ## Mean of hausman.test
  haus.result <- apply(haus.out,2,mean)

  ## Mean of fixed effects
  beta.mean <- apply(beta.est,2,mean)
  beta.mean.mle <- apply(beta.mle,2,mean)
  beta.mean.pql <- apply(beta.pql,2,mean)

  ## Bias of fixed effects
  bias.beta <- beta.mean-betat
  bias.beta.mle <- beta.mean.mle-betat
  bias.beta.pql <- beta.mean.pql-betat

  ## Empirical variance of fixed effects
  emp.var.beta <- apply(beta.est,2,var)
  emp.var.beta.mle <- apply(beta.mle,2,var)
  emp.var.beta.pql <- apply(beta.pql,2,var)

  ## Estimated variance of fixed effects
    beta.var.est <- t(apply(beta.var,1,diag))
    beta.mle.var.est <- t(apply(beta.mle.var,1,diag))
    beta.pql.var.est <- t(apply(beta.pql.var,1,diag))

  ## Average of estimated variance of fixed effects
  beta.var.mean <- diag(apply(beta.var,c(2,3),mean))
  beta.mle.var.mean <- diag(apply(beta.mle.var,c(2,3),mean))
  beta.pql.var.mean <- diag(apply(beta.pql.var,c(2,3),mean))


  ## Function to get Coverages
  get.coverages <- function(beta.est,betat,beta.var.est,nsimu){
    tmp1 <- abs(sweep(beta.est,MARGIN=2,betat,`-`))
    tmp2 <- tmp1/sqrt(beta.var.est)
    ##tmp2 <- sweep(tmp1,MARGIN=2,sqrt(beta.var.est),`/`)
    tmp3 <- apply(tmp2,1,stats::pnorm)<0.975

    if(length(betat)==1) {
      cov.prob <- sum(tmp3)/nsimu
    } else {
      cov.prob <- apply(tmp3,1,sum)/nsimu
    }
  }

  ## Get coverages
  cov.beta <- get.coverages(beta.est,betat,beta.var.est,nsimu)
  cov.beta.mle <- get.coverages(beta.mle,betat,beta.mle.var.est,nsimu)
  cov.beta.pql <- get.coverages(beta.pql,betat,beta.pql.var.est,nsimu)

  ## Organize results
  semi.out <- rbind(bias.beta,emp.var.beta,beta.var.mean,cov.beta)
  mle.out  <- rbind(bias.beta.mle,emp.var.beta.mle,beta.mle.var.mean,cov.beta.mle)
  pql.out  <- rbind(bias.beta.pql,emp.var.beta.pql,beta.pql.var.mean,cov.beta.pql)

  rownames.tmp <- c("bias","emp var","est var","coverage")
  rownames(semi.out) <- rownames.tmp
  rownames(mle.out) <- rownames.tmp
  rownames(pql.out) <- rownames.tmp

  colnames(semi.out) <- paste("beta",1:ncol(semi.out),sep="")
  colnames(mle.out) <- paste("beta",1:ncol(mle.out),sep="")
  colnames(pql.out) <- paste("beta",1:ncol(pql.out),sep="")

 list(semi.out=semi.out,mle.out=mle.out,pql.out=pql.out,
	haus.out=haus.result)
 } else {
   ###############
   ## Real data ##
   ###############

   ## parameter estimates
   beta.semi <- t(beta.est)
   beta.mle <- t(beta.mle)
   beta.pql <- t(beta.pql)

   ## estimated variances
    beta.var.est <- t(apply(beta.var,1,diag))
    beta.mle.var.est <- t(apply(beta.mle.var,1,diag))
    beta.pql.var.est <- t(apply(beta.pql.var,1,diag))


   ## 95\% confidence interval
   get.ci <- function(beta,var.tmp,alpha=0.05){
     hi <- beta + qnorm(1-alpha/2)*sqrt(var.tmp)
     lo <- beta - qnorm(1-alpha/2)*sqrt(var.tmp)
     list(hi=hi,lo=lo)
   }
   semi.ci <- get.ci(beta.semi,beta.var.est)
   mle.ci  <- get.ci(beta.mle,beta.mle.var.est)
   pql.ci  <- get.ci(beta.pql,beta.pql.var.est)

   semi.out <- rbind(beta.semi,beta.var.est,semi.ci$lo,semi.ci$hi)
   mle.out <- rbind(beta.mle,beta.mle.var.est,mle.ci$lo,mle.ci$hi)
   pql.out <- rbind(beta.pql,beta.pql.var.est,pql.ci$lo,pql.ci$hi)

   rownames.tmp <- c("est","var","lo","hi")
   rownames(semi.out) <- rownames.tmp
   rownames(mle.out) <- rownames.tmp
   rownames(pql.out) <- rownames.tmp

   colnames(semi.out) <- paste("beta",1:ncol(semi.out),sep="")
   colnames(mle.out) <- paste("beta",1:ncol(mle.out),sep="")
   colnames(pql.out) <- paste("beta",1:ncol(pql.out),sep="")

   ##################
   ## Hausman test ##
   ##################
     haus.out <- c(hausman.test(beta.semi,beta.mle,beta.var[1,,],beta.mle.var[1,,]),
   	       hausman.test(beta.semi,beta.pql,beta.var[1,,],beta.pql.var[1,,]))


   list(semi.out=semi.out,mle.out=mle.out,pql.out=pql.out,
	haus.out=haus.out)
  }

}

#' @importFrom Matrix rankMatrix
#' @importFrom stats pchisq
hausman.test <- function(beta.semi,beta.mle1,var.semi,var.mle){
   beta.diff <- beta.semi-beta.mle1
   var.diff  <- var.semi-var.mle
   k <- as.numeric(Matrix::rankMatrix(as.matrix(var.diff)))
   test.statistic <- as.matrix(beta.diff) %*%  solve(var.diff)
   test.statistic <- test.statistic %*% t(as.matrix(beta.diff))
   p.value <- 1 - stats::pchisq(test.statistic, k)
   reject.h0 <- p.value < 0.05

   haus.out <- cbind(test.statistic,k,p.value,reject.h0)
   colnames(haus.out) <- c("Statistic","df","p-value","reject_h0")
   rownames(haus.out) <- "Hausman.test"
   return(haus.out)
}





