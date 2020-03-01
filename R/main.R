
#################################
## Main function to run method ##
#################################

#' Robust generalized linear latent variable model
#'
#' Estimates parameters for repeated measures data using a generalized linear latent variable model
#' where the latent variable is a random slope.
#'
#'@param y.data (N by 2) matrix of the subject ID (first column) and repeated-measures outcome (second column).
#'The data is in long format, where each row is one time point per subject and N is the total number of repeated measures for all subjects.
#'@param x.data (N by (1+p)) matrix of the subject ID (first column) and p-many repeated-measure covariates that exert a fixed effect in the model.
#'The data is in long format where each row is one time point per subject and N is the total number of repeated measures for all subjects.
#'@param z.data (N by (1+q)) matrix of the subject ID (first column) and q-many repeated-measure covariates that have an associated random effect.
#'The data is in long format where each row is one time point per subject and N is the total number of repeated measures for all subjects. Currently,
#'only q=1 is implemented.
#'@param n scalar denoting the number of sujects in the study.
#'@param m n-dimensional vector containing the number of repeated measures per subject.
#'@param p scalar denoting the number of fixed effects in the model.
#'@param beta0 p-dimensional vector containing the initial estimate of the fixed effect parameters. Default is a p-dimensional vector of zeroes.
#'@param family character vector denoting the type of generalized linear model. Code currently only implements family="binomial"
#'for the logistic model.
#'
#'@section Details:
#'This function reports three ways to estimate parameters for repeated measures data using a generalized linear latent variable model.
#'The first is a semiparametric approach where the distribution of the latent variable is left unspecified. The second is the maximum likelihood
#'estimator which assumes the latent variable is normally distributed with mean zero.
#'The third is the penalized quasi-likihood estimator which does not assume the
#'latent variable is normally distributed, only that it has zero mean and finite variance.
#'
#'
#' @references
#' Breslow, N. E. & Clayton, D. G. (1993). Approximate inference in generalized linear mixed models.
#' Journal of the American Statistical Association 88, 9-25.
#'
#' Garcia, T.P. and Ma, Y. (2015). Optimal estimator for
#' logistic model with distribution-free random intercept.
#' Scandinavian Journal of Statistics, 43, 156-171.
#'
#' Ma, Y. & Genton, M. G. (2010). Explicit estimating equations for semiparametric generalized
#' linear latent variable models. Journal of the Royal Statistical Society, Series B 72, 475-495.
#'
#' Schall, R. (1991) Estimation in generalized linear models with random effects. Biometrika 78, 719-727.
#'
#' Wei, Y., Ma, Y., Garcia, T.P. and Sinha, S. (2019). Consistent estimator
#' for logistic mixed effect models. The Canadian Journal of Statistics, 47, 140-156. doi:10.1002/cjs.11482.
#'
#'
#'@return \code{rgllvm} returns a list containing
#' \itemize{
#'   \item{beta.est:}{fixed effect estimates from the semiparametric estimator
#'   that does not make any distributional assumptions on the latent variable.}
#'   \item{beta.var:}{estimated variances of the fixed effects from the  semiparametric
#'   estimator that does not make any distributional assumptions on the latent variable}
#'   \item{beta.mle:}{fixed effect estimates from the maximum likelihood estimator.}
#'   \item{beta.mle.var:}{estimated variances of the fixed effect estimates from the maximum likelihood estimator.}
#'   \item{beta.pql:}{fixed effect estimates from the penalized quasi-likelihood estimator.}
#'   \item{beta.pql.var:}{estimated variances of the fixed effects from the penalized quasi-likelihood estimator.}
#' }
#'@export
#'
#' @importFrom lqmm gauss.quad
#' @importFrom nleqslv nleqslv
rgllvm <- function(y.data,x.data,z.data,
                   n,m,p,
                   beta0=rep(0,p),
                   family="binomial"){


  ###########################
  ## Do checks on the data ##
  ###########################
  if(nrow(y.data)!= nrow(x.data)){
    stop("The number of rows of y.data and x.data should be the same.")
  }

  if(!is.null(z.data)){
    if(nrow(y.data)!=nrow(z.data)){
      stop("The number of row of y.data and z.data should be the same.")
    }
  }

  if(ncol(z.data)>2){
    stop("Current code only allows one random slope.")
  }

  ###############################################
  ## Set parameters for where to store results ##
  ###############################################
  maxm <- max(m)

  ######################################################################
  ## put data sets into the forms needed for estimation with f90 code ##
  ######################################################################
  y <- create.wide.form(y.data,n,m)
  x <- create.wide.form(x.data,n,m)
  if(!is.null(z.data)){
    z <- create.wide.form(z.data,n,m)
  }

    #######################################
    ## latent variable is a random slope ##
    #######################################

    ## set up variables for the analysis
    x1 <- create.list.form(x.data,n,m)
    c1 <- n
    u1 <- rep(0,c1)

    GLrule <- lqmm::gauss.quad(10,kind="hermite")
    t1 <- GLrule[[1]]
    w1 <- GLrule[[2]]

    mu <- 0; sigma <- 1 # used for MLE estimates, but not needed for our purposes.

    esteqf2 <- function(beta) {
      esteqf11=.Call('_rgllvm_esteqc', PACKAGE = 'rgllvm',y,x1,z,u1,mu,sigma,t1,w1,beta)
      esteqf11=Reduce(`+`, esteqf11)
      return (esteqf11/c1)
    }

    u2 <- allu(y,z)
    u1 <- getnewu(u2)

    ebeta <- nleqslv::nleqslv(beta0+0.1,esteqf2)

    simu1 <- ebeta[[1]]
    simu2 <- simu1

    ##To calculate the sandwich matrix of estiamted beta.
    ##sandmid is the middle part for calculating the sandwich matrix.


    temp1=.Call('_rgllvm_esteqc', PACKAGE = 'rgllvm',y,x1,z,u1,mu,sigma,t1,w1,simu2)
    temp2=rep(0,p)
    temp3=matrix(0,p,p)
    temp4=matrix(0,p,p)

    for (s in 1:c1) {
      temp2=temp1[[s]]
      temp2=as.matrix(temp2)
      temp3=temp2%*%t(temp2)
      temp4=temp4+temp3
    }


    A=matrix(0,p,p)
    temp5=matrix(0,p,p)
    for (s in 1:c1) {
        delta=simu2*0.001
        betal=simu2
        betar=simu2



        for (tt in 1:p) {
          betal[tt]=simu2[tt]-delta[tt]
          betar[tt]=simu2[tt]+delta[tt]
          yout1=.Call('_rgllvm_efficientscorefc', PACKAGE = 'rgllvm',y[s,],x1[[s]],z[s,],u1[[s]],mu,sigma,t1,w1,betal)
          yout2=.Call('_rgllvm_efficientscorefc', PACKAGE = 'rgllvm',y[s,],x1[[s]],z[s,],u1[[s]],mu,sigma,t1,w1,betar)
          A[,tt]=(yout2-yout1)/(2*delta[tt])
          betal=simu2
          betar=simu2
        }

        temp5=temp5+A

    }

    betaestvari=solve(temp5)%*%temp4%*%t(solve(temp5))

    ## for the output
    beta.est <- ebeta[[1]]
    beta.var <- betaestvari





  #####################################################################
  ## concatenate long-form data sets for estimation with MLE and PQL ##
  #####################################################################
#  data.set <- make.data.set(n,m,p,y.data,x.data,z.data)

#  ## get mle estimate
#  mle.out <- get.mle(data.set,p,z.data,nAGQ=20)
#  beta.mle <- mle.out$betas
#  beta.mle.var <- as.matrix(mle.out$betas.var)

#  ## get glmmPQL estimate
#  pql.out <- get.pql(data.set,p,z.data)
#  beta.pql <- pql.out$betas
#  beta.pql.var <- as.matrix(pql.out$betas.var)



  return(list(beta.est=beta.est,beta.var=beta.var#,
      # beta.mle=beta.mle,beta.mle.var=beta.mle.var,
      # beta.pql=beta.pql,beta.pql.var=beta.pql.var
      ))
}

################################################
## functions to transform data into wide form ##
################################################
create.wide.form <- function(x,n,m){
  subjectID <- x[,1]
  data.to.organize <- data.frame(x[,-1])
  lb <- ncol(data.to.organize)
  uniqueID <- unique(subjectID)

  if(length(uniqueID)!=n){
    stop("The number of unique subject ID's in the data does not match the sample size n.")
  }

  if(lb >1){
    out <- array(NA,dim=c(n,max(m),lb),
                     dimnames=list(1:n,paste("m",1:max(m),sep=""),1:lb))
  } else {
    out <- array(NA,dim=c(n,max(m)),
                 dimnames=list(1:n,paste("m",1:max(m),sep="")))
  }

  for(ll in 1:length(uniqueID)){
    data.tmp <- which(x[,1]==uniqueID[ll])
    if(lb >1){
      out[ll,1:length(data.tmp),] <- as.matrix(data.to.organize[data.tmp,])
    } else {
      out[ll,1:length(data.tmp)] <- data.to.organize[data.tmp,]
    }
  }

  return(out)
}



################################################
## functions to transform data into list form ##
################################################
create.list.form <- function(x,n,m){
  subjectID <- x[,1]
  data.to.organize <- data.frame(x[,-1])
  lb <- ncol(data.to.organize)
  uniqueID <- unique(subjectID)

  if(length(uniqueID)!=n){
    stop("The number of unique subject ID's in the data does not match the sample size n.")
  }

  out <- vector("list",n)

  for(ll in 1:length(uniqueID)){
    data.tmp <- which(x[,1]==uniqueID[ll])
    out[[ll]] <- as.matrix(data.to.organize[data.tmp,])
  }

  return(out)
}



###############################################
## Function to  make data suitable for glme4 ##
###############################################

make.data.set <- function(n,m,p,y.data,x.data,z.data){
  ## vector indicating which family each observation belongs
  family <- rep(1:n,times=m)

  data.set.out <- data.frame(cbind(family,y.data[,-1],x.data[,-1],z.data[-1]))
  colnames(data.set.out) <- c("family","Y",paste("X",1:p,sep=""),"Z")
  return(data.set.out)
}






#####################################
## Functions used for random slope ##
#####################################


#This is the function to obtain U for all i, here y is an n*m matrix, and z is an n*m matrix.#
allu <- function(y,z) {
  nofy=nrow(y)
  allu1=vector("list",nofy)
  for (i in 1:nofy) {
    allu1[[i]]=select.vv.random.slope(y[i,],z[i,])
  }
  return (allu1)
}



# function to find all possible u such that u_1+u_2+...+u_(m-1) equals to wi-y[i,1], here, y abd z are m dimentional vecotr for i-th cluster data#
select.vv.random.slope <- function(y,z){
  m=length(y)
  m.tmp <- m-1
  mat.tmp <- matrix(rep(c(0,1),each=m.tmp),ncol=m.tmp,nrow=2,byrow=TRUE)
  vv.combinations <-as.matrix(do.call(`expand.grid`,as.data.frame(mat.tmp)))
  n=nrow(vv.combinations)
  tmp=rep(0,n)

  for (i in 1:n) {
    tmp[i]=sum(vv.combinations[i,]*z[2:m])
  }
  wi=.Call('_rgllvm_getwic',PACKAGE = 'rgllvm',y,z)
  ##  index <- which(tmp<wi & tmp>=(wi-z[1]))
  index <- which(tmp==(wi-y[1]*z[1]))
  uall=vv.combinations[index,]
  uall1=unname(uall)

  return (uall1)
}


#This is the function to obtain suitable u, because previous one has vector inside and needs to be converted to matrix for later computation.#
getnewu <- function(u) {
  n=length(u)

  for (i in 1:n) {
    if (is.vector(u[[i]])) {
      u[[i]]=t(as.matrix(u[[i]]))
    }
  }
  return (u)
}








