
###############################
## Function to generate data ##
###############################

## n : sample size
## m : vector of number of people in family
## p : dimension(beta)
## fr : density f_R(r)
## fzr : density f_{Z|R}(Z|r)

#' Generate repeated measures data
#'
#' @param beta0 true fixed effects parameter.
#' @param n scalar denoting the number of sujects in the study.
#' @param m n-dimensional vector containing the number of repeated measures per subject.
#' @param fr.form character indicating the distribution of the latent variable. Options are "normal", "gamma", "uniform",
#' "t", "mixture.normal" and "mixture.norm.gamma" for a mixture of normal and gamma random variables.
#' @param fxr.form character indicating the distribution of the fixed covariates. Options are "uniform",
#' "uniform-plus-R2" where the covariate is uniformly distributed centered at the square of the latent variable r,
#' "normal","normal-mean-R" where the covariate has a normal distribution with the mean being the latent variable r,
#' "normal-sd-R" where the covariate has a normal distribution with the standard deviation being |r|.
#' @param latent.variable.type character indicating the type of latent variable generated. Options are "intercept" and "slope".
#' Default is "intercept".
#' @param center.x logical indicator. If TRUE, the fixed covariates are centered at zero.
#' @import stats
#' @export
gendata <- function(beta0=c(0.5,0,-0.5),
                    n=300,
                    m=rep(3,n),
                    fr.form = c("normal","gamma","uniform","t",
                                "mixture.normal","mixture.norm.gamma")[1],
                    fxr.form = c("uniform","uniform-plus-R2",
                                 "normal","normal-mean-R","normal-sd-R")[3],
                    latent.variable.type=c("intercept","slope")[1],
                    center.x=FALSE){
  p=length(beta0)

  r <- rep(0,n)
  x <- array(NA,dim=c(n,max(m),p),
             dimnames=list(1:n,paste("m",1:max(m),sep=""),paste("p",1:p,sep="")))
  y <- array(NA,dim=c(n,max(m)),
             dimnames=list(1:n,paste("m",1:max(m),sep="")))
  w <- array(NA,dim=c(n,max(m)),
             dimnames=list(1:n,paste("m",1:max(m),sep="")))
  z <- y
  ############################################
  ## Functions to generate random variables ##
  ############################################

  make.uniform.fzr <- function(min=-1,max=1){
    function(nr,r){
      stats::runif(nr,min=min,max=max)
    }
  }

  make.uniform.fzr.R <- function(min=1,max=2){
    ## Z= a R^2, a\neq 0.
    function(nr,r){
      runif(nr,min=min,max=max) + r^2
    }
  }

  make.normal.fzr <- function(mean=0,sd=1){
    function(nr,r){
      stats::rnorm(nr,mean=mean,sd=sd)
    }
  }

  make.normal.fzr.R <- function(sd=1){
    function(nr,r){
      stats::rnorm(nr,mean=r,sd=sd)
    }
  }

  make.normal.fzr.R2 <- function(mean=1){
    function(nr,r){
      stats::rnorm(nr,mean=mean,sd=abs(r))
    }
  }


  make.normal <- function(mean=0,sd=1){
    function(nr){
      stats::rnorm(nr,mean=mean,sd=sd)
    }
  }

  #' @import stats
  make.gamma <- function(shape=1.5,scale=2){
    function(nr){
      stats::rgamma(nr,shape=shape,scale=scale)
    }
  }

  make.uniform <- function(min=-2,max=2){
    function(nr){
      stats::runif(nr,min=min,max=max)
    }
  }

  make.t <- function(df=3,ncp=-2){
    function(nr){
      stats::rt(nr,df=df,ncp=ncp)
    }
  }

  make.cauchy <- function(location=0,scale=1){
    function(nr){
      stats::rcauchy(nr,location=location,scale=scale)
    }
  }

  #' @import stats
  make.mixture.normals <- function(mean1=0,sd1=1,mean2=1,sd2=2,mix=0.5){
    function(nr){
      if(runif(1) < mix){
        stats::rnorm(nr,mean=mean1,sd=sd1)
      } else {
        stats::rnorm(nr,mean=mean2,sd=sd2)
      }
    }
  }

  make.mixture.norm.gamma <- function(mean=0,sd=1,shape=1.5,scale=2,mix=0.5){
    function(nr){
      if(stats::runif(1) < mix){
        stats::rnorm(nr,mean=mean,sd=sd)
      } else {
        stats::rgamma(nr,shape=shape,scale=scale)
      }
    }
  }

  if(fr.form=="normal"){
    fr <-  make.normal(mean=0,sd=1)
  } else if(fr.form=="gamma"){
    fr <-  make.gamma(shape=1,scale=2)
  } else if(fr.form=="uniform"){
    fr <-  make.uniform(min=-1,max=1)
  } else if(fr.form=="t"){
    fr <- make.t(df=3,ncp=0)
  } else if(fr.form=="mixture.normal"){
    fr <- make.mixture.normals(mean1=3,sd1=1,mean2=-3,sd2=0.25,mix=0.9)
  } else if(fr.form=="mixture.norm.gamma"){
    fr <- make.mixture.norm.gamma(mean=0,sd=sqrt(1),
                                  shape=1.5,scale=0.25,mix=0.1)
  }

  if(fxr.form=="uniform"){
    fxr <-  make.uniform.fzr(min=-1,max=1)
  } else if(fxr.form=="uniform-plus-R2"){
    # X = a R^2, a \neq 0
    fxr <-  make.uniform.fzr.R(min=1,max=2)
  } else if(fxr.form=="normal"){
    fxr <-   make.normal.fzr(mean=0,sd=1)
  } else if(fxr.form=="normal-mean-R"){
    fxr <- make.normal.fzr.R(sd=1)
  } else if(fxr.form=="normal-sd-R"){
    fxr <-   make.normal.fzr.R2(mean=1)
  }

  for(i in 1:n){
    ## generate deviate for random effect
    r[i] <- fr(1)
  }

  ## center the random effects
  r <- r-mean(r)

  for(i in 1:n){
    for(j in 1:m[i]){
      ##print(j)
      ## compute z_{ij}^T \beta + r_i
      x[i,j,] <- fxr(p,r=r[i])

      if(latent.variable.type=="intercept"){
        z[i,j] <- 1
      }
    }
  }

  ## center the covariates?
  if(center.x==TRUE){
    x <- x-mean(x)
  }

  for(i in 1:n){
    for(j in 1:m[i]){
      w[i,j] <- x[i,j,]  %*% beta0 + z[i,j]* r[i]
      prob   <- 1 / ( 1 + exp(-w[i,j]) )
      if(stats::runif(1) < prob){
        y[i,j] <-1
      } else {
        y[i,j] <-0
      }
      ##y[i,j] <- rbinom(1,1,prob)
    }
  }

  ###########################
  ## Put data in long form ##
  ###########################
  family <- rep(1:n,times=m)

  ## make  y into a vector of values
  new.y <- as.vector(t(y))
  new.y <- as.vector(na.omit(new.y))

  ## make z into a vector of values
  new.z <- as.vector(t(z))
  new.z <- as.vector(na.omit(new.z))

  ## make x into a matrix of values
  x.values <- NULL
  for (s in 1:p){
    new.x <- x[,,s]
    new.x <- as.vector(t(new.x))
    new.x <- as.vector(na.omit(new.x))
    x.values <- cbind(x.values,new.x)
  }

  data.set.out <- data.frame(cbind(family,new.y,x.values,new.z))
  colnames(data.set.out) <- c("family","Y",paste("X",1:p,sep=""),"Z")

  return(data.set.out)
}


