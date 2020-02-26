
###############################
## Function to generate data ##
###############################

## n : sample size
## m : vector of number of people in family
## p : dimension(beta)
## fr : density f_R(r)
## fzr : density f_{Z|R}(Z|r)


gendata <- function(beta,n,m=rep(2,n),p,
                    fr,
                    fzr,
                    fr.form = c("normal","gamma","uniform","t",
                                "mixture.normal","mixture.norm.gamma")[1],
                    fzr.form = c("normal","normalR","uniform")[3],
                    center.z=FALSE,ma.orig=FALSE,
                    ytmp=NULL,ztmp=NULL,real.data=FALSE){
  r <- rep(0,n)
  z <- array(NA,dim=c(n,max(m),p),
             dimnames=list(1:n,paste("m",1:max(m),sep=""),paste("p",1:p,sep="")))
  y <- array(NA,dim=c(n,max(m)),
             dimnames=list(1:n,paste("m",1:max(m),sep="")))
  w <- array(NA,dim=c(n,max(m)),
             dimnames=list(1:n,paste("m",1:max(m),sep="")))


  if(real.data==FALSE){
    for(i in 1:n){
      ## generate deviate for random effect
      r[i] <- fr(1)
    }

    ## center the random effects
    if(fzr.form=="depcor0"){
      ## make R such that E(R)=0 and E(R^3)=0
      ## We do this by finding a,b in aR+b
      ##    such that E(aR+b)=0 and E((aR+b)^3)=0.

      get.root.function <- function(R){
        function(theta){
          a1=theta[1]
          b1=theta[2]
          c1=theta[3]

          Rnew <- a1*R + b1
          Xnew <- Rnew^2 + c1

          c(F1=mean(Xnew * Rnew),
            F2=mean(Xnew)-1.5,
            F3=mean(Rnew))
        }
      }

      root.function <- get.root.function(r)
      myroot <- multiroot(root.function,start=c(-1,1,2))
      solution <- myroot$root
      r <- (solution[1]*r+solution[2])

    } else {
      ## center the random effects
      r <- r-mean(r)
    }

    for(i in 1:n){
      for(j in 1:m[i]){
        ##print(j)
        ## compute z_{ij}^T \beta + r_i
        if(ma.orig==TRUE){
          z[i,j,j] <-1
        } else {
          z[i,j,] <- fzr(p,r=r[i])
        }
      }
    }

    ## center the covariates?
    if(center.z==TRUE){
      z <- z-mean(z)
    }

    for(i in 1:n){
      for(j in 1:m[i]){
        w[i,j] <- z[i,j,]  %*% beta + r[i]
        prob   <- 1 / ( 1 + exp(-w[i,j]) )
        if(runif(1) < prob){
          y[i,j] <-1
        } else {
          y[i,j] <-0
        }
        ##y[i,j] <- rbinom(1,1,prob)
      }
    }

  } else {
    ## real data
    tmp <- 0
    for(i in 1:n){
      for(j in 1:m[i]){
        tmp <- tmp + 1
        # print(tmp)

        y[i,j] <- ytmp[tmp]
        z[i,j,] <- as.numeric(ztmp[tmp,])
      }
    }

  }
  list(yout=y,zout=z,r=r)
}



## Functions to generate data
#' @import stats
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

#' @import stats
make.normal.fzr <- function(mean=0,sd=1){
  function(nr,r){
    stats::rnorm(nr,mean=mean,sd=sd)
  }
}

#' @import stats
make.normal.fzr.R <- function(sd=1){
  function(nr,r){
    stats::rnorm(nr,mean=r,sd=sd)
  }
}

#' @import stats
make.normal.fzr.R2 <- function(mean=1){
  function(nr,r){
    stats::rnorm(nr,mean=mean,sd=abs(r))
  }
}

#' @import stats
make.normal.fzr.R <- function(sd=1){
  function(nr,r){
    stats::rnorm(nr,mean=r,sd=sd)
  }
}

#' @import stats
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

#' @import stats
make.uniform <- function(min=-2,max=2){
  function(nr){
    stats::runif(nr,min=min,max=max)
  }
}

#' @import stats
make.t <- function(df=3,ncp=-2){
  function(nr){
    stats::rt(nr,df=df,ncp=ncp)
  }
}

#' @import stats
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

#' @import stats
make.mixture.norm.gamma <- function(mean=0,sd=1,shape=1.5,scale=2,mix=0.5){
  function(nr){
    if(stats::runif(1) < mix){
      stats::rnorm(nr,mean=mean,sd=sd)
    } else {
      stats::rgamma(nr,shape=shape,scale=scale)
    }
  }
}
