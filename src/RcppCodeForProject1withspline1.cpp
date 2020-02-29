# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <math.h>
#include <RcppArmadilloExtensions/sample.h>
#include <vector>

using namespace Rcpp;

//This cpp function,sumofrowbycol, is used to obtain the dot product sum of a row in a matrix and a vector.
//The input x is an m*p matrix, the input beta is a p length vector.
//The output temp is an m*1 matrix, each element is the corresponding dot product.

// [[Rcpp::export]]
arma::mat sumofrowbycol(const arma::mat & x, const arma::vec & beta){
  int m = x.n_rows;
  arma::mat temp(m,1);
  temp.fill(0.0);

  for (int j=0;j<m;j++) {
    temp.submat(j,0,j,0) = x.row(j)*beta;
  }
  return (temp);
}



// Obtainning the score function of i-th cluster
// Here, I omit the i subscript, indicating it is the i-th cluster data we are working on. y is m dimentional, x is (m,p) dimentional, z is m dimentional, mu,sigma,t is one
// dimentional real number,beta is p dimentional,kbeta is the parameter indicating which beta
// we are taking derivate with respect to
// In this program, noticing that I use kbeta-1 instead of kbeta, because in Cpp, the index starts from 0, while in R, the index starts from 1.

//This function, scoref1c,is used to calculate the j-th part of the derived score function.


// [[Rcpp::export]]
List scoref1c(const arma::vec & y, const arma::mat & x, const arma::vec & z, const double & mu, const double & sigma, const double & t, const arma::vec & beta,int kbeta) {
  int m=y.size();
  arma::mat sum1(m,1);
  sum1.fill(0.0);
  arma::mat product1(m,1);
  product1.fill(0.0);


  for (int j=0;j<m;j++) {
    sum1.submat(j,0,j,0)=y(j)*x(j,kbeta-1)-(exp(sumofrowbycol(x,beta)(j,0)+(sqrt(2)*sigma*t+mu)*z(j)))*x(j,kbeta-1)/(1+exp(sumofrowbycol(x,beta)(j,0)+(sqrt(2)*sigma*t+mu)*z(j)));
    product1.submat(j,0,j,0)=exp(y(j)*(sumofrowbycol(x,beta)(j,0)+(sqrt(2)*sigma*t+mu)*z(j))-log(1+exp(sumofrowbycol(x,beta)(j,0)+(sqrt(2)*sigma*t+mu)*z(j))));
  }
  return Rcpp::List::create(
    Rcpp::Named("sum part for each j in score function.") = sum1,
    Rcpp::Named("product part for each j in score function.") = product1
  ) ;
}


//This scoref22 function is used to extract sum1 and product1 in scoref1c function result, and produce the numerator and denominator part for later calculating score function.
//In this function, y is m dimentional vector, x is (m,p) dimentional matrix, z is m dimentional vector,mu,sigma,t is one dimentional real number. beta is p dimentional, kbeta.
//is the parameter indicating which beta we are taking derivative with respect to.
// [[Rcpp::export]]
List scoref2c(const arma::vec & y, const arma::mat & x, const arma::vec & z, const double & mu, const double & sigma, const double & t, const arma::vec & beta,int kbeta) {
  arma::mat sum1=scoref1c(y,x,z,mu,sigma,t,beta,kbeta)(0);
  arma::mat product1=scoref1c(y,x,z,mu,sigma,t,beta,kbeta)(1);
  double sum2=accu(sum1);
  int m=y.size();
  double product2=1;
  double numerator;
  double denominator;

  for (int i=0;i<m;i++) {
    product2=product2*product1(i,0);
  }

  numerator=sum2*product2;
  denominator=product2;

  return Rcpp::List::create(
    Rcpp::Named("numerator.") = numerator,
    Rcpp::Named("denominator.")=denominator
  ) ;
}

//This function, scoreffinalc, is used to calculate the score function, we use Gauss Hermite Quadrature to calculate the numerator and denominator separately, then taking the quotient. t and w are from GHPW. For example, t=unlist(GHPW(30)[[1]]) abd w=unlist(GHPW(30)[[2]]).
// The input y should m dimensional vector, x should (m,p) dimensional matrix, z should be m dimensional vector, mu,sigma is one dimensional real number,beta should be p dimensional vector,kbeta is the parameter indicating which beta we are taking derivative with respect to.


// [[Rcpp::export]]
double scoreffinalc(const arma::vec & y, const arma::mat & x, const arma::vec & z, const double & mu, const double & sigma, const arma::vec & t, const arma::vec & w, const arma::vec & beta,int kbeta) {
  int len=t.size();
  arma::mat numerator2(len,1);
  arma::mat denominator2(len,1);
  numerator2.fill(0.0);
  denominator2.fill(0.0);

  for (int l=0;l<len;l++) {
    double tempw=w(l);
    double tempt=t(l);

    arma::vec tempscoref2c1=scoref2c(y,x,z,mu,sigma,tempt,beta,kbeta)(0);
    arma::vec tempscoref2c2=scoref2c(y,x,z,mu,sigma,tempt,beta,kbeta)(1);

    double tempscoref2c1d=accu(tempscoref2c1);
    double tempscoref2c2d=accu(tempscoref2c2);

    numerator2(l)=tempw*tempscoref2c1d;
    denominator2(l)=tempw*tempscoref2c2d;

  }
  double numerator2sum=accu(numerator2);
  double denominator2sum=accu(denominator2);
  double scorefunctionfinal=numerator2sum/denominator2sum;

  return (scorefunctionfinal);
}

//This function, getwic, is used to get W_i in equation (5), since now our z_ij is one dimention, wi is just a real number.

// [[Rcpp::export]]
double getwic(const arma::vec & y, const arma::vec & z) {
  arma::vec temp1=y.t()*z;
  double wi=accu(temp1);

  return (wi);
}

//This function,getzilrc,is used to get Z_iL Z_iR, Z_iL inverse in equation (5), here Z_iL and Z_iL inverse is a 1*1 matrix as q now is 1. Z_iR is 1*(m-1) matrix.
//The input z should be a m dimentional vector.


// [[Rcpp::export]]

List getzilrc(const arma::vec & z) {
  int dimz=z.size();
  arma::vec zil=z.subvec(0,0);
  arma::vec zilinv=1/zil;
  arma::vec zir=z.subvec(1,dimz-1);
  arma::mat minvl(dimz,1);
  arma::mat minvlbottom(dimz-1,1);
  minvlbottom.fill(0.0);
  minvl.fill(0.0);
  minvl.submat(0,0,0,0)=zilinv;
  minvl.submat(1,0,dimz-1,0)=minvlbottom;
  arma::mat minvr(dimz,dimz-1);
  minvr.fill(0.0);
  minvr.submat(0,0,0,dimz-2)=-zilinv*zir.t();
  arma::mat eye1(dimz-1,dimz-1);
  eye1.eye();
  minvr.submat(1,0,dimz-1,dimz-2)=eye1;
  arma::mat minv(dimz,dimz);
  minv.fill(0.0);
  minv.submat(0,0,dimz-1,0)=minvl;
  minv.submat(0,1,dimz-1,dimz-1)=minvr;
  arma::mat zirtrans(1,dimz-1);
  zirtrans=zir.t();

  return Rcpp::List::create(
    Rcpp::Named("z_iL.") = zil,
    Rcpp::Named("z_iL inverse.") = zilinv,
    Rcpp::Named("z_iR.") = zirtrans,
    Rcpp::Named("M inverse.")=minv
  ) ;
}


//This function, xbetamatrixc, is used to obtain following matrices: xbmatrix1 is (x_1t*beta,x_2t*beta,...,x_qt*beta),xbmatrix2 is (x_(q+1)t*beta,....,x_mt*beta).
//The input x is a m*p matrix and beta is a p dimentional vector.

// [[Rcpp::export]]
List xbetamatrixc(const arma::mat & x, const arma::vec & beta) {
  int m=x.n_rows;
  arma::mat xbmatrix1(1,1);
  xbmatrix1.fill(0.0);
  arma::mat xbmatrix2(m-1,1);
  xbmatrix2.fill(0.0);
  xbmatrix1.submat(0,0,0,0)=sumofrowbycol(x,beta)(0,0);

  for (int j=1;j<m;j++) {
    xbmatrix2.submat(j-1,0,j-1,0)=sumofrowbycol(x,beta)(j,0);

  }

  return Rcpp::List::create(
    Rcpp::Named("xbmatrix1.") = xbmatrix1,
    Rcpp::Named("xbmatrix2.") = xbmatrix2.t()
  );
}

//This function,exppartin5c,is used to calculate the exponential part in (5).
//The input x is an m*p matrix, z is an m dimentional vector, u is an m-1 vector, indicating each choice of u.
//in (5) for i-th cluster data.

// [[Rcpp::export]]
arma::mat exppartin5c(const arma::mat & x, const arma::vec & z, const arma::vec & u, const arma::vec & beta) {
  arma::mat temp1=xbetamatrixc(x,beta)(0);
  arma::mat temp2=xbetamatrixc(x,beta)(1);
  arma::mat temp3=getzilrc(z)(1);
  arma::mat temp4=getzilrc(z)(2);

  arma::mat exppin5=exp(temp1*(-temp3*temp4*u)+temp2*u);

  return (exppin5);
}

//This function, yin5c, is used to get the first argument formula in (5) score function, here we require y and z to be just m dimentional vector, u is a 1*(m-1) dimentional matrix. This function will return a m length vector#
//later to be used in score function in (5). We require the input u to be uall[[m-1]][o,], where o is the index from 1 to nrow(uall[[m-1]])#

// [[Rcpp::export]]

arma::mat yin5c(const arma::vec & y, const arma::vec & z, const arma::vec & u) {
  int dimz=z.size();
  arma::mat temp1=getzilrc(z)(3);
  double wi=getwic(y,z);
  arma::mat wimat(1,1);
  wimat=wi;
  arma::mat temp2(1,dimz);
  temp2.fill(0.0);
  temp2.submat(0,0,0,0)=wimat;
  arma::mat temp3=u;
  temp2.submat(0,1,0,dimz-1)=temp3.t();
  arma::mat temp4=getzilrc(z)(3);
  arma::mat yinfive=temp4*temp2.t();

  return (yinfive.t());
}

//This function, formula5denoc, is used to calculate the denominator in (5).Here, u should be uall[[m-1]].
//Thus u here is a matrix, whose row number is number of choice of u, and column number is m-1.

// [[Rcpp::export]]

double formula5denoc( const arma::mat & x, const arma::vec & z, const arma::mat & u, const arma::vec & beta) {
  int nofu=u.n_rows;
  int mofu=u.n_cols;
  arma::mat denotemp(1,nofu);
  denotemp.fill(0.0);
  double denosum=0;

  for (int l=0;l<nofu;l++) {
    arma::mat utemp=u.submat(l,0,l,mofu-1);
    utemp.reshape(mofu,1);
    double temp2=exppartin5c(x,z,utemp,beta)(0,0);
    denotemp(0,l)=temp2;
  }
  denosum=accu(denotemp);

  return (denosum);
}


//This function, formula5numerc, is used to calculate the numeritor in (5). Here, u should be uall[[m-1]].
//The input y is an m dimentional vector, x is an m*p matrix, z is an m dimentional vector.
//u here is a matrix, whose row number is number of choice of u, and column number is m-1.
//mu and sigma is an one dimentional real number.
//t and w are from GHPW. For example, t=unlist(GHPW(30)[[1]]) abd w=unlist(GHPW(30)[[2]]).
//beta is a p dimentional vector, kbeta is the integer parameter indicating which beta we are taking derivate with respect to.

// [[Rcpp::export]]

double formula5numerc(const arma::vec & y, const arma::mat & x, const arma::vec & z, const arma::mat & u, const double & mu, const double & sigma,  const arma::vec & t, const arma::vec & w, const arma::vec & beta, int kbeta) {
  int nofu=u.n_rows;
  int mofu=u.n_cols;
  arma::mat numertemp(1,nofu);
  numertemp.fill(0.0);
  double numersum=0;

  for (int l=0;l<nofu;l++) {
    arma::mat utemp=u.submat(l,0,l,mofu-1);
    arma::mat mwu=yin5c(y,z,utemp.t());
    numertemp.submat(0,l,0,l)=scoreffinalc(mwu.t(),x,z,mu,sigma,t,w,beta,kbeta)*exppartin5c(x,z,utemp.t(),beta);
  }
  numersum=accu(numertemp);

  return (numersum);
}

//This function, formula5forallkc, is used to calculate the formula (5) for all kbeta.

// [[Rcpp::export]]

arma::vec formula5forallkc(const arma::vec & y, const arma::mat & x, const arma::vec & z, const arma::mat & u, const double & mu, const double & sigma, const arma::vec & t, const arma::vec & w, const arma::vec & beta) {
  double deno1=formula5denoc(x,z,u,beta);
  int k=beta.size();
  arma::vec formula5fork(k);
  formula5fork.fill(0.0);
  arma::vec formula5forallk(k);
  formula5forallk.fill(0.0);

  for (int i=0;i<k;i++) {
    formula5fork(i)=formula5numerc(y,x,z,u,mu,sigma,t,w,beta,i+1);
  }
  formula5forallk=formula5fork/deno1;

  return (formula5forallk);
}

//This function,efficientscorefc, is used to obtain the locally efficient score function for i-th cluster data.
//The input y should m dimensional vector, x should (m,p) dimensional matrix, z should be m dimensional vector, mu,sigma is one dimensional real number,beta should be p dimensional vector.
//t and w are from GHPW. For example, t=unlist(GHPW(30)[[1]]) abd w=unlist(GHPW(30)[[2]]).

// [[Rcpp::export]]

arma::vec efficientscorefc(const arma::vec & y, const arma::mat & x, const arma::vec & z, const arma::mat & u, const double & mu, const double & sigma, const arma::vec & t, const arma::vec & w, const arma::vec & beta) {
  int q=beta.size();
  arma::vec scoref(q);
  scoref.fill(0.0);
  arma::vec efficientf(q);
  efficientf.fill(0.0);

  for (int i=0;i<q;i++) {
    scoref(i)=scoreffinalc(y,x,z,mu,sigma,t,w,beta,i+1);
  }
  efficientf=scoref-formula5forallkc(y,x,z,u,mu,sigma,t,w,beta);

  return (efficientf);
}


//This function,esteqc, is used to obtain the estimating equation.
//u should be a List, where number of matrix inside is row number of y.
//The input y is a n*m matrix, x is a List of n element, where each element is an m*p matrix.
//

// [[Rcpp::export]]

List esteqc(const arma::mat & y, const List & x, const arma::mat & z, const List & u, const double & mu, const double & sigma,  const arma::vec & t, const arma::vec & w, const arma::vec & beta) {
  int k=y.n_rows;
  List esteqsum(k);
  esteqsum.fill(0.0);

  for (int i=0;i<k;i++) {
    arma::mat temp1=u(i);
    arma::vec temp2=(y.row(i)).t();
    arma::vec temp3=(z.row(i)).t();
    arma::mat temp4=x(i);

    esteqsum(i)=efficientscorefc(temp2,temp4,temp3,temp1,mu,sigma,t,w,beta);
  }

  return (esteqsum);
}

//This function is used to extrac element of list from esteq to formulate f_beta_i for system of non linear equations solver.

// [[Rcpp::export]]

arma::vec esteqfc(const arma::mat & y, const List & x, const arma::mat & z, const List & u, const double & mu, const double & sigma,  const arma::vec & t, const arma::vec & w, const arma::vec & beta) {
  List esteqsum=esteqc(y,x,z,u,mu,sigma,t,w,beta);
  int k=y.n_rows;
  int l=beta.size();
  arma::vec sum1(l);
  sum1.fill(0.0);

  for (int j=0;j<l;j++) {
    for (int i=0;i<k;i++) {
      arma::vec temp1=esteqsum(i);
      sum1(j)=sum1(j)+temp1(j);
    }
  }

  return (sum1);
}

//This function, PCLijklc, is used to prepared the integral function suitable for the Gauss Hermite Quadrature Method.

// [[Rcpp::export]]

double PCLijklc(const arma::mat & xij, const arma::mat & xkl, const arma::vec & beta, const double & sigmaq, const double & t) {
  int p=beta.size();
  double temp1=0;


  for (int i=0;i<p;i++) {
    double temp2=xkl(0,i)*beta(i)-xij(0,i)*beta(i);
    temp1=temp1+temp2;
  }

  double exppart=exp(temp1+sqrt(2)*sigmaq*t);
  double intpart=(1/sqrt(3.14159265358979323846))*(1/(1+exppart));

  return (intpart);
}

//This function, PCLijklalltc, is used to calculate the integral for each X_ij,X_kl, using the Gauss Hermite Quadrature Method.
//The input xij,xkl is vector of length p,where p is the length of beta, sigmaq is a real number.
//t and w are from GHPW. For example, t=unlist(GHPW(30)[[1]]) abd w=unlist(GHPW(30)[[2]]).

// [[Rcpp::export]]

double PCLijklalltc(const arma::mat & xij, const arma::mat & xkl, const arma::vec & beta, const double & sigmaq, const arma::vec & t, const arma::vec & w) {
  int lt=t.size();
  double temp1=0;


  for (int i=0;i<lt;i++) {
    double temp2=w(i)*PCLijklc(xij,xkl,beta,sigmaq,t(i));
    temp1=temp1+temp2;
  }

  return (temp1);
}

//This function, PCLallc, is used to calculate the log of the pseudo-conditional likelihood.
//The input y is the y matirx of size n*m, the input x is a list of size n, where each element in the list is a m*p matrix.
//The input beta is a length p vector. sigmaq is a real number, whose square is the variance of distribution of Q_(ij),(kl).
//t and w are from GHPW. For example, t=unlist(GHPW(30)[[1]]) abd w=unlist(GHPW(30)[[2]]).

// [[Rcpp::export]]

double PCLallc(const arma::mat & y,const List & x, const arma::vec & beta, const double & sigmaq, const arma::vec & t,const arma::vec & w) {
  int n=y.n_rows;
  int m=y.n_cols;
  double temp4=0;

  for (int i=0;i<n;i++) {
    for (int j=0;j<m;j++) {
      for (int k=i+1;k<n;k++) {
        for (int l=0;l<m;l++) {
          double temp1=y(i,j);
          double temp2=y(k,l);

          if (temp1>0.5 && temp2<0.5) {
            arma::mat tempxi=x(i);
            arma::mat tempxij=tempxi.row(j);
            arma::mat tempxk=x(k);
            arma::mat tempxkl=tempxk.row(l);

            double temppcl1=log(PCLijklalltc(tempxij,tempxkl,beta,sigmaq,t,w));

            double temp3=temp1*(1-temp2)*temppcl1;
            temp4=temp4+temp3;
            continue;

          }

          if (temp1<0.5 && temp2>0.5) {
            arma::mat tempxi=x(i);
            arma::mat tempxij=tempxi.row(j);
            arma::mat tempxk=x(k);
            arma::mat tempxkl=tempxk.row(l);

            double temppcl2=log(PCLijklalltc(tempxkl,tempxij,beta,sigmaq,t,w));

            double temp3=temp2*(1-temp1)*temppcl2;
            temp4=temp4+temp3;
            continue;
          }
        }
      }
    }
  }

  return (2*temp4/(n*(n-1)));
}

//This function, PCLaikc, is used to calculate the a(D_i,D_k,theta) part of the log of the pseudo-conditional likelihood in the draft.
//The input i indicates the i-th data point and k indicates the k-th data point.
//The input y is the y matirx of size n*m, the input x is a list of size n, where each element in the list is a m*p matrix.
//The input beta is a length p vector. sigmaq is a real number, whose square is the variance of distribution of Q_(ij),(kl).
//t and w are from GHPW. For example, t=unlist(GHPW(30)[[1]]) abd w=unlist(GHPW(30)[[2]]).

// [[Rcpp::export]]

double PCLaikc(int i, int k, const arma::mat & y,const List & x,const arma::vec & beta, const double & sigmaq,const arma::vec & t,const arma::vec & w) {
  int m=y.n_cols;
  double temp4=0;

  for (int j=0;j<m;j++) {
    for (int l=0;l<m;l++) {
      double temp1=y(i,j);
      double temp2=y(k,l);
      arma::mat tempxi=x(i);
      arma::mat tempxij=tempxi.row(j);
      arma::mat tempxk=x(k);
      arma::mat tempxkl=tempxk.row(l);

      double temp3=temp1*(1-temp2)*log(PCLijklalltc(tempxij,tempxkl,beta,sigmaq,t,w))+temp2*(1-temp1)*log(PCLijklalltc(tempxkl,tempxij,beta,sigmaq,t,w));
      temp4=temp4+temp3;
    }
  }

  return (temp4);
}

//This function,firderPCLaikc, is used to calculate the first derivative of a(D_i,D_k,theta) part of the log of the pseudo-conditional likelihood in the draft.
//The input i indicates the i-th data point and k indicates the k-th data point.
//The input y is the y matirx of size n*m, the input x is a list of size n, where each element in the list is a m*p matrix.
//The input beta is a length p vector. sigmaq is a real number, whose square is the variance of distribution of Q_(ij),(kl).
//t and w are from GHPW. For example, t=unlist(GHPW(30)[[1]]) abd w=unlist(GHPW(30)[[2]]).


// [[Rcpp::export]]

arma::vec firderPCLaikc(int i, int k, const arma::mat & y,const List & x,const arma::vec & beta, const double & sigmaq,const arma::vec & t,const arma::vec & w) {
  int p=beta.size();

  double delta0=sigmaq*0.001;

  arma::vec betader(p);
  betader.fill(0.0);

  for (int o=0;o<p;o++) {
    arma::vec deltatemp(p);
    deltatemp.fill(0.0);
    deltatemp(o)=beta(o)*0.001;

    betader(o)=(PCLaikc(i,k,y,x,beta+deltatemp,sigmaq,t,w)-PCLaikc(i,k,y,x,beta-deltatemp,sigmaq,t,w))/(2*deltatemp(o));

  }

  double sigmaqder=(PCLaikc(i,k,y,x,beta,sigmaq+delta0,t,w)-PCLaikc(i,k,y,x,beta,sigmaq-delta0,t,w))/(2*delta0);

  arma::vec firstderPCLaik(p+1);
  firstderPCLaik.subvec(0,0)=sigmaqder;
  firstderPCLaik.subvec(1,p)=betader;

  return (firstderPCLaik);
}



//This function,secderPCLaikc, is used to calculate the second derivative of a(D_i,D_k,theta) part of the log of the pseudo-conditional likelihood in the draft.
//The input i indicates the i-th data point and k indicates the k-th data point.
//The input y is the y matirx of size n*m, the input x is a list of size n, where each element in the list is a m*p matrix.
//The input beta is a length p vector. sigmaq is a real number, whose square is the variance of distribution of Q_(ij),(kl).
//t and w are from GHPW. For example, t=unlist(GHPW(30)[[1]]) abd w=unlist(GHPW(30)[[2]]).


// [[Rcpp::export]]

arma::mat secderPCLaikc(int i, int k, const arma::mat & y,const List & x,const arma::vec & beta, const double & sigmaq,const arma::vec & t,const arma::vec & w) {
  int p=beta.size();

  double delta2=sigmaq*0.001;

  arma::mat secderrestcol(p+1,p+1);
  secderrestcol.fill(0.0);

  arma::vec secdercol1=(firderPCLaikc(i,k,y,x,beta,sigmaq+delta2,t,w)-firderPCLaikc(i,k,y,x,beta,sigmaq-delta2,t,w))/(2*delta2);
  secderrestcol.submat(0,0,p,0)=secdercol1;

  for (int s=0;s<p;s++) {
    arma::vec deltatemp(p);
    deltatemp.fill(0.0);

    deltatemp(s)=beta(s)*0.001;

    arma::vec coltemp=(firderPCLaikc(i,k,y,x,beta+deltatemp,sigmaq,t,w)-firderPCLaikc(i,k,y,x,beta-deltatemp,sigmaq,t,w))/(2*deltatemp(s));

    secderrestcol.submat(0,s+1,p,s+1)=coltemp;
  }

  return (secderrestcol);

}


//This function,secderPCLallc, is used to calculate the A inverse part consists of the second derivative of a(D_i,D_k,theta) part of the log of the pseudo-conditional likelihood in the draft for the sandwich matrix.
//The input y is the y matirx of size n*m, the input x is a list of size n, where each element in the list is a m*p matrix.
//The input beta is a length p vector. sigmaq is a real number, whose square is the variance of distribution of Q_(ij),(kl).
//t and w are from GHPW. For example, t=unlist(GHPW(30)[[1]]) abd w=unlist(GHPW(30)[[2]]).

// [[Rcpp::export]]

arma::mat secderPCLaikallc(const arma::mat & y,const List & x,const arma::vec & beta, const double & sigmaq,const arma::vec & t,const arma::vec & w) {
  int n=y.n_rows;
  int p=beta.size();
  arma::mat temp1(p+1,p+1);
  temp1.fill(0.0);

  for (int i=0;i<n;i++) {
    for (int k=i+1;k<n;k++) {
      arma::mat temp2=secderPCLaikc(i,k,y,x,beta,sigmaq,t,w);
      temp1=temp1+temp2;
    }
  }

  return (inv(temp1));

}

//This function,varofsthetac, is used to calculate the B part of the log of the pseudo-conditional likelihood in the draft for the sandwich matrix.
//The input y is the y matirx of size n*m, the input x is a list of size n, where each element in the list is a m*p matrix.
//The input beta is a length p vector. sigmaq is a real number, whose square is the variance of distribution of Q_(ij),(kl).
//t and w are from GHPW. For example, t=unlist(GHPW(30)[[1]]) abd w=unlist(GHPW(30)[[2]]).


// [[Rcpp::export]]
arma::mat varofsthetac(const arma::mat & y,const List & x,const arma::vec & beta, const double & sigmaq,const arma::vec & t,const arma::vec & w) {
  int n=y.n_rows;
  int p=beta.size();
  arma::mat temp1(p+1,p+1);
  temp1.fill(0.0);

  arma::cube allfirderPCLaik(p+1,n,n);
  allfirderPCLaik.fill(0.0);

  for (int i=0;i<n;i++) {
    for (int k=i+1;k<n;k++) {
      allfirderPCLaik.subcube(0,k,i,p,k,i)=firderPCLaikc(i,k,y,x,beta,sigmaq,t,w);
    }
  }

  for (int i=0;i<n;i++) {
    for (int k=i+1;k<n;k++) {
      arma::vec temp2=allfirderPCLaik.subcube(0,k,i,p,k,i);
      temp1=temp1+temp2*trans(temp2);
    }

  }


  for (int i=0;i<n;i++) {
    for (int k=i+1;k<n;k++) {
      for (int kp=i+1;kp<n;kp++) {
        if (k!=kp) {
          arma::vec temp3=allfirderPCLaik.subcube(0,k,i,p,k,i);
          arma::vec temp4=allfirderPCLaik.subcube(0,kp,i,p,kp,i);
          temp1=temp1+temp3*trans(temp4);
        }
      }
    }
  }

  for (int k=1;k<n;k++) {
    for (int i=0;i<k-1;i++) {
      for (int ip=0;ip<k-1;ip++) {
        if (i!=ip) {
          arma::vec temp5=allfirderPCLaik.subcube(0,k,i,p,k,i);
          arma::vec temp6=allfirderPCLaik.subcube(0,k,ip,p,k,ip);
          temp1=temp1+temp5*trans(temp6);
        }

      }
    }

  }

  return (temp1);
}

//This function,varofsthetac, is used to calculate the sandwich matrix of the log of the pseudo-conditional likelihood in the draft.
//The input y is the y matirx of size n*m, the input x is a list of size n, where each element in the list is a m*p matrix.
//The input beta is a length p vector. sigmaq is a real number, whose square is the variance of distribution of Q_(ij),(kl).
//t and w are from GHPW. For example, t=unlist(GHPW(30)[[1]]) abd w=unlist(GHPW(30)[[2]]).

// [[Rcpp::export]]

arma::mat PCLsandwich(const arma::mat & y,const List & x,const arma::vec & beta, const double & sigmaq,const arma::vec & t,const arma::vec & w) {
  arma::mat a=secderPCLaikallc(y,x,beta,sigmaq,t,w);
  arma::mat b=varofsthetac(y,x,beta,sigmaq,t,w);

  arma::mat sandwich=a*b*trans(a);

  return (sandwich);
}

//The following is to test if my previous derivative function is correct.

// [[Rcpp::export]]

double test1(arma::vec beta) {
  double temp1=3*beta(0)*beta(0)+4*beta(1)*beta(1)*beta(1)+5*beta(0)*beta(0)*beta(2);

  return (temp1);
}

// [[Rcpp::export]]

arma::vec firdertest1(arma::vec beta) {
  arma::vec temp2(3);
  temp2.fill(0.0);

  for (int i=0;i<3;i++) {
    arma::vec deltatemp(3);
    deltatemp.fill(0.0);
    deltatemp(i)=beta(i)*0.001;

    temp2(i)=(test1(beta+deltatemp)-test1(beta-deltatemp))/(2*deltatemp(i));
  }

  return (temp2);
}

// [[Rcpp::export]]

arma::mat secdertest1(arma::vec beta) {
  arma::mat temp3(3,3);
  temp3.fill(0.0);

  for (int j=0;j<3;j++) {
    arma::vec deltatemp1(3);
    deltatemp1.fill(0.0);
    deltatemp1(j)=beta(j)*0.001;

    arma::vec coltemp=(firdertest1(beta+deltatemp1)-firdertest1(beta-deltatemp1))/(2*deltatemp1(j));

    temp3.submat(0,j,2,j)=coltemp;
  }

  return (temp3);

}

//End of test.


//This function ,pic, is to calculate Pi_i in the derivative of log likelihood of logistic regression wrt
//beta_k(where k is from 0 to K, K is product of number of basis we need for X and Z, plus the number of beta
//we want to estimate, plus 1(intecept), for example, when we have 2 beta to estimate, and number of basis for X and Z
// is 6, K is then (6*6)*6+2+1, where (6*6) stands for combination of basis for x1 and x1, 6 is number of basis for Z,2 is
// the dimension of beta,1 is intercept.)
//The input xi is the i-th row of output of b-spline basis function, where i is from 1 to n*m, where n is number of subject
//m is the number of observation for each subject.
//The input beta is all the parameter we want to estimate, the first component should be intercept alpha in the paper,
//here we use beta_0 to denote alpha, beta_1 to beta_p is the fix effect we want, beta_p+1 to beta_p+(K-1-p) is the coefficient
//for basis we need to estimate but later we will throw them away.


// [[Rcpp::export]]

double pic(const arma::vec & xi, const arma::vec & beta) {
  double sum1=exp(accu(xi%beta));

  double exppartofpii=sum1/(sum1+1);

  return (exppartofpii);
}


//This function, esteqbetaic, is used to calculate the derivative of log likelihood of logistic regression wrt
//beta_k(where k is from 0 to K, K is product of number of basis we need for X and Z, plus the number of beta
//we want to estimate, plus 1(intecept), for example, when we have 2 beta to estimate, and number of basis for X and Z
// is 6, K is then (6*6)*6+2+1, where (6*6) stands for combination of basis for x1 and x1, 6 is number of basis for Z,2 is
// the dimension of beta,1 is intercept.)
//the input l is from 0 to K, K is defined the same as above, to determine which beta we are taking the derivative wrt.
//the input y is a vector of length n*m, where we put the entries of original y matrix of n*m into a vector, where the first m
//component is the first row of original y matrix.
//the input x is a matrix of nm*K, using first row as an example, the first component is always 1, the second component to p+1 component is
//x value corresponding to beta_1 to beta_p we would like to estiamte, p+2 component to p+(K-p-1) is the product of basis function I wrote in
//the draft.
//The input beta is all the parameter we want to estimate, the first component should be intercept alpha in the paper,
//here we use beta_0 to denote alpha, beta_1 to beta_p is the fix effect we want, beta_p+1 to beta_p+(K-1-p) is the coefficient
//for basis we need to estimate but later we will throw them away.

// [[Rcpp::export]]

double esteqbetaic(const int & l, const arma::vec & y, const arma::mat & x, const arma::vec & beta) {
  int n=y.size();
  int m=beta.size();
  arma::vec xl=x.col(l);
  double sum1=0;
  arma::vec xrowi(m);
  xrowi.fill(0.0);

  for (int i=0;i<n;i++) {
    double yi=y(i);
    double xil=xl(i);
    xrowi=trans(x.row(i));

    double pii=pic(xrowi,beta);

    sum1=sum1+yi*xil-pii*xil;

  }

  return (sum1);
}
//This function, esteqallbetaic,is used to calucalte all the derivative of log likelihood of logistic regression wrt.
//The output is a vector of same length as beta, so later we could feed this output to nleqslve function in r to solve this non-linear
//system of equations.
//The input is the same as in function esteqbetaic.


// [[Rcpp::export]]

arma::vec esteqallbetaic(const arma::vec & y, const arma::mat & x, const arma::vec & beta) {
  int s=beta.size();
  arma::vec allesteq(s);
  allesteq.fill(0.0);

  for (int betasub=0;betasub<s;betasub++) {

    allesteq(betasub)=esteqbetaic(betasub,y,x,beta);

  }

  return (allesteq);
}

//This function, sinforcell, is used to calculate the cell of fisher information matrix of spline estimator.
//Note that the return value acutally negative second derivate of log-likelihood function. So we can use it directly
//to form the fisher information matrix, no need to take additional negative.
//The input k is the subscript of the first variable we want to take derivative wrt, kprime is the subscript of the second
// variable we want to take derivative wrt.
//The input matrix x is the transformed matrix of n*m row and p columns from the spline settings.



// [[Rcpp::export]]

double sinforcell(const int & k, const int & kprime, const arma::mat & x, const arma::vec & beta, const arma::vec & pimat1) {
  int n=x.n_rows;
  //int p=beta.size();
  double sum1=0;

  for (int i=0;i<n;i++) {
    double pii=pimat1(i);
    double xik=x(i,k);
    double xikp=x(i,kprime);

    sum1=sum1+xik*pii*(1-pii)*xikp;
  }

  return (sum1);
}

// [[Rcpp::export]]

arma::mat sfisherinforinv(const arma::mat & x, const arma::vec & beta) {
  int p=beta.size();
  arma::mat fisher(p,p);
  fisher.fill(0.0);
  int n=x.n_rows;
  arma::vec xrowi(p);
  xrowi.fill(0.0);
  arma::vec pimat(n);
  pimat.fill(0.0);

  for (int i=0;i<n;i++) {
    arma::vec xrowi=trans(x.row(i));
    double pii=pic(xrowi,beta);
    pimat(i)=pii;
  }

  for (int i=0;i<p;i++) {
    for (int j=0;j<p;j++) {
      fisher(i,j)=sinforcell(i,j,x,beta,pimat);
    }
  }

  return (fisher);
}

