! 5-12-2014: Main functions used in COHORT ext code

subroutine getvv_terms(n,m,maxm)
  implicit none
  integer, intent(in) :: n,maxm
  integer,dimension(n),intent(in) :: m
  integer, dimension(maxm,maxm,int(2.**maxm)) :: index
  integer, dimension(maxm,maxm) :: index_count
  double precision, dimension(maxm,int(2.**maxm),maxm) :: vv_combinations

  call getvv(n,m,maxm,vv_combinations)
  call getindex(maxm,vv_combinations,index,index_count)
  !  do mm=1,maxm
  !     print*,'mm=',mm
  !     do jj=1,int(2.**mm)
  !        print*,'jj=',jj,'vv=',vv_combinations(mm,jj,1:mm)
  !     end do
  !  end do
  !
  !  do mm=1,maxm
  !     print*,'mm=',mm
  !     do kkk=1,mm
  !        print*,'kkk=',kkk,'index=',index(mm,kkk,:),'index_count=',index_count(mm,kkk)
  !     end do
  !  end do
  return
end subroutine getvv_terms

! subroutine to get vv_combinations
subroutine getvv(n,m,maxm,vv_combinations)
  implicit none
  integer, intent (in) :: n,maxm
  integer, dimension(n),intent(in) :: m
  double precision,dimension(maxm,int(2.**maxm),maxm),intent(out) :: vv_combinations
  logical :: one
  integer :: index1,index2,mm,j
  double precision :: tmp
  ! set logical indicator
  do mm=1,maxm
     do j=1,mm
        one=.false.
        index1=1
        index2=int(2.**mm/2.**j)
        do while(index2.le.int(2.**mm))
           if(one) then
              tmp=1.
              one=.false.
           else
              tmp=0.
              one=.true.
           end if
           vv_combinations(mm,index1:index2,j)=tmp
           index1=index2+1
           index2=index2+int(2.**mm/2.**j)
        end do
     end do
  end do
  return
end subroutine getvv

! subroutine to get indices for different w_i values
!   based on vv_combinations
subroutine getindex(maxm,vv_combinations,index,index_count)
  implicit none
  integer, intent (in) :: maxm
  double precision,dimension(maxm,int(2.**maxm),maxm),intent(in) :: vv_combinations
  integer,dimension(maxm,maxm,int(2.**maxm)),intent(out) :: index
  integer,dimension(maxm,maxm),intent(out) :: index_count
  double precision,allocatable,dimension(:) :: mysum
  integer :: k,j,mm
  do mm=1,maxm
     allocate(mysum(int(2.**mm)))
     do k=1,int(2**mm)
        mysum(k)=sum(vv_combinations(mm,k,:))
     end do
     do k=1,mm
        index_count(mm,k)=0
        do j=1,int(2.**mm)
           if(mysum(j).eq.dble(k)) then
              index_count(mm,k)=index_count(mm,k)+1
              index(mm,k,index_count(mm,k))=j
           end if
        end do
     end do
     deallocate(mysum)
  end do
  return
end subroutine getindex


!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! Functions used to compute expectation
!-------------------------------------------------

! subroutine to compute E(V_i|W_i,Z_i;beta)
subroutine expect_v_logistic(lb,n,beta,ynew,z,ev)
  use commondata, only :
       vv_combinations,index,index_count,&
       n1,m1
  implicit none
  integer, intent(in) :: lb, n, maxm
  double precision,dimension(lb),intent(in) :: beta
  double precision,dimension(n,maxm),intent(in) :: ynew
  double precision,dimension(n,maxm,lb),intent(in) :: z
  double precision,dimension(n,maxm),intent(out) :: ev
  double precision,allocatable,dimension(:,:) :: Ainv
  double precision,dimension(lb) :: ztmp
  double precision :: wi,main_tmp,denominator,tmp2
  double precision,allocatable,dimension(:) :: theta,thetaz,Vtmp,vv_new,numerator,tmp1,vv_new2,kappa,kappaz
  integer,dimension(n) :: m
  integer,allocatable,dimension(:) :: other_index
  integer :: i,j,kk,uu,ii,mi,jj
  tol=1e-6
  ev=0.
  m=m1
  do i=1,n1
     !print*,'i=1',i

     mi=m(i)

     if(m1(i).gt.1.) then
        ! allocate dimensions
        allocate(Ainv(m(i),m(i)))
        allocate(theta(m(i)))
        allocate(thetaz(m(i)))
        allocate(Vtmp(m(i)))
        allocate(vv_new(m(i)-1))
        allocate(vv_new2(m(i)))
        allocate(numerator(m(i)-1))
        allocate(tmp1(m(i)))
        allocate(other_index(m(i)-1))
        allocate(kappa(m(i)))
        allocate(kappaz(m(i)))

        ! get Ainv
        call get_Ainv(m(i),Ainv)

        ! compute thetaz=beta^T* z
        thetaz=matmul(z(i,1:m(i),:),beta)
        !print*,'thetaz=',thetaz

        ! theta corresponds to \xi in main text
        theta=thetaz


        ! get wi
        wi = sum(ynew(i,1:m(i)))
        !print*,'wi=',wi

        ! get kappa(j)
        do j=1,m(i)
           ztmp=z(i,j,:)-z(i,1,:)
           kappaz(j)=dot_product(ztmp,beta)
           kappa(j)=kappaz(j)
        end do

        ! compute ev
        if(wi<tol) then
           ev(i,2:m(i))=0.
        else if( (wi.ge.(m(i)-tol)).and.(wi.le.(m(i)+tol))) then
           ev(i,2:m(i))=1.
        else
           !---------------------!
           ! no missing values   !
           !---------------------!
           numerator=0.
           denominator=0.

           do uu=1,index_count(m(i),int(wi))
              !print*,'uu=',uu
              vv_new=vv_combinations(m(i),index(m(i),int(wi),uu),2:m(i))
              !print*,'vv_new=',vv_new
              Vtmp(1)=wi
              Vtmp(2:m(i))=vv_new
              tmp1=matmul(Ainv,Vtmp)
              !print*,'Vtmp=',Vtmp
              !print*,'tmp1(:)=',tmp1(:)
              tmp2=dot_product(theta,tmp1)
              !print*,'theta*tmp1=',tmp2
              main_tmp=exp(tmp2)
              numerator=numerator + vv_new *main_tmp
              denominator=denominator + main_tmp
           end do
           if(abs(denominator).gt.tol) then
              ev(i,2:m(i))=numerator/denominator
           else
              ev(i,2:m(i))=0.
           end if
        end if

        ! de-allocate size of Ainv
        deallocate(Ainv)
        deallocate(theta)
        deallocate(thetaz)
        deallocate(Vtmp)
        deallocate(vv_new)
        deallocate(vv_new2)
        deallocate(numerator)
        deallocate(tmp1)
        deallocate(other_index)
        deallocate(kappa)
        deallocate(kappaz)
     end if
  end do
  return
end subroutine expect_v_logistic



! subroutine to get terms in seff_beta
subroutine seff_beta_terms_logistic(beta,ynew,z,sbeta_terms)
  use commondata, only : n,maxm,lb,&
       n1,m1
  implicit none
  double precision,dimension(n,maxm),intent(in) :: ynew
  double precision,dimension(n,maxm,lb),intent(in) :: z
  double precision,dimension(lb),intent(in) :: beta
  double precision,dimension(lb,n),intent(out) :: sbeta_terms
  double precision,allocatable,dimension(:,:) :: Ainv
  double precision,allocatable,dimension(:) :: Vi,tmp1
  double precision,dimension(lb) :: tmp2
  double precision,dimension(n,maxm) :: ev
  integer :: i
  sbeta_terms=0.

  ! compute E(V_i|W_i,Z_i;beta)
  call expect_v_logistic(beta,ynew,z,ev)
  !print*,'ev=',ev

  ! form terms
  do i=1,n1
     if(m1(i).gt.1.) then
        ! allocate dimensions
        allocate(Ainv(m1(i),m1(i)))
        allocate(Vi(m1(i)))
        allocate(tmp1(m1(i)))

        ! get Ainv
        call get_Ainv(m1(i),Ainv)

        ! set Vi
        Vi = ynew(i,1:m1(i))
        Vi(1)=0.
        !print*,'Vi=',Vi
        !print*,'ev=',ev(i,1:m1(i))
        tmp1=matmul(Ainv,(Vi(1:m1(i))-ev(i,1:m1(i))))
        tmp2=matmul(transpose(z(i,1:m1(i),:)),tmp1)
        !print*,'tmp1=',tmp1
        !print*,'i=',i,'tmp2=',tmp2
        sbeta_terms(:,i)=tmp2

        !deallocate dimensions
        deallocate(Ainv)
        deallocate(Vi)
        deallocate(tmp1)
     end if
  end do
  return
end subroutine seff_beta_terms_logistic



! subroutine to estimate beta
subroutine estbeta(lb,betaest,yout,extflag)
  use commondata, only : tol,n,wv,nw,&
       ynew_use_logistic,z_use
  implicit none
  integer,intent(in) :: lb
  double precision,dimension(lb),intent(inout) :: betaest
  double precision,dimension(lb),intent(out) :: yout
  integer, intent(out) :: extflag
  double precision,dimension(lb,n) :: sbeta_terms
  double precision :: checka
  integer :: i,kk,eflaga

  if(maxval(abs(betaest)).ge.2.) then
     yout=999
     return
  end if

  !-----------------------------------------------------
  ! we use ynew_use as "y" values.
  ! when there is no censored data, ynew_use= original y.
  !--------------------------------------------------------

  call seff_beta_terms_logistic(betaest,ynew_use_logistic,&
       z_use,sbeta_terms)

  !print*,'beta=',betaest
  yout=sum(sbeta_terms,dim=2)

10 return
end subroutine estbeta



! subroutine to compute Ainv
subroutine get_Ainv(m,Ainv)
  implicit none
  integer, intent(in) :: m
  double precision,dimension(m,m),intent(out) :: Ainv
  call get_diag(dble(1),m,Ainv)
  Ainv(1,2:m) = -1.
  return
end subroutine get_Ainv

! subroutine to create a diagonal matrix
subroutine get_diag(d,n,diag_matrix)
  implicit none
  integer, intent(in) :: n
  double precision,intent(in) :: d
  double precision,dimension(n,n),intent(out) :: diag_matrix
  integer :: i
  diag_matrix = 0.
  do i=1,n
     diag_matrix(i,i)=1.
  end do
  diag_matrix=diag_matrix * d
  return
end subroutine get_diag

! subroutine to create a diagonal matrix
subroutine diag_elements(n,matrix,diag)
  implicit none
  integer, intent(in) :: n
  double precision,dimension(n,n),intent(in) :: matrix
  double precision,dimension(n),intent(out) :: diag
  integer :: i
  do i=1,n
     diag(i)=matrix(i,i)
  end do
  return
end subroutine diag_elements




! subroutine to do estimation (using imputation)
subroutine estimation_logistic(betat,ynew,z,eflag,betaest,var,var_beta)
  use commondata, only : tol,n,maxm,lb,nw,&
       ynew_use_logistic,z_use,n1
  implicit none
  double precision,dimension(lb),intent(in) :: betat
  double precision,dimension(n,maxm),intent(in) :: ynew
  double precision,dimension(n,maxm,lb),intent(in) :: z
  integer,intent(out) :: eflag
  double precision,dimension(lb),intent(out) :: betaest,var_beta
  double precision,dimension(lb,lb),intent(out) :: var
  double precision,dimension(lb) :: yout
  double precision,dimension(nw) :: wv
  integer :: tt,rr,ll,eflaga,lb1,jj,ii
  external estbeta

  !------------------
  ! terms to use
  !------------------
  ynew_use_logistic=ynew
  z_use=z

  !--------------------------
  ! do estimation procedure
  !--------------------------
  betaest=betat
  call hybrd10(estbeta,lb,betaest,yout,tol,eflag,wv,nw)

  if(eflag.ne.1) then
     goto 20
  end if

  !--------------------
  ! estimate variance
  !--------------------
  call getvar(lb,betaest,var,var_beta)

20 return
end subroutine estimation_logistic


subroutine getvar(lb,betaest,var,var_beta)
  use commondata, only : tol,n,n1,maxm,wv,nw,&
       ynew_use_logistic,z_use
  implicit none
  integer,intent(in) :: lb
  double precision,dimension(lb),intent(in) :: betaest
  double precision,dimension(lb),intent(out) :: var_beta
  double precision,dimension(lb,lb),intent(out) :: var
  double precision,dimension(lb,lb) :: tmp,var_tmp
  double precision,dimension(lb,lb) :: varmat_beta
  double precision,dimension(lb,n) :: sbeta_terms
  integer :: i,xx

  call seff_beta_terms_logistic(betaest,ynew_use_logistic,&
       z_use,sbeta_terms)

  !----------------------
  ! get variance matrix
  !---------------------
  var_tmp=0.
  do i=1,n1
     call muly(lb,sbeta_terms(:,i),tmp)
     var_tmp=var_tmp+tmp
  end do
  call inv(lb,var_tmp,var)
  !print*,'var_tmp=',var_tmp

  !-----------------------
  ! get variance for beta
  !-----------------------
  call diag_elements(lb,var,var_beta)
  !print*,'varmat_beta=',varmat_beta

  return
end subroutine getvar



! multiply vector with itself: y=tt'
subroutine muly(lb,t,y)
  integer,intent(in) :: lb
  double precision,dimension(lb),intent(in) :: t
  double precision,dimension(lb,lb),intent(out) :: y
  integer :: i,j
  y=0
  do i=1,lb
     do j=1,lb
        y(i,j)=t(i)*t(j)
     end do
  end do
  return
end subroutine muly

! invert a matrix
subroutine inv(lb,A,A1)
  integer,intent(in) :: lb
  double precision,dimension(lb,lb),intent(in) :: A
  double precision,dimension(lb,lb),intent(out) :: A1
  integer :: i
  integer,dimension(lb) :: ipvt
  double precision,dimension(lb,lb) :: B
  double precision :: cond
  double precision,dimension(lb) :: wv
  A1=0
  do i=1,lb
     A1(i,i)=1
  end do
  call decomp(lb,lb,A,cond,ipvt,wv)
  do i=1,lb
     call solvels(lb,lb,A,A1(:,i),ipvt)
  end do
  return
end subroutine inv


! multiply vectors: y=st'
subroutine mul(lb,s,t,y)
  integer,intent(in) :: lb
  double precision,dimension(lb),intent(in) :: s,t
  double precision,dimension(lb,lb),intent(out) :: y
  integer :: i,j
  do i=1,lb
     do j=1,lb
        y(i,j)=s(i)*t(j)
     end do
  end do
  return
end subroutine mul

! set all elements of matrix to a particular value
subroutine setB(p,q,a,B)
  implicit none
  integer, intent(in) :: p,q
  double precision,intent(in) :: a
  double precision,dimension(p,q),intent(out) :: B
  integer :: ii,jj
  do ii=1,p
     do jj=1,q
        B(ii,jj) = a
     end do
  end do
  return
end subroutine setB


