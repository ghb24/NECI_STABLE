   program grid1d
    implicit none
    integer, parameter :: ngrid1 = 40,maxq = 1000000
    double precision :: xygrid(ngrid1),matrix1e(ngrid1,ngrid1),overlap(ngrid1,ngrid1)
    double precision :: orbital(ngrid1,ngrid1),ywork(maxq),energies1e(ngrid1)
    double precision, dimension(ngrid1,ngrid1,ngrid1) :: wavefunction1,wavefunction2,wavefunction3
    double precision:: tempint(ngrid1,ngrid1)
    double precision, allocatable :: integral3e(:,:,:,:,:,:)
    double precision, allocatable :: int33(:,:,:,:,:,:)
    double precision, parameter :: pi = 3.14159265358979d0
    double precision gridmin,gridmax,midpoint,step,x,x1,x2,x3,f12(2),f13(2),test,test2,test3,integral_value
    double precision try,a,b,c,d,diff,try1,try2
   
    integer :: i,j,k,l,m,n,info,ii,jj,kk,ll,ip,iq,nmax,ic,ic1,ir,is,it,iu,m1,m2,m3,p,q,r,s,t,u,i3(3),j3(3)
    character*10 cdist,rfile


   gridmin =  0.d0
   gridmax =  10.d0
   midpoint = (gridmax+gridmin)/2.d0
   step = (gridmax-gridmin)/(ngrid1)
   do i = 1,ngrid1
    xygrid(i) = gridmin+(i-1)*step+step/2.d0
   enddo
   matrix1e(:,:) = 0.d0
   do i = 1,ngrid1
    x = xygrid(i)
    matrix1e(i,i) = 1.d0/step**2
    if(i.ne.1) matrix1e(i-1,i) = -1.d0/(2*step**2)
    if(i.ne.1) matrix1e(i,i-1) = -1.d0/(2*step**2)
    if(i.ne.ngrid1) matrix1e(i+1,i) = -1.d0/(2*step**2)
    if(i.ne.ngrid1) matrix1e(i,i+1) = -1.d0/(2*step**2)
! here is the 1e potential, change here to generate different orbitals
! He atom at midpoint, this generates symmetric orbitals
    matrix1e(i,i) = matrix1e(i,i) - 2.0/dsqrt((x-midpoint)**2+1)
! He atom moved from  midpoint, this generates orbitals which don't have symmetry
!    matrix1e(i,i) = matrix1e(i,i) - 2.0/dsqrt((x-midpoint+1.0)**2+1)
   enddo
   call dsyev('V','U',ngrid1,matrix1e,ngrid1,energies1e,ywork,maxq,info)
   do i = 1,ngrid1
    do j = 1,ngrid1
     orbital(i,j) = matrix1e(i,j)
    enddo
   enddo


! prints out the 7 lowest orbitals
   open(7,file='orbitals.dat')
   do i = 1,ngrid1
     write(7,*) xygrid(i),orbital(i,1:7)
   enddo
   close(7)
   write(6,*) "orbitals printed to file 'orbitals.dat', to visualise orbitals, e.g. plot 'orbitals.dat' u 1:2 w l ,'' u 1:3 w l"

   nmax = 7
   allocate(integral3e(nmax,nmax,nmax,nmax,nmax,nmax))
   allocate(int33(nmax,nmax,nmax,nmax,nmax,nmax))
   integral3e = 0.d0
   do i = 1,ngrid1
    x1 = xygrid(i)
    tempint(:,:) = 0.d0
    do j = 1,ngrid1
     x2 = xygrid(j)
     f12 = tau_r12(x1,x2)
     do ip = 1,nmax
      do iq = 1,nmax
       tempint(ip,iq) = tempint(ip,iq) + f12(2)*orbital(j,ip)*orbital(j,iq)
      enddo
     enddo
    enddo
    do ip = 1,nmax
     do iq = 1,nmax
      do ir = 1,nmax
       do is = 1,nmax
        do it = 1,nmax
         do iu = 1,nmax
           integral3e(ip,ir,it,iq,is,iu) = integral3e(ip,ir,it,iq,is,iu) +   &
          tempint(ir,is)*tempint(it,iu)*orbital(i,ip)*orbital(i,iq)
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
    do ip = 1,nmax
     do iq = 1,nmax
      do ir = 1,nmax
       do is = 1,nmax
        do it = 1,nmax
         do iu = 1,nmax
           int33(ip,ir,it,iq,is,iu) =  1.d0/3.d0*integral3e(ip,ir,it,iq,is,iu) & 
                                     + 1.d0/3.d0*integral3e(it,ip,ir,iu,iq,is) &
                                     + 1.d0/3.d0*integral3e(ir,it,ip,is,iu,iq)
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo

! integral3e treats electron 1 as special
! int33 is the symmetrized integral 
! order of indices is the same for both and corresponds to physicist ordering (chemists notation doesnt exist yet!)
! e.g. 1L = orbital_1 left , 2R = orbital_2 right
! int33(1L,2L,3L,1R,2R,3R)



!   open(16,file='ints3e.dat',form='unformatted')
!   write(6,*) 'all 3e integrals written unformatted to ints3e.dat'
!   write(16) nmax
   open(16,file='ints3e.dat',form='formatted')
   write(6,*) 'all 3e integrals written formatted to ints3e.dat'
   do ip = 1,nmax
    do iq = 1,nmax
     do ir = 1,nmax
      do is = 1,nmax
       do it = 1,nmax
        do iu = 1,nmax
!          write(16) ip,iq,ir,is,it,iu,int33(ip,iq,ir,is,it,iu)
          write(16,*) int33(ip,iq,ir,is,it,iu),ip,iq,ir,is,it,iu
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo


!  Change specification of indices of the determinants here
   i3 = (/2,1,3/)
   j3 = (/4,5,7/)
!
   do i = 1,3
    if(i3(i).gt.nmax) then
      write(6,*) 'index of left determinant chosen greater than nmax,',i3(i),nmax
      write(6,*) 'either increase nmax or change choice of determinant!'
      stop
    endif
    if(j3(i).gt.nmax) then
      write(6,*) 'index of right determinant chosen greater than nmax,',j3(i),nmax
      write(6,*) 'either increase nmax or change choice of determinant!'
      stop
    endif
   enddo

   do m1 = 1,ngrid1
    do m2 = 1,ngrid1
     do m3 = 1,ngrid1
      wavefunction1(m1,m2,m3) =  (orbital(m1,i3(1))*(orbital(m2,i3(2))*orbital(m3,i3(3))-orbital(m2,i3(3))*orbital(m3,i3(2))) &
                               -  orbital(m2,i3(1))*(orbital(m1,i3(2))*orbital(m3,i3(3))-orbital(m1,i3(3))*orbital(m3,i3(2))) &
                    +  orbital(m3,i3(1))*(orbital(m1,i3(2))*orbital(m2,i3(3))-orbital(m1,i3(3))*orbital(m2,i3(2))) )/sqrt(6.d0)
      wavefunction2(m1,m2,m3) =  (orbital(m1,j3(1))*(orbital(m2,j3(2))*orbital(m3,j3(3))-orbital(m2,j3(3))*orbital(m3,j3(2))) &
                               -  orbital(m2,j3(1))*(orbital(m1,j3(2))*orbital(m3,j3(3))-orbital(m1,j3(3))*orbital(m3,j3(2))) &
                    +  orbital(m3,j3(1))*(orbital(m1,j3(2))*orbital(m2,j3(3))-orbital(m1,j3(3))*orbital(m2,j3(2))) )/sqrt(6.d0)
     enddo
    enddo 
   enddo

   integral_value  = 0.d0
   do i = 1,ngrid1
    x1 = xygrid(i)
    do j = 1,ngrid1
     x2 = xygrid(j)
     do k = 1,ngrid1
      x3 = xygrid(k)
      f12 = tau_r12(x1,x2)
      f13 = tau_r12(x1,x3)
      integral_value = integral_value + f12(2)*f13(2)*wavefunction1(i,j,k)*wavefunction2(i,j,k)
     enddo
    enddo
   enddo
   
! value as I would calculate slater-condon rule in 3e case from <Phi_Left|L |sum_u sigma_u P_u Phi_Right>
! where Phi_Left is the string of orbitals on the left 
! and u is all possible permutations and sigma_u is the sign of permutation eg +E, -P_12, + P_123  (E is identity)
!   try  = 0
!   p = i3(1)
!   q = i3(2)
!   r = i3(3)
!   s = j3(1)
!   t = j3(2)
!   u = j3(3)
!
!   try = try +  int33(p,q,r,s,t,u)
!   try = try -  int33(p,q,r,t,s,u)
!   try = try -  int33(p,q,r,u,t,s)
!   try = try -  int33(p,q,r,s,u,t)
!   try = try +  int33(p,q,r,u,s,t)
!   try = try +  int33(p,q,r,t,u,s)
!   write(6,*) 'Check of integral, difference =', integral_value -try

   write(6,*) '3 indices of Psi_L',i3
   write(6,*) '3 indices of Psi_R',j3
   write(6,*) 'Value of <Psi_L|L_123|Psi_R>',integral_value





contains

  pure function tau_r12(x1,x2)
    double precision, intent(in) :: x1,x2
    double precision tau_r12(2)
    double precision alpha

    alpha = 1.d0
    tau_r12(1) = exp(-alpha*(x1-x2)**2)
    tau_r12(2) = -2*alpha*(x1-x2)*exp(-alpha*(x1-x2)**2)

  end function tau_r12

end program
