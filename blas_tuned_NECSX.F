#undef __UNROLL
c#define __SAFE
#undef __SAFE
c----------------------------------------------------------------------
      subroutine  dscal(n,da,dx,incx)
c----------------------------------------------------------------------
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
#ifdef __SAFE
      if( n.gt.0 .and. incx.gt.0 ) then 
#endif
         if(incx.eq.1) then
c
c           code for increment equal to 1
c
c
c           clean-up loop
c
#ifdef __UNROLL
            m = mod(n,5)
            if( m .ne. 0 ) then
               do i = 1,m
                 dx(i) = da*dx(i)
               end do
	    end if
            if( n .ge. 5 ) then
               mp1 = m + 1
               do i = mp1,n,5
                 dx(i) = da*dx(i)
                 dx(i + 1) = da*dx(i + 1)
                 dx(i + 2) = da*dx(i + 2)
                 dx(i + 3) = da*dx(i + 3)
                 dx(i + 4) = da*dx(i + 4)
               end do
	    endif
#else
            do i = 1,n
              dx(i) = da*dx(i)
            end do
#endif
         else
c
c           code for increment not equal to 1
c
            nincx = n*incx
            do i = 1,nincx,incx
              dx(i) = da*dx(i)
            end do
         endif
#ifdef __SAFE
      endif
#endif
      return
      end
      
c----------------------------------------------------------------------
      subroutine  zdscal(n,da,zx,incx)
c----------------------------------------------------------------------
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double complex zx(*)
      double precision da
      integer i,incx,ix,n
c
#ifdef __SAFE
      if( n.le.0 .or. incx.le.0 )return
#endif
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      do 10 i = 1,n
        zx(ix) = dcmplx(da,0.0d0)*zx(ix)
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        zx(i) = dcmplx(da,0.0d0)*zx(i)
   30 continue
      return
      end
      
c----------------------------------------------------------------------
      subroutine  dcopy(n,dx,incx,dy,incy)
c----------------------------------------------------------------------
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
#ifdef __SAFE
      if(n.gt.0)then
#endif
        if(incx.ne.1.or.incy.ne.1) then
c
c        code for unequal increments or equal increments
c          not equal to 1
c
	   ix = 1
	   iy = 1
	   if(incx.lt.0)ix = (-n+1)*incx + 1
	   if(incy.lt.0)iy = (-n+1)*incy + 1
	   do i = 1,n
             dy(iy) = dx(ix)
             ix = ix + incx
             iy = iy + incy
	   end do
         else
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
#ifdef __UNROLL
            m = mod(n,7)
	    if( m .ne. 0 ) then
	      do i = 1,m
        	dy(i) = dx(i)
	      end do
	    end if 
	    if( n .ge. 7 ) then
              mp1 = m + 1
	      do i = mp1,n,7
        	dy(i) = dx(i)
        	dy(i + 1) = dx(i + 1)
        	dy(i + 2) = dx(i + 2)
        	dy(i + 3) = dx(i + 3)
        	dy(i + 4) = dx(i + 4)
        	dy(i + 5) = dx(i + 5)
        	dy(i + 6) = dx(i + 6)
	      end do
            end if
#else
	    do i = 1,n
              dy(i) = dx(i)
	    end do
#endif
        end if
#ifdef __SAFE
      end if
#endif
      return
      end

c----------------------------------------------------------------------
      subroutine  zcopy(n,dx,incx,dy,incy)
c----------------------------------------------------------------------
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      complex*16 dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
#ifdef __SAFE
      if(n.gt.0)then
#endif
        if(incx.ne.1.or.incy.ne.1) then
c
c        code for unequal increments or equal increments
c          not equal to 1
c
	   ix = 1
	   iy = 1
	   if(incx.lt.0)ix = (-n+1)*incx + 1
	   if(incy.lt.0)iy = (-n+1)*incy + 1
	   do i = 1,n
             dy(iy) = dx(ix)
             ix = ix + incx
             iy = iy + incy
	   end do
         else
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
#ifdef __UNROLL
            m = mod(n,7)
	    if( m .ne. 0 ) then
	      do i = 1,m
        	dy(i) = dx(i)
	      end do
	    end if 
	    if( n .ge. 7 ) then
              mp1 = m + 1
	      do i = mp1,n,7
        	dy(i) = dx(i)
        	dy(i + 1) = dx(i + 1)
        	dy(i + 2) = dx(i + 2)
        	dy(i + 3) = dx(i + 3)
        	dy(i + 4) = dx(i + 4)
        	dy(i + 5) = dx(i + 5)
        	dy(i + 6) = dx(i + 6)
	      end do
            end if
#else
	    do i = 1,n
              dy(i) = dx(i)
	    end do
#endif
        end if
#ifdef __SAFE
      end if
#endif
      return
      end
      
c----------------------------------------------------------------------
      double precision function ddot(n,dx,incx,dy,incy)
c----------------------------------------------------------------------
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
#ifdef __SAFE
      if(n.le.0)return
#endif
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      
c----------------------------------------------------------------------
      subroutine  drot (n,dx,incx,dy,incy,c,s)
c----------------------------------------------------------------------
c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp,c,s
      integer i,incx,incy,ix,iy,n
c
#ifdef __SAFE
      if(n.le.0)return
#endif
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
   20 do 30 i = 1,n
        dtemp = c*dx(i) + s*dy(i)
        dy(i) = c*dy(i) - s*dx(i)
        dx(i) = dtemp
   30 continue
      return
      end
      
                  
