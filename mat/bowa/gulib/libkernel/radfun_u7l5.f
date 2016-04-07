      subroutine radfun_u7l5(radius, chebyshev)
      implicit double precision (a-h, o-z)
      dimension cheby(14)
      dimension chebyshev(1)
      parameter (r0=6371.d0)
      parameter (rmoho=6371.d0-24.4d0)
      parameter (r670=6371.d0-670.d0)
      parameter (rcmb=3480.d0)
      logical upper, lower
c
c---- evaluate model at this depth
c
      do i=1, 14
	 chebyshev(i)=0.d0
	 cheby(i)=0.d0
      enddo 
      upper=.false.
      lower=.false.
      if(radius.gt.rcmb.and.radius.lt.r670) then
        lower=.true.
      else if(radius.ge.r670.and.radius.lt.rmoho) then
        upper=.true.
      endif
      if(upper) then
              u=(radius+radius-rmoho-r670)/(rmoho-r670)
              call chebyfun(u,13,cheby)
      else if(lower) then
              u=(radius+radius-r670-rcmb)/(r670-rcmb)
              call chebyfun(u,13,cheby)
      else
	      return
      endif
      if(upper) then
	      do i=1, 8
		 chebyshev(i) = cheby(i)
	      enddo
      else if(lower) then
	      do i=1, 6
		 chebyshev(i+8) = cheby(i)
	      enddo
      endif
      return
      end

      subroutine chebyfun(u,kmax,f)
      implicit double precision (a-h, o-z)
      dimension chebycoeff(0:13),f(0:kmax)
      data  chebycoeff /
     . 0.70710678118655,1.2247448713916,1.0350983390135,1.0145993123918,
     . 1.00803225754840,1.0050890913907,1.0035149493262,1.0025740068320,
     . 1.00196657023780,1.0015515913133,1.0012554932754,1.0010368069141,
     . 1.00087070107920,1.0007415648034 /
       
      if(kmax.gt.13)then
         write(*,"(' kmax exceeds the limit in chebyfun')")
         stop
      endif
       
      f(0)=1.d0
      f(1)=u
      twou=2.d0*u
       
      do k=2,kmax
         f(k) = twou*f(k-1)-f(k-2)
      enddo
       
      do k=0,kmax
         f(k)=f(k)*chebycoeff(k)
      enddo
       
      return
      end
