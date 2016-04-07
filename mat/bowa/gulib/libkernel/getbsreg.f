      subroutine getbsreg(r,rk,i0,i1,ii)
c
c  Get the spline which contains the given radius
c  Assuming the splines are ordered from bottom to top. 
c  Input:	r --- radius
c		rk --- array of radial spline nodes
c		i0 --- beginning index
c		i1 --- ending index
c  Output:
c		ii --- index for the given radius
c
      implicit double precision (a-h,o-z)
      dimension rk(1)
      
      ii=0

      if(i0.ge.i1) stop 'wrong spline knot index!!'
      if(r.le.rk(i0).or.r.gt.rk(i1)) then
	 return
      endif	

      j=1
      do i=i0, i1-1
	 if(r.gt.rk(i).and.r.le.rk(i+1)) then
		ii=j
		goto 10
	 endif
	 j=j+1
      enddo
10    return
      end
