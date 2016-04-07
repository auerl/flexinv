      subroutine dot_splines(a,inda,n,coef,dotp)
      dimension a(1),coef(1)
      integer inda(1)
c
      dotp=0.0
      do i=1, n
	 ind=inda(i)
         dotp = dotp + a(i) * coef(ind)
      enddo
      return
      end
