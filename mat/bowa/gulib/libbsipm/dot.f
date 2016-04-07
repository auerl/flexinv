      subroutine dot(a,c,ncol1,ncol2,dotp)
      dimension a(1),c(1)
c
      dotp = 0.0
      do  i=ncol1,ncol2
c	write(*,*) i, a(i), c(i)
         dotp = dotp + a(i) * c(i)
      enddo
c
      return
      end
