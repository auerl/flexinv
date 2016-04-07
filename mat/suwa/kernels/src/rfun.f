      real*4 function rfun(iopt,k,rr)
c this provides a real*4 function shell for calling
c the radfun subroutine.
c
c  input:
c    iopt = choice of radial function
c    k    = order of function desired
c    rr   = normalized radius
c
c  output:
c    rfun = radial function value
c
c  calls:
c    radfun
c
      real*4 rr
      real*8 r,func(15)
c
      r = rr * 6371.d0
c
      call radfun_split(iopt,func,r)
      rfun = func(k+1)
c
      return
      end
