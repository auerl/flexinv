      subroutine xlinearint(x,x0,x1,y0,y1,xx)
c
c Linear interpolation, good enough for the purpose here, though cubic
c interpolation may be more accurate.  But this is time saving.
c
      implicit double precision (a-h, o-z)

      if(x0.eq.x1) then
	 xx=0.5d0*(y0+y1)
      else
	 xx=y0+(y1-y0)*(x-x0)/(x1-x0)
      endif
      return
      end
