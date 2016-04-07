      real*8 function fdelta(del, h)
c  computes the spherical B-spline fucntion for a given distance del.
c  The normalization is 1 in this version.  The equations have been
c  simplified in this version to save computation time, check notes, 
c  J.G, 2002.
c  Input:	del ---- distance from a node
c		h   ---- spline radius (or average between nodes)
c
      implicit real*8 (a-h, o-z)
      deltil=del/h
      if(deltil.gt.2.d0) then
	  fdelta=0.d0
      else if(deltil.le.1.d0) then
	  fdelta=(4.d0+3.d0*deltil*deltil*(deltil-2.d0))
      else
	  fdelta=(2.d0-deltil)**3.d0
      endif

      fdelta=fdelta*0.25
      return
      end
