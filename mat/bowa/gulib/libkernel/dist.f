      real*8 function dist(scol, slon, rcol, rlon)
c
c  find distance between source and receiver, return in radians
c
	implicit double precision (a-h,o-z)
	dist=acos(cos(scol)*cos(rcol)+sin(scol)*sin(rcol)*cos(rlon-slon))
	return
	end
