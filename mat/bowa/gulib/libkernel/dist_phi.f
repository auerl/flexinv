      real*8 function dist_phi(th, ph1, ph2)
c
c  find distance between two points assuming point 2 has theta=90 deg, 
c  return in radians
c
	implicit double precision (a-h,o-z)
	dist_phi=acos(sin(th)*cos(ph1-ph2))
	return
	end
