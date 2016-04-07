	subroutine conv2geocen(xlat, xlon, theta, phi)
c
c convert a geographic coordinate to geocentric coordinate
c
	implicit double precision (a-h,o-z)
	data fac /0.993305621334896/
	data rad /0.017453292519943/
	theta = (90.0-xlat)*rad
	if(theta.gt.1.0e-20) then 
		theta = 1.5707963-atan2(fac*cos(theta),sin(theta))
	else
		theta = 1.5707963-atan2(fac*cos(theta),1.0e-20)
	endif
	phi = xlon
	if(phi.lt.0.d0) phi = phi+360.0
	phi = phi*rad
	return
	end
 
