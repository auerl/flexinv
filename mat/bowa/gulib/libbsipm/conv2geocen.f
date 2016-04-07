	subroutine conv2geocen(lat, lon, theta, phi)
c
c convert a geographic coordinate to geocentric coordinate
c
	implicit double precision (a-h,o-z)
	real*8 lat, lon, theta, phi
	data fac /0.993305621334896/
	data rad /0.017453292519943/
	write(*,*) 'lat=', lat,' lon=', lon
	theta = (90.0-lat)*rad
	write(*,*) 'theta = ', theta
	if(theta.gt.1.0e-20) then 
		theta = 1.5707963-atan2(fac*cos(theta),sin(theta))
	else
		theta = 1.5707963-atan2(fac*cos(theta),1.0e-20)
	endif
	phi = lon
	if(phi.lt.0d0) phi = phi+360.0
	phi = phi*rad
	return
	end
 
