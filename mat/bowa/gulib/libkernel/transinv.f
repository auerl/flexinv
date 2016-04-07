	subroutine transinv(theta, phi, transmat, theta1, phi1)
c  for a given theta and phi, apply the transformation matrix
c  to the cartesian components of the location.  Output theta1
c  and phi1 are rotated locations in the path coordinates.
c  J.G, 2002.
	implicit double precision (a-h,o-z)
	dimension transmat(*)
	x = sin(theta)
	y = x*sin(phi)
	x = x*cos(phi)
	z = cos(theta)
c
	x1 = transmat(1)*x + transmat(4)*y + transmat(7)*z
	y1 = transmat(2)*x + transmat(5)*y + transmat(8)*z
	z1 = transmat(3)*x + transmat(6)*y + transmat(9)*z
c
	phi1 = 0d0
	theta1 = acos(z1)
	if(dabs(theta1).gt.1.0e-6) then
		if(x1.ne.0d0.or.y1.ne.0d0) phi1 = atan2(y1,x1)
		if(phi.lt.0d0) phi = phi+6.28318530
	endif
	return
	end
