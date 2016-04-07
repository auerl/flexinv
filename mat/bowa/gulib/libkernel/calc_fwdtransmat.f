	subroutine calc_fwdtransmat(eqtheta, eqphi, sttheta, stphi, transmat)
c
c calculates the forward transformation matrix which rotates
c the given path to the equator.
c ax, ay, az, bx, by, bz, cx, cy, cz are the 9 elements in the matrix
c input  :   eqtheta, eqphi, sttheta, stphi
c output :   transmat(9)
c
	implicit double precision (a-h,o-z)
	dimension transmat(9)
c
	ax = sin(eqtheta)
	ay = ax*sin(eqphi)
	ax = ax*cos(eqphi)
	az = cos(eqtheta)
c	
	bx = sin(sttheta)
	by = bx*sin(stphi)
	bx = bx*cos(stphi)
	bz = cos(sttheta)
c
	call outprod(ax, ay, az, bx, by, bz, cx, cy, cz)
	call outprod(cx, cy, cz, ax, ay, az, bx, by, bz)
c	
	call normalize(ax, ay, az)
	call normalize(bx, by, bz)
	call normalize(cx, cy, cz)
	transmat(1) = ax
	transmat(2) = ay
	transmat(3) = az
	transmat(4) = bx
	transmat(5) = by
	transmat(6) = bz
	transmat(7) = cx
	transmat(8) = cy
	transmat(9) = cz
	return
	end
	
	subroutine outprod(a1, a2, a3, b1, b2, b3, c1, c2, c3)
	implicit double precision (a-h,o-z)
c
	c1 = a2*b3-a3*b2
	c2 = a3*b1-a1*b3
	c3 = a1*b2-a2*b1
	return
	end
		
	subroutine normalize(x, y, z)
	implicit double precision (a-h,o-z)
c	
	temp = x*x+y*y+z*z
	if(temp.gt.1.0e-20) then
		temp = dsqrt(temp)
		if(temp.eq.0.d0) stop 'error in calc_fwdtransmat.f!'
		temp1 = 1.0/temp
		x = x*temp1
		y = y*temp1
		z = z*temp1
	else
		x = 2.d0
		y = 0.d0
		z = 0.d0
	endif
	return
	end	
