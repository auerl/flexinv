c
c   drsple returns the value of the function y(x) evaluated at point s
c   using the cubic spline coefficients computed by drspln and saved in
c   q. assuming x(1)<=x(2)... <=x(n). If s is outside the interval
c   (x(i1),x(i2)) rsple extrapolates
c   using the first or last interpolation polynomial.  the arrays must
c   be dimensioned at least - x(i2), y(i2), and q(i2)(3).
c
	double precision function drsple(i1, i2, x, y, q, s)
	implicit real*8 (a-h,o-z)
	real*8 x(i2), y(i2), q(3, i2), s, h
	integer i1, i2, i
	
	i = 2
	ii = i2-1
c make sure within bounds
	i = max(i, i1)
	i = min(i, ii)
	if(s.lt.x(i1)) i=i1
	if(s.gt.x(i2)) i=ii
	if(s.ge.x(i1).and.s.le.x(i2)) then
10		if(s.lt.x(i+1)) goto 20
			i = i+1
			goto 10
20		continue
30		if(s.ge.x(i)) goto 40
			i = i-1
			goto 30
40		continue
	endif
	h = s-x(i)
	drsple = y(i)+h*(q(1,i)+h*(q(2,i)+h*q(3,i)))
	return
	end
	        
