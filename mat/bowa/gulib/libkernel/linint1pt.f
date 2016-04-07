	subroutine linint1pt(x1,x2,y1,y2,x,yout)
	implicit real*8 (a-h,o-z)
	a=(y2-y1)/(x2-x1)
	b=y1-a*x1
	yout=a*x+b
	return
	end
