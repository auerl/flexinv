	subroutine dbspfun(ib, x, xknt, nknt, b)
c  calculates radial B-spline function using double precision
c  Input:	ib = index of B-splines which x falls
c	         x = radius of the value to be calculates
c	      xknt = array of radial knots
c	      nknt = number of radial knots
c  Output:
c	       b[] = B-spline values, only 4 non-zero ones
c
	implicit double precision (a-h,o-z)
	dimension b(1),xknt(-1:nknt+2)
	data ib0/0.d0/
	save ib0,a1,a2,a3,a4
	
	if(ib.ne.ib0)then
	   temp=1.d0/(xknt(ib+1)-xknt(ib))
	   a1=temp/(xknt(ib+1)-xknt(ib-1))
	   a3=temp/(xknt(ib+2)-xknt(ib))
           a4=a3/(xknt(ib+3)-xknt(ib))
	   temp=xknt(ib+2)-xknt(ib-1)
	   a2=a1/temp
	   a3=a3/temp
	   a1=a1/(xknt(ib+1)-xknt(ib-2))
	   ib0=ib
	endif
	  
	x1=x-xknt(ib-1)
	x2=x-xknt(ib)
	x3=x-xknt(ib+1)
	x4=x-xknt(ib+2)
  
	temp=x3*x3*a1
	temp1=x3*x1*a2+x2*x4*a3
	b(1)=-temp*x3
	b(2)=(x-xknt(ib-2))*temp+x4*temp1
	temp=x2*x2*a4
	b(4)=x2*temp
	b(3)=-x1*temp1-(x-xknt(ib+3))*temp 
c
c ... for 1st derivative = 0 case
c
c	if(ib.eq.1)then
c	  b(2)=b(2)+b(1)
c	endif
c	if(ib.eq.(nknt-1))then
c	   b(3)=b(3)+b(4)
c	endif
c
c...for 2nd derivative = 0
c
	if(ib.eq.1)then
	  b(3)=b(3)+b(2)/3.d0
	  b(2)=b(2)/1.5d0+b(1)
	endif
	if(ib.eq.2)then
	  b(2)=b(2)+b(1)/3.d0
	  b(1)=b(1)/1.5d0
	endif
	if(ib.eq.(nknt-1))then
	   b(2)=b(2)+b(3)/3.d0
	   b(3)=b(3)/1.5d0+b(4)
	endif
	
	if(ib.eq.(nknt-2))then
	   b(3)=b(3)+b(4)/3.d0
	   b(4)=b(4)/1.5d0
	endif

	return
	end

