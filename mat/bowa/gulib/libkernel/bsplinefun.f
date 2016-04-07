	subroutine bsplinefun(ib,x,xknt,ibot,itop,b)
c calculates radial B-spline function.
c   input:   ib = the index of B-splines in which x falls, begins with 1
c            x  = the radius of the value to be calculated
c	     ibot = index of the lower knots (only splines between
c		    ibot and itop are considered, designed to
c		    deal with either continuous or discontinuous models
c		    Y.G. 2002.
c	     itop = index of top knots 
c   output: 
c	     b = nonzero splines saved in b(1) to b(4)
c   	    
c
        implicit double precision (a-h,o-z)

c *** basic parameters, naming is obsolete
      parameter(maxfreq=8)
      parameter(maxnode=maxfreq*maxfreq*10+2)
      parameter(maxrknot=50)
      parameter(maxparm=maxnode*maxrknot)
      data pi /3.141592653579/
      data rad/57.2957795130824/
      data reprad /57.295779579/    



	dimension xknt(maxrknot), tempxknt(0:maxrknot)
	dimension b(4)
	data isetknt/1/
	data ib0/0/,iumlm0/0/
	save ib0,a1,a2,a3,a4

c
c  copy input array into an temp array that start with the 0th element
c  and work with that for the rest of the indexing.
c
	b(1)=0.d0
	b(2)=0.d0
	b(3)=0.d0
	b(4)=0.d0
	
	if(ibot.gt.itop) stop 'wrong index in bsplinefun.f'
        if(x.lt.xknt(ibot).or.x.gt.xknt(itop)) goto 10
	knt=itop-ibot+1	! number of knots
	do i=0, knt-1
		tempxknt(i)=xknt(i+ibot)
	enddo
	kntm1=knt-1
	if(ib.eq.1.or.ib.eq.2) then
		x0=tempxknt(0)
		x1=tempxknt(0)
		x4=tempxknt(ib+1)
		x5=tempxknt(ib+2)
	else if(ib.eq.(kntm1-1).or.ib.eq.kntm1) then
    		x0=tempxknt(ib-3)
    		x1=tempxknt(ib-2)
    		x4=tempxknt(kntm1)
    		x5=tempxknt(kntm1)
	else
    		x0=tempxknt(ib-3)
    		x1=tempxknt(ib-2)
    		x4=tempxknt(ib+1)
    		x5=tempxknt(ib+2)
	endif
	x2=tempxknt(ib-1)
	x3=tempxknt(ib)
c
c Always computes the following	
	temp=1.d0/(x3-x2)
	a1=temp/(x3-x1)
	a3=temp/(x4-x2)
	a4=a3/(x5-x2)
	temp=x4-x1
	a2=a1/temp
	a3=a3/temp
	a1=a1/(x3-x0)
	
	xx1=x-x1
	xx2=x-x2
	xx3=x-x3
	xx4=x-x4
	
	temp=xx3*xx3*a1
	temp1=xx3*xx1*a2+xx2*xx4*a3
	b(1)=-temp*xx3
	b(2)=(x-x0)*temp+xx4*temp1
	temp=xx2*xx2*a4
	b(4)=xx2*temp
	b(3)=-xx1*temp1-(x-x5)*temp
c
c ...for 2nd derivative = 0
c
	if(ib.eq.1) then
		b(3)=b(3)+b(2)/3.d0
		b(2)=b(2)/1.5d0+b(1)
		b(1)=0.d0
	endif
	if(ib.eq.2) then
		b(2)=b(2)+b(1)/3.d0
		b(1)=b(1)/1.5d0
	endif
	if(ib.eq.kntm1) then
		b(2)=b(2)+b(3)/3.d0
		b(3)=b(3)/1.5d0+b(4)
		b(4)=0.d0
	endif
	if(ib.eq.(kntm1-1)) then
		b(3)=b(3)+b(4)/3.d0
		b(4)=b(4)/1.5d0
	endif

10	return
	end
