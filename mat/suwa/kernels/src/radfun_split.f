
      subroutine radfun_split(iopt,func,r)

C radial function for a split B-Spline model.
c Returns values of radial functions at specified radius. 
c For some iopt's, it computes all radial orders (e.g. 4 for M84A) 
c even though the calling code may only use a single radial order 
c at a time.
c
c input:
c  iopt = choice of radial functions
c         1 = Legendre polynomials in upper mantle, as in M84A
c         2 =    "         "       "  lower   "   , as in LO2.56
c         3 = constants, 3 layer upper mantle
c         4 = constants, crustal portion of M84C
c         5 = constant, single layer upper mantle
c         6 = single layer crust, from Smith and Masters
c         7 = Tanimoto model, 10 depth knots, lin. interp. between knots
c         8 = 4 knots, linear interpolation between knots
c         9 = Tanimoto whole mantle model, 11 layers
c         10= Inoue et al. whole mantle model
c         11= 4 layer (5 knot) upper mantle (lin. interp. between knots)
c         12= Legendre polynomial l. m. with D" layer
c         13= 4 knots, linear interpolation between knots (PREM crust)
c         14= Chebyshev upper mantle k=0 to 4
c         15= Chebyshev lower mantle k=0 to 8
c         16= Chebyshev whole mantle k=0 to 8
c         17= 14 layer whole mantle
c         18= 14 layer upper mantle (Zhang's model - reduced from 26 layers)
c         19= 9 knots, B-spline  upper mantle
c         20= 12 knots, B-sline  lower mantle
c         21= 14 knots, B-sline  whole mantle
c   r    = radius (km)
c
c output:
c   func() = radial function values
c
      implicit real*8(a-h,o-z)
      dimension func(15),con(8),rtan(10),rtanw(12),rinoue(16),conc(14),
     &          rlmod(15),rz(15)
c
      data rcmb,r670,rmoho/3480.d0,5701.d0,6346.6d0/
      data rsol,rmidc,r220,r400 /6368.d0, 6356.d0, 6151.d0, 5971.d0/
c
c Legendre polynomial normalizing coefficients
c
      data con/0.7071067811865476d0,1.224744871391589d0,
     +  0.7905694150420949d0,0.9354143466934853d0,0.2651650429449553d0,
     +  0.29315098498896d0,0.15934435979977d0,0.17116329922036d0/
c
c Chebyshev polynomial normalizing coefficients
c
      data conc/0.70710678118655d0,1.2247448713916d0,1.0350983390135d0,
     &         1.0145993123918d0,1.0080322575484d0,1.0050890913907d0,
     &         1.0035149493262d0,1.0025740068320d0,1.0019665702378d0,
     &         1.0015515913133d0,1.0012554932754d0,1.0010368069141d0,
     &         1.0008707010792d0,1.0007415648034d0/
c
c knots for option 7
      data rtan/6338.d0,6258.d0,6155.d0,6053.d0,5950.d0,5825.d0,5700.d0,
     &          5613.d0,5483.d0,5396.d0/
c
c layer boundaries for option 9
      data rtanw/6371.d0,6151.d0,5971.d0,5701.d0,5349.d0,5087.d0,
     &           4816.d0,4555.d0,4283.d0,4012.d0,3741.d0,3480.d0/
c
c layer boundaries for option 10
      data rinoue/6342.d0,6293.d0,6223.d0,6133.d0,6023.d0,
     &            5893.d0,5742.d0,5571.d0,5380.d0,5168.d0,4936.d0,
     &            4683.d0,4411.d0,4118.d0,3805.d0,3471.d0/
c
c layer boundaries for option 17
      data rlmod/6346.6d0,6261.d0,6151.d0,6021.d0,5871.d0,
     &           5701.d0, 5471.d0,5221.d0,4971.d0,4721.d0,
     &           4471.d0, 4221.d0,3971.d0,3721.d0,3479.999d0/
c
c layer boundaries for option 18
      data rz/6371.0d0,6356.0d0,6346.6d0,6318.8d0,6291.0d0,
     &        6251.0d0,6211.0d0,6171.0d0,6131.0d0,6091.0d0,
     &        6051.0d0,6011.0d0,5971.0d0,5931.0d0,5870.999d0/      
c
c boundary of D'' layer of option 12
      data rlay/3700.d0/
c
c
      do 20 i=1,15
   20 func(i)=0.d0
c
c Legendre polynomial upper mantle, as in M84C
c
      if (iopt.eq.1) then
        if (r.gt.rmoho.or.r.lt.r670) return
        u = (r+r-rmoho-r670)/(rmoho-r670)
c
c Legendre polynomial lower mantle, as in L02.56
c
      elseif (iopt.eq.2) then
        if (r.gt.r670.or.r.lt.rcmb) return
        u = (r+r-r670-rcmb)/(r670-rcmb)
c
c 3 layer upper mantle, using 20 km Moho
c
      elseif (iopt.eq.3) then
        if (r.le.r670) return
        if (r.le.6351.d0.and.r.gt.r220) func(1) = 1.d0
        if (r.le.r220.and.r.gt.r400) func(2) = 1.d0
        if (r.le.r400.and.r.gt.r670) func(3) = 1.d0
        return
c
c M84C crustal function for shear velocity
c
      elseif (iopt.eq.4) then
        if (r.le.rmoho.or.r.gt.rsol) return
        if (r.gt.rmidc.and.r.le.rsol)  func(1) = -0.005254
        if (r.gt.rmoho.and.r.le.rmidc) func(1) = 0.001751
        return
c
c single layer upper mantle, PREM crustal thickness
c
      elseif (iopt.eq.5) then
        if (r.le.r670) return
        if (r.le.rmoho.and.r.gt.r670) func(1) = 1.d0
        return
c
c  single layer crust, moho at 20 km
c
      elseif (iopt.eq.6) then
        if (r.le.6351.d0) return
        if (r.le.6368.400000000001d0.and.r.gt.6351.d0) func(1) = 1.d0
        return
c
c 10 knot upper mantle, as in Tanimoto's DSXRG?
c
      elseif (iopt.eq.7) then
        if (r.gt.6338.d0.or.r.lt.5396.d0) return
        do 30 i=1,9
        j = i
  30    if (r.ge.rtan(i+1)) goto 40
        return
  40    u = (r - rtan(j+1)) / (rtan(j) - rtan(j+1))
        func(j)   = u
        func(j+1) = 1.d0 - u
        return
c
c  4 knot upper mantle, using Moho at 20 km 
c
      elseif (iopt.eq.8) then
        if (r.le.r670) return
        if (r.le.6351.d0.and.r.gt.r220) then
          u = (r - r220) / (6351.d0 - r220)
          func(1) = u
          func(2) = 1.d0 - u
        endif
        if (r.le.r220.and.r.gt.r400) then
          u = (r - r400) / (r220 - r400)
          func(2) = u
          func(3) = 1.d0 - u
        endif
        if (r.le.r400.and.r.gt.r670) then
          u = (r - r670) / (r400 - r670)
          func(3) = u
          func(4) = 1.d0 - u
        endif
        return
c
c 11 layer mantle, as in MDLSH
c
      elseif (iopt.eq.9) then
        if (r.gt.6371.d0.or.r.lt.3480.d0) return
        do 50 i=1,11
        j = i
  50    if (r.gt.rtanw(i+1)) goto 60
        return
  60    func(j) = 1.d0
        return
c
c 16 layer mantle, as in Inoue et al.
c
      elseif (iopt.eq.10) then
        if (r.gt.6342.d0.or.r.lt.3480.d0) return
        do 70 i=1,16
        j = i
  70    if (r.gt.rinoue(i+1)) goto 80
        return
  80    func(j) = 1.d0
        return
c
c 5 knot upper mantle, using PREM Moho depth
c
      elseif (iopt.eq.11) then
        if (r.le.r670) return
        if (r.le.rmoho.and.r.gt.6291.d0) then
          u = (r - 6291.d0) / (rmoho - 6291.d0)
          func(1) = u
          func(2) = 1.d0 - u
        elseif (r.le.6291.d0.and.r.gt.r220) then
          u = (r - r220) / (6291.d0 - r220)
          func(2) = u
          func(3) = 1.d0 - u
        elseif (r.le.r220.and.r.gt.r400) then
          u = (r - r400) / (r220 - r400)
          func(3) = u
          func(4) = 1.d0 - u
        elseif (r.le.r400.and.r.gt.r670) then
          u = (r - r670) / (r400 - r670)
          func(4) = u
          func(5) = 1.d0 - u
        endif
        return
c
c Legendre polynomial lower mantle, with a constant D'' layer
c
      elseif (iopt.eq.12) then
        if (r.gt.r670.or.r.lt.rcmb) return
        if (r.gt.rlay) then
          u = (r+r-r670-rcmb)/(r670-rcmb)
        else
          func(6) = 1.d0
          return
        endif
c
c 4 knot (3 layer) upper mantle, PREM crust radius
c
      elseif (iopt.eq.13) then
        if (r.le.r670) return
        if (r.le.rmoho.and.r.gt.r220) then
          u = (r - r220) / (rmoho - r220)
          func(1) = u
          func(2) = 1.d0 - u
        endif
        if (r.le.r220.and.r.gt.r400) then
          u = (r - r400) / (r220 - r400)
          func(2) = u
          func(3) = 1.d0 - u
        endif
        if (r.le.r400.and.r.gt.r670) then
          u = (r - r670) / (r400 - r670)
          func(3) = u
          func(4) = 1.d0 - u
        endif
        return
c
c Chebyshev polynomial in upper mantle
c
      elseif (iopt.eq.14) then
        if (r.gt.rmoho.or.r.lt.r670) return
        u = (r+r-rmoho-r670)/(rmoho-r670)
        goto 300
c
c Chebyshev polynomial in lower mantle
c
      elseif (iopt.eq.15) then
        if (r.gt.r670.or.r.lt.rcmb) return
        u = (r+r-r670-rcmb)/(r670-rcmb)
        goto 300
c
c Chebyshev polynomial for whole mantle
c
      elseif (iopt.eq.16) then
        if (r.gt.rmoho.or.r.lt.rcmb) return
        u = (r+r-rmoho-rcmb)/(rmoho-rcmb)
        goto 300
c
c 14 layer mantle
c
      elseif (iopt.eq.17) then
        if (r.gt.rmoho.or.r.lt.rcmb) return
        do 90 i=1,14
        j = i
  90    if (r.gt.rlmod(i+1)) goto 95
        return
  95    func(j) = 1.d0
        return
c
c 14 layer UPPER mantle
c
      elseif (iopt.eq.18) then
        if (r.gt.6371.d0.or.r.lt.5871.d0) return
        do 100 i=1,14
          j = i
  100     if (r.gt.rz(i+1)) goto 105
        return
  105   func(j) = 1.d0
        return
c
c B-spline, upper mantle, 9 knots
c
      elseif (iopt.eq.19) then
        if (r.gt.rmoho.or.r.lt.r670) return
           call bsrfun_um(r,func)
        return
c
c B-spline, lower mantle, 12 knots
c
      elseif (iopt.eq.20) then
        if (r.gt.r670.or.r.lt.rcmb) return
           call bsrfun_lm(r,func)
        return
c
c B-spline, whole mantle, 21 knots
c
      elseif (iopt.eq.21) then
        if (r.gt.rmoho.or.r.lt.rcmb) return
           call bsrfun_wm(r,func)
        return
      endif
c
c only get here if iopt = 1, 2 or 12 (i.e. using Legendre polynomials)
c
      usq=u*u
      func(1) = con(1)
      func(2) = con(2)*u
      func(3) = con(3)*(3.d0*usq-1.d0)
      func(4) = con(4)*u*(5.d0*usq-3.d0)
      func(5) = con(5)*(usq*(35.d0*usq-30.d0)+3.d0)
      func(6) = con(6)*u*(usq*(63.d0*usq-70.d0)+15.d0)
      func(7) = con(7)*(usq*(usq*(231.d0*usq-315.d0)+105.d0)-5.d0)
      func(8) = con(8)*u*(usq*(usq*(429.d0*usq-693.d0)+315.d0)-35.d0)
      return
c
c compute Chebyshev polynomials (iopt=14 or 15)
c
 300  u2 = 2.d0 * u
      func(1) = 1.d0
      func(2) = u
      do 310 i=3,14
  310 func(i) =  u2 * func(i-1) - func(i-2)
      do 320 i=1,14
  320 func(i) = conc(i) * func(i)
c
c explicit calculation of Chebyshev polynomials      
c      
c 300  usq     = u * u
c      func(1) = conc(1)
c      func(2) = conc(2) * u
c      func(3) = conc(3)*(2.d0*usq-1.d0)
c      func(4) = conc(4)*u*(4.d0*usq-3.d0)
c      func(5) = conc(5)*(usq*(8.d0*usq-8.d0)+1.d0)
c      func(6) = conc(6)*u*(usq*(16.d0*usq-20.d0)+5.d0)
c      func(7) = conc(7)*(usq*(usq*(32.d0*usq-48.d0)+18.d0)-1.d0)
c      func(8) = conc(8)*u*(usq*(usq*(64.d0*usq-112.d0)+56.d0)-7.d0)
c      func(9) = conc(9)*(usq*(usq*(usq*(128.d0*usq-256.d0)+160.d0)
c     &                   -32.d0)+1.d0)
c      func(10) = conc(10)*u*(usq*(usq*(usq*
c     &           (256.d0*usq-576.d0)+432.d0)-120.d0)+9.d0)
c      func(11) = conc(11)*(usq*(usq*(usq*(usq*
c     &          (512*usq-1280.d0)+1120.d0)-400.d0)+50.d0)-1.d0)
c      func(12) = conc(12)*u*(usq*(usq*(usq*(usq*
c     &           (1024.d0*usq-2816.d0)+2816.d0)-1232.d0)+220.d0)-11.d0)
c      func(13) = conc(13)*(usq*(usq*(usq*(usq*(usq*
c     &  (2048.d0*usq-6144.d0)+6912.d0)-3584.d0)+840.d0)-72.d0)+1.d0)
c      func(14) = conc(14)*u*(usq*(usq*(usq*(usq*(usq*
c     &(4096.d0*usq-13312.d0)+16640.d0)-9984.d0)+2912.d0)-364.d0)+13.d0)
c
      return
      end

	subroutine bsrfun_um(r,func)
	parameter (maxknt=30)
	dimension xknt(-1:maxknt+2)
	dimension b(4)
	real*8 r,func(1)		! func is initialized in call routine
	data isetknt/1/
	save knt,xknt
	
	if(knt.gt.maxknt)stop "too many knots in bsrfun_um"

	if(isetknt.eq.1)then	! set up um knots (only once)
	   call umknots(xknt,knt,maxknt)
	   isetknt=0
	endif
	
	x=6371.0-r		! knots are described as depths
	
	do ib=1,knt
	  do i=1,4
	     if((ib+i-3).ge.1.and.(ib+i-2).le.knt)then
	        if(x.ge.xknt(ib+i-3).and.x.le.xknt(ib+i-2))then
                  call bsplker(b,x,xknt,knt,ib+i-3,1)
                  func(ib)=b(5-i)
                endif
             endif
          enddo
        enddo

        return
        end

	subroutine bsrfun_lm(r,func)
	parameter (maxknt=30)
	dimension xknt(-1:maxknt+2)
	dimension b(4)
	real*8 r,func(1)		! func is initialized in call routine
	data isetknt/1/
	save knt,xknt
	
	if(knt.gt.maxknt)stop "too many knots in bsrfun_lm"

	if(isetknt.eq.1)then	! set up um knots (only once)
	   call lmknots(xknt,knt,maxknt)
	   isetknt=0
	endif
	
	x=6371.0-r		! knots are described as depths
	
	do ib=1,knt
	  do i=1,4
	     if((ib+i-3).ge.1.and.(ib+i-2).le.knt)then
	         if(x.ge.xknt(ib+i-3).and.x.le.xknt(ib+i-2))then
                  call bsplker(b,x,xknt,knt,ib+i-3,2)
                  func(ib)=b(5-i)
                 endif
             endif
          enddo
        enddo

        return
        end

	subroutine bsrfun_wm(r,func)
	parameter (maxknt=14)
	dimension xknt(-1:maxknt+2)
	dimension b(4)
	real*8 r,func(1)		! func is initialized in call routine
	data isetknt/1/
	save knt,xknt
	
	if(knt.gt.maxknt)stop "too many knots in bsrfun_wm"

	if(isetknt.eq.1)then	! set up wm knots (only once)
	   call wmknots(xknt,knt,maxknt)
	   isetknt=0
	endif
	
	x=6371.0-r		! knots are described as depths
	
	do ib=1,knt
	  do i=1,4
	     if((ib+i-3).ge.1.and.(ib+i-2).le.knt)then
	         if(x.ge.xknt(ib+i-3).and.x.le.xknt(ib+i-2))then
                  call bsplker(b,x,xknt,knt,ib+i-3,3)
                  func(ib)=b(5-i)
               endif
           endif
        enddo
      enddo

      return
      end

	subroutine bsplker(b,x,xknt,knt,ib,iumlm)
	dimension b(1),xknt(-1:knt+2)
	data ib0/0/,iumlm0/0/
	save ib0,a1,a2,a3,a4
	
	if(ib.ne.ib0.or.iumlm0.ne.iumlm)then	! diff. region, upper or lower mantle
	   temp=1.0/(xknt(ib+1)-xknt(ib))
	   a1=temp/(xknt(ib+1)-xknt(ib-1))
	   a3=temp/(xknt(ib+2)-xknt(ib))
         a4=a3/(xknt(ib+3)-xknt(ib))
	   temp=xknt(ib+2)-xknt(ib-1)
	   a2=a1/temp
	   a3=a3/temp
	   a1=a1/(xknt(ib+1)-xknt(ib-2))
	   ib0=ib
	   iumlm0=iumlm
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
c	if(ib.eq.(knt-1))then
c	   b(3)=b(3)+b(4)
c	endif
c
c...for 2nd derivative = 0
c
	if(ib.eq.1)then
	  b(3)=b(3)+b(2)/3.0
	  b(2)=b(2)/1.5+b(1)
	endif
	if(ib.eq.2)then
	  b(2)=b(2)+b(1)/3.0
	  b(1)=b(1)/1.5
	endif
	if(ib.eq.(knt-1))then
	   b(2)=b(2)+b(3)/3.0
	   b(3)=b(3)/1.5+b(4)
	endif
	
	if(ib.eq.(knt-2))then
	   b(3)=b(3)+b(4)/3.0
	   b(4)=b(4)/1.5
	endif

	return
	end

C Upper mantle knots for a split B-Spline model.
C The number of knot value in the upper mantle is 6.
	subroutine umknots(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	
	knt=6 
	print *, ' set upper mantle knots to ', knt
	xknt(1)=24.4
	xknt(2)=100.0
	xknt(3)=225.0
	xknt(4)=350.0
	xknt(5)=500.0
	xknt(6)=670.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end


C Lower mantle knots for a split B-Spline model.
C The number of radial knots in the lower mantle is 8.
	subroutine lmknots(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	
	knt=8 
	print *, ' set lower mantle knots to ', knt
	xknt(1)=670.0
	xknt(2)=820.0
	xknt(3)=1320.0
	xknt(4)=1820.0
	xknt(5)=2320.0
	xknt(6)=2550.0
	xknt(7)=2791.0
	xknt(8)=2891.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)

	return
	end


	subroutine wmknots(xknt,knt,nmax)
	dimension xknt(-1:nmax)
	
	print *, ' set whole mantle knots'
	
	knt=14  
	xknt(1)=24.4
	xknt(2)=75.0
	xknt(3)=125.0
	xknt(4)=225.0
	xknt(5)=350.0
	xknt(6)=500.0
	xknt(7)=670.0
	xknt(8)=820.0
	xknt(9)=1320.0
	xknt(10)=1820.0
	xknt(11)=2320.0
	xknt(12)=2550.0
	xknt(13)=2791.0
	xknt(14)=2891.0
  
        xknt(0)=xknt(1)
        xknt(-1)=xknt(0)
        xknt(knt+1)=xknt(knt)
        xknt(knt+2)=xknt(knt)
  
	return
	end
