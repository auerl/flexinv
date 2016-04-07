	subroutine getkIn_int(iradfun,qray,lmax,kmax,xtu)
c  The follownig routine uses Gauss-Legendre 10-point integration
c  to evaluate integrals of k_I_n of Dziewonski, 1984.
c  Also see Wei-Jia Su's thesis for reference, YG, 2000.
c  *** NOTE:*************************************************************
c       This program is only valid for mantle inversion only,
c       and it is not designed for phases like SP, PS etc.
c       For phases such as SKS, the calculation of dr/dzeta is
c       not reasonable.  But since that part of the kernel is
c       not used for mantle inversion, it is ok.
c ************************************************************************
c  Input:	qray --- model derivative of qtau, this is g kernel
c	     iradfun --- choice of radial function
c		lmax --- maximum spherical harmonic (or surface B-splines)
c		kmax --- radial order (number of radial B-spline)
c		xtu  --- turning point of ray
c  Output:
c		sum  --- kIn kernel saved in /kernel/ 
c
c  
	implicit double precision(a-h,o-z)
	integer kmax, iradfun,lmax
	dimension qrad(3, 50000), fwork(3,50000)
	dimension funcp(15), funcm(15)
	dimension x(5),w(5),dxp(10), dxm(10), dxpf(10), dxmf(10)
	dimension valp(10), valm(10)
	dimension qray(4,5,20)
        data x/.1488743389d0,.4333953941d0,.6794095682d0,
     *       .8650633666d0,.9739065285d0/
        data w/.2955242247d0,.2692667193d0,.2190863625d0,
     *       .1494513491d0,.0666713443d0/
        parameter (rmoho=6371.d0-24.4d0)
        parameter (rcmb=3480.d0)
        common/coeff/coef(4,8,20)
	common/path$/delt(2,1600,2),nfin(2),kdep
	common/layr/nl(20),xb(20),xt(20),ifanis,nplay
	common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
	common/ray$/pmin,pmax,nps(20,2),npsrc(2),nn,iup,jref,npleg,nsleg
	common/anisopath/theta_path(3200),rad_path(3200),ntheta
        common/kernel/sum(362,20,10),ifsplit,ifdiff
	common/premvel/vphp,vpvp,vshp,vsvp,etap,rhop,vpiso,vsiso
        common/rayinfo/theta_bt(20),rad_bt(20),nbot,ithbot,thetamax,
     &     rayseg_th(20),nrayseg_type(20),nseg,delreq
	common/xnorm/xnorm
        common/parm1/xa,xc,xn,xl,xf,tdif(5)
	common/isot/isotrp
	ttiso = 0.d0 	! isotropic (so called) travel time to test integration
c============================real code after test =============================
c initialize sum() and 
c spherical Harmonics
c

	sumker=0.d0
	write(*,*) ' ------------ray segments---------------------'
	do i=1, ntheta
		write(101,*) theta_path(i)*180.0/3.141592, rad_path(i)
	enddo
c
c initialization
c
	if(iradfun.ge.14.and.iradfun.le.18.or.iradfun.eq.22) then
		do i=1, 5 
			do j=1, kmax+1
				do k=1, 2*lmax+1
					sum(k,j,i) = 0.d0
				enddo
			enddo
		enddo
	else if(iradfun.gt.18.and.iradfun.lt.22) then
c
c B-splines
c
		do i=1, 5
			do j=1, kmax+1
				do k=1, lmax
					sum(k,j,i) = 0.d0
				enddo
			enddo
		enddo
	endif
c----------------------------------------------------------------------
	rad_path(ntheta+1) = rad_path(ntheta)
	theta_path(ntheta+1) = theta_path(ntheta)
	
	call drspln(1, ntheta, theta_path, rad_path, qrad, fwork)  ! cubic spline
	do i=1, ntheta-1
		write(99,*) theta_path(i)*180.0/3.14159265, 6371.-rad_path(i)
		xm = 0.5d0*(theta_path(i+1)+theta_path(i))
      		xr = 0.5d0*(theta_path(i+1)-theta_path(i))
		do j=1, 5
c find corresponding theta values to pole locations of weighting function
			dx = xr*x(j)
			xmpdx = xm + dx
			xmmdx = xm - dx
          		rp=drsple(1,ntheta,theta_path,rad_path,qrad,xmpdx)
          		rm=drsple(1,ntheta,theta_path,rad_path,qrad,xmmdx)
c
c... below to ensure to avoid the sigularity-----------------------------
			if(abs(rp-xtu).lt.0.01) rp=xtu+0.01
			if(abs(rm-xtu).lt.0.01) rm=xtu+0.01
c------------------------------------------------------------------------
      			do ilay=1, numlyr
	 			if(rp.gt.xb(ilay).and.rp.le.xt(ilay)) then
					n1 = ilay
					goto 1011
				endif
      			enddo
1011			continue
      			do ilay=1, numlyr
	 			if(rm.gt.xb(ilay).and.rm.le.xt(ilay)) then
					n2 = ilay
					goto 1012
				endif
      			enddo
1012			continue
c
c  --------qq1, qq2 contain qtau for rm1 and rm2-----------------------
c
     			call qtau(rp,n1,q1)
			if(q1.ge.0.0) then
c
c  if SV or SH reflect at fluid boundary, then q2=-1.d0
c
      				qq1=dsqrt(q1)
      				call tder(rp,qq1)
      				do ipar=1, 5
					dxp(ipar)=tdif(ipar)
      				enddo
			else
      				do ipar=1, 5
					dxp(ipar)=0.d0
      				enddo				
			endif
     			call qtau(rm,n2,q2)
			if(q2.ge.0.0) then
      				qq2=dsqrt(q2)
      				call tder(rm,qq2)
      				do ipar=1, 5
					dxm(ipar)=tdif(ipar)
      				enddo
			else
      				do ipar=1, 5
					dxm(ipar)=0.d0
      				enddo				
			endif
			write(102, *) xmpdx, tdif(3)
			write(103, *) xmpdx, tdif(4)
c----------------------------------------------------------------------
c
c  the following finds (dr/dzeta)*M0, if 1
c
			call drdzeta(rp,valp,1,xmpdx,vviso)
			rpviso = vviso
			call drdzeta(rm,valm,1,xmmdx,vviso)
			rmviso = vviso
c
c Instead of using absolute perturbation for inversion this uses
c dv/v0 or deta/eta0 as a model parameters.  valp and valm contains (dr/dzeta)*v0 terms.
c drdzeta = (dr/dzeta)
			do ip=1, 5
			  dxp(ip) = dxp(ip)*valp(ip)*xr
			  dxm(ip) = dxm(ip)*valm(ip)*xr
			enddo
c
c Chebyshev polynomial (whole mantle) and spherical harmonics (horizontal)
c			
   5			if(iradfun.ge.14.and.iradfun.le.18) then
			    call radfun(iradfun, funcp, rp)
			    call radfun(iradfun, funcm, rm)
			else if(iradfun.eq.22) then
c Upper mantle and Lower mantle split Chebyshev model (U7L5)
			    call radfun_u7l5(rp,funcp)
			    call radfun_u7l5(rm,funcm)
			else
			    write(*,*) 'Not implemented yet, quit...'
			    stop
			endif
c
c----------------------------------------------------------------------
c uses sh velocity for the ray to calculate the travel times
c the equation is r*r/(v*v*p)
			if(rpviso.gt.0.0001) then
				ttisop = xr*rp*rp/(rpviso*rpviso*p)
			else
				ttisop=0.d0
			endif
			if(rmviso.gt.0.0001) then
				ttisom = xr*rm*rm/(rmviso*rmviso*p)
			else
				ttisom=0.d0
			endif
c			print*, rm, rmviso, rpviso
			ttiso=ttiso+w(j)*(ttisop+ttisom)
c
c----------------------------------------------------------------------
c
			do m = 1, lmax+1
      				ym=m-1
      				cp=dcos(ym*xmpdx)
     		 		cm=dcos(ym*xmmdx)
      				sp=dsin(ym*xmpdx)
      				sm=dsin(ym*xmmdx)
			        do k=1, kmax+1
c ****loop over 5 anisotropic parameters
				    do ip=1, 5
				       dxpf(ip)=dxp(ip)*funcp(k)
			       	       dxmf(ip)=dxm(ip)*funcm(k)
				       if(m.eq.1) then
                     			  sum(1,k,ip)=sum(1,k,ip)+w(j)*(dxpf(ip)+dxmf(ip))
				       else
c ****integration in the order of Vph, Vpv, Vsh, Vsv, eta
                     			  sum(2*m-2,k,ip)=sum(2*m-2,k,ip)+w(j)*(dxpf(ip)*cp+dxmf(ip)*cm)
                     		          sum(2*m-1,k,ip)=sum(2*m-1,k,ip)+w(j)*(dxpf(ip)*sp+dxmf(ip)*sm)
				       endif
			            enddo
				enddo	
			enddo
		enddo
c				if(rp.lt.3550.or.rm.lt.3550) then
c				   do jjj=1, 14
c					write(*,*) jjj, 'rp=', rp, sum(1,1,3),sum(1,1,4)
c				   enddo
c				endif
c		print*, sum(1,1,1),sum(1,1,2),sum(1,1,3),sum(1,1,4)
c		print*, valm(1), valm(2), valm(3), valm(4)
c		write(33, *) 6371.0-rp, sum(1,1,3)
c		write(23,*) 6371.0-rp, sum(1,1,4)  
c		write(24,*) xm*180.0/3.1415926, -ttiso
	enddo
c *** S wave
	write(*,*) 'integrated S derivative = ', sum(1,1,3)+sum(1,1,4)
	write(*,*) 'integrated P derivative = ', sum(1,1,1)+sum(1,1,2)
	write(*,*) 'Horizontally integrated travel time =', ttiso
	return
	end

