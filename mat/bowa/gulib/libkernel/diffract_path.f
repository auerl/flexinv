	subroutine diffract_path(th)
c  this subroutine will add the flat portion of the diffracted
c  wave into the original path, YG, 1998.
c  
	implicit real*8(a-h,o-z)
        common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
	common/anisopath/theta_path(50000),rad_path(50000),ntheta !la 3200->50000
	common/kernel/sum(362,20,10),ifsplit,ifdiff
        common/rayinfo/theta_bt(20),rad_bt(20),nbot,ithbot,thetamax,
     &     rayseg_th(20), nrayseg_type(20),nseg,delreq
	dimension  tempr(50000), tempt(50000) !la 1600->50000

	thdif = delreq*3.141592653579/180.d0-th ! diffracted distance
	write(*,*) 'diffracted distance =', thdif*180.d0/3.141592653579
	dth = 1.0*3.14159265/180.d0
	write(*,*) 'dth =', dth
	if(thdif.lt.(dth*0.5d0)) then
      		ndiff = 0
	else if(thdif.lt.(dth*2.d0)) then
		ndiff = 1
		dth = thdif/2.d0
	else
		ndiff = thdif/dth
		d1 = thdif-ndiff*dth
		if(d1.lt.dth*0.2d0) ndiff = ndiff-1			
	endif 
	write(*,*) 'ndiff =', ndiff, ' dth=', dth
	j=0
	do i=ithbot+1, ntheta
c use 0.5 degree spacing for diffracted portion
		j=j+1
		tempr(j) = rad_path(i)
		tempt(j) = theta_path(i)
	enddo 
	j=1
	df = 0.d0
	n1 = ithbot+1
	n2 = ithbot+ndiff
	do i=n1, n2
		df = df + dth
		rad_path(i) = rad_path(ithbot)
		theta_path(i) = df+theta_path(ithbot)
c		write(*,*) 'i=', i, theta_path(i)*180./3.1415927, rad_path(i)
	enddo
	do i=1, ntheta-ithbot
		rad_path(i+n2) = tempr(i)
		theta_path(i+n2) = tempt(i)+thdif
c		write(*,*) i+n2, theta_path(i+n2)*180./3.1415927, rad_path(i+n+2)
	enddo
	ntheta = ntheta + ndiff
	return
	end
