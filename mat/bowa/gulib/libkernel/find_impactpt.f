      subroutine find_impactpt(rdis,arytrans,idim1,ntrans,aryrefl,idim2,nrefl)
      implicit double precision (a-h,o-z)
c
c  This routine finds the points where a ray either transmit through
c  a discontinuity or reflect at a discontinuity.  It allows us to 
c  locate these impact points and compute surface kernels for topography
c  inversion using SS/PP precursors.
c  Input:
c	rdis ---  radius of a given discontinuity
c	raypath in /anisopath/
c  Output:
c	arytrans  --- array containing transmission points through a discontinuity
c	ntrans    --- number of transmission points
c	aryrefl   --- array containing reflection points at a discontinuity
c	nrefl     --- number of reflection points
c
c
      common/anisopath/theta_path(3200),rad_path(3200),ntheta
      dimension arytrans(1), aryrefl(1), tempary(20)

      do i=1, idim1
	 arytrans(i)=0.d0
	 tempary(i)=0.d0
      enddo
      do i=1, idim2
	 aryrefl(i)=0.d0
      enddo
      ntrans=0
      nrefl=0
      rlast=rad_path(1)
      rdifflast=0.d0
c
c
      do i=2, ntheta
	rdiff1=rad_path(i)-rdis
	rdiff2=rlast-rdis
	rdiffnew=rad_path(i)-rlast
	if(rdifflast.gt.0.d0.and.rdiffnew.lt.0.d0) then
c*** The following is a crude method to find where ray turns.  Then check
c*** if it turns near the given discontinuity, if so, precursor turning point!
		if(abs(rlast-rdis).lt.0.1) then
			nrefl=nrefl+1
			aryrefl(nrefl)=xlast
		endif
	endif 
c
c*** transmission through a discontinuity
	if(abs(rdiff1).lt.0.01.and.abs(rdiff2).gt.0.01) then
c*** if there is a ray point exactly at discontinuity, accounted once only
		ntrans=ntrans+1
		tempary(ntrans)=theta_path(i)
		
	else if((rdiff1.gt.0.d0.and.rdiff2.lt.0.d0).or.
     &		(rdiff1.lt.0.d0.and.rdiff2.gt.0.d0)) then
c*** if nothing close to discontinuity, linearly interpolate
		xx=(rdis-rlast)*(theta_path(i)-xlast)/(rad_path(i)-rlast)
		ntrans=ntrans+1
		tempary(ntrans)=xx
	endif
	rlast=rad_path(i)
	xlast=theta_path(i)
	rdifflast=rdiffnew
      enddo
c
c The following is integrity check, for SS precursors, the reflection point
c usually gets accounted also as a transmission point
c
      k=0
      do i=1, ntrans
	 do j=1, nrefl
c...  use .01 degree accuracy to pick them out
	 	if(abs(tempary(i)-aryrefl(j)).lt.0.00017) then
			goto 10
		endif
	 enddo
	 k=k+1
	 arytrans(k)=tempary(i)
10    enddo
      ntrans=k
c      do i=1, nrefl
c		print*,i,  '  reflection =', aryrefl(i)*180/3.1415
c      enddo
c      do i=1, ntrans
c		print*,i,  '  transmission =', arytrans(i)*180/3.1415
c      enddo
      return
      end
 
