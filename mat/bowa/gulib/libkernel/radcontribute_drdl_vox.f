      subroutine radcontribute_drdl_vox(xleng,xlseg,nseg,qrad_l,qth,xtu,ierr)	
c     +		rknt,nknt,numsplit,isplayer,xtu)
c
c This routine integrates over dlength using  Gauss-Legendre 5-point integration
c
c input: th, phi     --- geocentric locations of splines
c 	   h           --- size of the spherical splines
c 	   xlseg       --- incremental epicentral distance
c	   xleng       --- distance along ray
c	   nseg        --- total number of points along segment
c 	   rknt, nknt  --- radius and number of radial knots
c	   numsplit    --- number of spliting depth in model (e.g., 1 if split at 670)
c	   isplayer    --- 2-D array to save spliting layers
c	   		   e.g., U6L8 will be saved as 
c			   isplayer(1,1)=1 (layer 1 bottom index)
c			   isplayer(1,2)=8 (layer 1 top index)
c			   isplayer(2,1)=9 (layer 2 bottom index)
c			   isplayer(2,2)=14 (layer 2 top index)
c	   xtu         --- turning point of ray
c output:
c	   arowk       --- array of radial contributions for the present horizontal
c			   spline
c Notes: slight inacuracies near interface at Moho.
c        modified by lapo boschi from radcontribute_drdl by yu gu
c
c
      implicit double precision (a-h,o-z)

c  basic parameters, naming is obsolete
      parameter(maxfreq=8)
      parameter(maxnode=maxfreq*maxfreq*10+2)
      parameter(maxrknot=50)
      parameter(maxparm=maxnode*maxrknot)
      data pi /3.141592653579/
      data rad/57.2957795130824/
      data reprad /57.295779579/    

      common/amatrix/arow(maxparm),indarow(maxparm),nonzero_a
      common/invopt/numinvvar,invparlist(5)
      common/radcontrib/qseg(5,50000),arowk(5,maxrknot) !la 3200->50000
      common/anisopath/theta_path(50000), rad_path(50000),ntheta
      dimension xlseg(1),xleng(1),qrad_l(3,1),qth(3,1)
      dimension rknt(1), isplayer(maxrknot,1) ! arowk already defined in common - lapo 22.1.2010
      dimension b1(4)
      dimension x(5),w(5)
      dimension valp(5),qp(5),sum(5)
      data x/.1488743389d0,.4333953941d0,.6794095682d0,
     *       .8650633666d0,.9739065285d0/
      data w/.2955242247d0,.2692667193d0,.2190863625d0,
     *       .1494513491d0,.0666713443d0/


      ierr=0
      do i=1,5
	arowk(i,1)=0.d0
	sum(i)=0.d0
	qp(i)=0.d0
      enddo


c... Gauss-Legendre integration over the running angle.
c... I need: incremental length along ray (xlseg) starting from entry point into voxel
c... and radius and delta along ray


      do i=1,nseg-1 ! Loop over all (nptint=10) segs of ray within voxel
	   xr=xlseg(i+1)-xlseg(i) ! Compute length (km) of current segment
	   do j=1,5 ! LOOP OVER WHAT? 


c... find theta values at pole locations of weighting function
		  dx=xr*x(j) ! WHAT IS DONE HERE? SOME TYPE OF WEIGHTING?
		  xmpdx=xlseg(i)+dx

c... interpolate r with respect to xleng using cubic splines
c... will not work near receiver because of error in yu's definitions
		  rp=drsple(1,ntheta,xleng,rad_path,qrad_l,xmpdx) 

c... Some exceptions and escape sequences
		  if(abs(rp-xtu).lt.0.01) rp=xtu+0.01  ! ensures the singularity is avoided
      	  if(rp.gt.6371.)then ! escape when rp is too large
	           ierr=1
	           return
	        endif

c... extract theta value at current location (segment)
		  tp=drsple(1,ntheta,xleng,theta_path,qth,xmpdx)

c... qp will contain have kernels for Vph, Vpv, Vsh, Vsv and eta respectively
c... note that calcqvec implicitly also uses p (see tder.f)!
		  call calcqvec(rp,qp)

c... Gu: valp,valm = (dr/dzeta)*v0, drdl=(dr/ds) 
c... Lapo is not sure what this means and assumes valp=dr/ds)
		  call calcdrdl(rp,valp,1,tp,viso)
	
		  do ipar=1,numinvvar ! Loop over actually selected inversion parameters

			ind=invparlist(ipar)
			qp(ind)=qp(ind)*valp(ind)*xr

c------------------------------------------------------------------------
c ... By setting B(i)*S(j) to 1, we effectively assume a dv/v=1.0, then
c ... the sum of the int(dq/dv1+dq/dv2)=-(travel time) if this works well.
c ... To test that, 
c ... the array sum here can be activated to compare travel time with sum 
c ... of kernels, YG, 2002.
c
c			sum(ind)=sum(ind)+w(j)*qp(ind)
c------------------------------------------------------------------------

			arowk(ind,1)=arowk(ind,1)+w(j)*qp(ind) ! save row for current measurement

		  enddo
	   enddo
      enddo

c ... Activate the following print for testing (make sure to uncomment the 
c ... corresponding test codes in voxel_int.f
c     print*, 'The sum of S kernel =', sum(3)+sum(4)
c-----------------------------------------------------------------------

      return
      end
