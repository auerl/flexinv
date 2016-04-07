      subroutine radcontribute_drdl(th,phi,h,xleng,xlseg,nseg,qrad_l,qth,
     +		rknt,nknt,numsplit,isplayer,xtu)
c
c compute radial contribution and multiply to the lateral contribution of
c splines.  Note fdelta computes horizontal spline values for a given distance
c away from spline center, and bsplinefun computes radial spline function value.
c This routine integrates over dlength using  Gauss-Legendre 5-point integration
c
c input:   th, phi     --- geocentric locations of splines
c 	   h           --- size of the spherical splines
c 	   xlseg       --- keep track of absolute length of ray
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
c Note: slight inacuracies near interface at Moho.
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


      common/amatrix/arow(maxparm),indarow(maxparm),nonzero_a
      common/invopt/numinvvar,invparlist(5)
      common/radcontrib/qseg(5,3200),arowk(5,maxrknot)
      common/anisopath/theta_path(3200), rad_path(3200),ntheta
      dimension xlseg(1),xleng(1),qrad_l(3,1),qth(3,1)
c      dimension rknt(1), arowk(5,maxrknot), isplayer(maxrknot,1)
      dimension rknt(1), isplayer(maxrknot,1) ! arowk already defined in common - lapo 22.1.2010
      dimension b1(4)
      dimension x(5),w(5)
      dimension valp(5),qp(5),sum(5)
      data x/.1488743389d0,.4333953941d0,.6794095682d0,
     *       .8650633666d0,.9739065285d0/
      data w/.2955242247d0,.2692667193d0,.2190863625d0,
     *       .1494513491d0,.0666713443d0/
      data rmoho/6346.62891/
      data rcmb/3479.96/

	
c	print*, 'what goes on??'
      ierr=0
      do i=1, 5
      	do k=1, nknt
		arowk(i,k)=0.d0
	enddo
	sum(i)=0.d0
	qp(i)=0.d0
      enddo
c
c...Gauss-Legendre integration over the running angle.
c
      do i=1, nseg-1
	 xr = xlseg(i+1)-xlseg(i)
cTEST
c	print*,i,xlseg(i)

	 do j=1, 5
c... find theta values at pole locations of weighting function
		dx=xr*x(j)
		xmpdx=xlseg(i)+dx
c... interpolate r using cubic splines
		rp=drsple(1,ntheta,xleng,rad_path,qrad_l,xmpdx)
		if(abs(rp-xtu).lt.0.01) rp=xtu+0.01  ! ensures the singularity is avoided
		tp=drsple(1,ntheta,xleng,theta_path,qth,xmpdx)
c... compute horizontal splines
	 	dist=dist_phi(th,phi,tp)
	 	bsp=fdelta(dist,h)

c... compute radial contribution from splines
		if(numsplit.eq.0) then
			call getbsreg(rp,rknt,1,nknt,iregp)  !continuous model
      			call bsplinefun(iregp,rp,rknt,1,nknt,b1)
		else
			do is=1, numsplit+1
c... loop through radial layers
			   rbot=rknt(isplayer(is,1))
			   rtop=rknt(isplayer(is,2))
			   if(rp.gt.rbot.and.rp.le.rtop) then
				call getbsreg(rp,rknt,isplayer(is,1),isplayer(is,2),iregp)
				goto 5
			   endif
			enddo
5      			call bsplinefun(iregp,rp,rknt,isplayer(is,1),isplayer(is,2),b1)
		endif
c
c... compute dq/dv for the selected velocities, the orders are
c... Vph, Vpv, Vsh, Vsv and eta.
		call calcqvec(rp,qp)
cTEST
	write(99,"(6(e14.7,1x))")rp,(qp(k),k=1,5)
c
c ...Instead of using absolute perturbation for inversion the option ``1'' uses
c ...dv/v0 or deta/eta0 as a model parameters.  valp and valm contains 
c ...(dr/dzeta)*v0 terms.     drdl=(dr/ds)
		call calcdrdl(rp,valp,1,tp,viso)

		do ipar=1, numinvvar
			ind=invparlist(ipar)
			qp(ind)=qp(ind)*valp(ind)*xr
c------------------------------------------------------------------------
c ...By setting B(i)*S(j) to 1, we effectively assume a dv/v=1.0, then
c ...the sum of the int(dq/dv1+dq/dv2)=-(travel time) if this works well.
c ...To test that, 
c ...the array sum here can be activated to compare travel time with sum 
c... of kernels, YG, 2002.
c
c			sum(ind)=sum(ind)+w(j)*qp(ind)
cTEST
c			write(199,"(6(e14.7,1x))")rp,(sum(k),k=1,5) ! lapo
c------------------------------------------------------------------------
c
c... multiply by horizontal function value
			qp(ind)=qp(ind)*bsp
			iregpm2=iregp+isplayer(is,1)-3
			do jb=1, 4
			   kp=iregpm2+jb
			   if(kp.ge.isplayer(is,1).and.kp.le.isplayer(is,2)) then
				 arowk(ind,kp)=arowk(ind,kp)+w(j)*qp(ind)*b1(jb)
			   endif
			enddo
		enddo
10		continue
	 enddo
      enddo
c   Activate the following print for testing (make sure to uncomment the corresponding
c   test codes in spline_int.f------------------------------------------
c      print*, 'The sum of S kernel =', sum(3)+sum(4)
c-----------------------------------------------------------------------
cTEST
c	pause
      return
      end
