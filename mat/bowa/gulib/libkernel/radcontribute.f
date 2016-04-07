
      subroutine radcontribute(th,phi,h,tseg,rseg,nseg,qrad,
     +		rknt,nknt,xtu)
c
c compute radial contribution and multiply to the lateral contribution of
c splines.  Note fdelta computes horizontal spline values for a given distance
c away from spline center, and bsplinefun computes radial spline function value.
c This routine integrates over dlength using  Gauss-Legendre 10-point integration/
c Note:  I suggest not to use this routine since the accuracy is slightly
c  	 off (not bad for S and ScS, but 8 seconds or so off of SS).  This is 
c	 mainly due to the instability of dr/dzeta which has a singularity problem.
c	 But for the purpose of inversion, this is fine.
c	 A better version:  radcontribute_drdl.f (integrate along ds instead of dzeta)
c	 J.Gu, 2002.
c
c input:   th, phi     --- geocentric locations of splines
c 	   h           --- size of the spherical splines
c 	   tseg,rseg   --- running angle, radius along segment
c	   nseg        --- total number of points along segment
c 	   rknt, nknt  --- radius and number of radial knots
c	   xtu         --- turning point of ray
c output:
c	   arowk       --- array of radial contributions for the present horizontal
c			   spline
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


      common/amatrix/arow(maxparm),indarow(maxparm),nonzero_a
      common/invopt/numinvvar,invparlist(5)
      common/radcontrib/qseg(5,3200),arowk(5,maxrknot)
      dimension tseg(1),rseg(1),qrad(3,1)
c      dimension rknt(1), arowk(5,maxrknot)
      dimension rknt(1) ! arowk already defined in common block - lapo 22.1.2010
      dimension func(3200),b0(4),b1(4)
      dimension x(5),w(5),dxp(10),dxm(10),dxpf(10),dxmf(10)
      dimension valp(5), valm(5), qp(5), qm(5), sum(5)
      data x/.1488743389d0,.4333953941d0,.6794095682d0,
     *       .8650633666d0,.9739065285d0/
      data w/.2955242247d0,.2692667193d0,.2190863625d0,
     *       .1494513491d0,.0666713443d0/

	
      ierr=0
      do i=1, 5
      	do k=0, nknt
		arowk(i,k)=0.d0
	enddo
	sum(i)=0.d0
	qp(i)=0.d0
	qm(i)=0.d0
      enddo

c
c...Gauss-Legendre integration over the running angle.
c
      do i=1, nseg-1
	 write(12, *) tseg(i)*reprad, rseg(i)
	 xm = 0.5d0*(tseg(i+1)+tseg(i))
	 xr = 0.5d0*(tseg(i+1)-tseg(i))
	 do j=1, 5
c... find theta values at pole locations of weighting function
		dx = xr*x(j)
		xmpdx = xm + dx
		xmmdx = xm - dx
c... interpolate r use cubic splines
c		rp=drsple(1,nseg,tseg,rseg,qrad,xmpdx)
c		rm=drsple(1,nseg,tseg,rseg,qrad,xmmdx)
		rp=drsple(i,i+1,tseg,rseg,qrad,xmpdx)
		rm=drsple(i,i+1,tseg,rseg,qrad,xmmdx)
c... compute horizontal splines
	 	dist=dist_phi(th,phi,xmpdx)
	 	bsp=fdelta(dist,h)
	 	dist=dist_phi(th,phi,xmmdx)
	 	bsm=fdelta(dist,h)
		write(11, *) xmmdx*reprad, rm
c
c... compute radial contribution from splines
		call getbsreg(rm,rknt,nknt,iregm)
		call getbsreg(rp,rknt,nknt,iregp)
      		call bsplinefun(iregm,rm,rknt,nknt,b0)
      		call bsplinefun(iregp,rp,rknt,nknt,b1)
c... below ensures that the bottoming singularity is avoided
c
		if(abs(rp-xtu).lt.0.01) rp=xtu+0.01
		if(abs(rm-xtu).lt.0.01) rm=xtu+0.01
c
c... compute dq/dv for the selected velocities, the orders are
c... Vph, Vpv, Vsh, Vsv and eta.
		call calcqvec(rp,qp)
		call calcqvec(rm,qm)
c
c... error checking to make sure that singularity does not 
c    result in unresonable integration.
		call drdzeta(rp,valp,1,xmpdx,viso)
		rpviso=viso
		call drdzeta(rm,valm,1,xmmdx,viso)
		rmviso=viso
		write(13, *) xmmdx*reprad, valm(3)
		write(14, *) xmmdx*reprad, valm(4)
c
c ...Instead of using absolute perturbation for inversion the option ``1'' uses
c ...dv/v0 or deta/eta0 as a model parameters.  valp and valm contains 
c ...(dr/dzeta)*v0 terms.     drdzeta = (dr/dzeta)
		do ipar=1, numinvvar
			ind=invparlist(ipar)
c ...By setting B(i)*S(j) to 1, we effectively assume a dv/v=1.0, then
c ...the sum of the int(dq/dv1+dq/dv2)=-(travel time) if this works well.
c
			qp(ind) = qp(ind)*valp(ind)*1.0*xr
			qm(ind) = qm(ind)*valm(ind)*1.0*xr

c			qp(ind) = qp(ind)*valp(ind)*bsp*xr
c			qm(ind) = qm(ind)*valm(ind)*bsm*xr
c			print*, qp(ipar),valp(ind), xr
			iregmm2=iregm-2
			iregpm2=iregp-2
			do jb=1, 4
			   km=iregmm2+jb
			   if(km.gt.0.and.km.le.nknt) then
				 arowk(ipar,km)=arowk(ipar,km)+w(j)*(qp(ind)+qm(ind))*b0(jb)
				print*, 'arowk=', arowk(ipar,km)
			   endif
			   kp=iregpm2+jb
			   if(kp.gt.0.and.kp.le.nknt) then
				 arowk(ipar,kp)=arowk(ipar,kp)+w(j)*(qp(ind)+qm(ind))*b1(jb)		
			   endif
			enddo
			sum(ind)=sum(ind)+w(j)*(qp(ind)+qm(ind))
		enddo
10		continue
	 enddo
      enddo
	print*, 'The travel time =', sum(3)+sum(4)
	stop
      return
      end
