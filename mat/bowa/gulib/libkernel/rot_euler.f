c  The subroutine rot_euler rotates coordinate frames based on the Euler
c  angles alpha, beta and gamma (see Edmonds page 7 and page 53).
c  The input theta and phi are the coordinates in the geographic
c  coordinate system. The angle zeta is measured counterclockwise
c  from due east. The returned values of theta, phi, and zeta are
c  in the equatorial coordinate system.
c  The returned value of theta is between 0 and pi, and the returned
c  value of phi is between 0 and 2*pi.
c  The inverse operation, i.e. from the rotated (equatorial)
c  frame back to the original (geographic) frame is performed
c  by replacing alpha by -gamma, beta by -beta, and gamma by -alpha
c  in the call to the subroutine rot_euler.

c The following are called by the rot_euler subroutine
	function icheck_if_even(i)
	n=(i-2*((i)/2))
	if(n.eq.0) then
		icheck_if_even=1
	else 
		icheck_if_even=0
	endif
	return
	end

	function icheck_if_odd(i)
	n=(i-2*((i)/2))
	if(n.eq.0) then
	 	icheck_if_odd=1
	else 
		icheck_if_odd=0
	endif
	return
	end

	subroutine collon2euler(scol,slon,rcol,rlon,delta,alpha,beta,gamma)
	implicit double precision (a-h,o-z)
 	parameter(PI=3.141592653579)
 	parameter(PI_2=PI/2.0)	
	tiny=1e-6
 	
	delta=acos(cos(scol)*cos(rcol)+sin(scol)*sin(rcol)*cos(rlon-slon))
	sina=sin(rlon-slon)*sin(rcol)/sin(delta)
	cosa=(cos(rcol)-cos(scol)*cos(delta))/(sin(scol)*sin(delta))
  	az=atan2(sina,cosa)
  	sina=sin(az)
	cosa=cos(az)
  	if(cosa.eq.0) then
      		if(cos(scol).gt.0) then
			eta=PI_2
      		else 
			eta=-PI_2
		endif
 	else
      		eta=atan(cos(scol)/(cosa*sin(scol)))
	endif
  	if(eta.lt.0.0) eta=eta+PI
  	gamma=PI_2+eta
  	mu=atan2(sina*sin(scol),cos(scol)/sin(eta))
  	beta=PI_2-mu
  	ep=atan2(sina*sin(eta),cos(eta)/sin(scol))
  	alpha=-PI_2+slon-ep
	return
	end

	subroutine reduce(pth, pph)
c force theta between 0 and pi, and phi between 0 and 2*pi
	implicit double precision (a-h,o-z)
 	parameter(PI=3.141592653579)
 	parameter(PI_2=PI/2.0)

	PI2=PI*2
  	th=pth
	ph=pph
  	i=abs(ifix(ph/PI2))
	if(ph.lt.0.0) then
		ph=ph+(i+1)*PI2
  	else
		if(ph.gt.PI2) ph=ph-i*PI2
	endif
  	pph=ph
  	if(th.lt.0.0.or.th.gt.PI) then
    		i=ifix(th/PI)
    		if(th.gt.0.0) then
      			if(icheck_if_odd(i).eq.1) then
        			th=(i+1)*PI-th
				if(ph.lt.PI) then
					ph=ph+PI
				else
					ph=ph-PI
				endif
      			else 
				th=th-i*PI
			endif
    		else
      			if(icheck_if_even(i).eq.1) then
        			th=-th+i*PI
				if(ph.lt.PI) then
					ph=ph+PI
				else
					ph=ph-PI
				endif
      			else
     		 		th=th-i*PI
			endif
  		endif
    		pth=th
		pph=ph
  	endif
	return
	end


	subroutine rot_euler(alpha,beta,gamma,theta,phi,zeta)
	implicit double precision (a-h,o-z)
  	dimension r_hat(0:3),rp_hat(0:3),k_hat(0:3),kp_hat(0:3),rotation(0:3,0:3)

  	sint=sin(theta)
	cost=cos(theta)
  	sinp=sin(phi)   
	cosp=cos(phi)
  	sina=sin(alpha)
	cosa=cos(alpha)
  	sinb=sin(beta)
	cosb=cos(beta)
  	sing=sin(gamma)
  	cosg=cos(gamma)
  	sinz=sin(zeta)
	cosz=cos(zeta)
	r_hat(0)=sint*cosp
	r_hat(1)=sint*sinp
	r_hat(2)=cost
	k_hat(0)=-(sinz*cost*cosp+cosz*sinp)
  	k_hat(1)=cosz*cosp-sinz*cost*sinp
  	k_hat(2)=sinz*sint
  	rotation(0,0)=cosg*cosb*cosa-sing*sina
  	rotation(0,1)=cosg*cosb*sina+sing*cosa
  	rotation(0,2)=-cosg*sinb
  	rotation(1,0)=-sing*cosb*cosa-cosg*sina
  	rotation(1,1)=-sing*cosb*sina+cosg*cosa
  	rotation(1,2)=sing*sinb
  	rotation(2,0)=sinb*cosa
  	rotation(2,1)=sinb*sina
  	rotation(2,2)=cosb
  	do i=0, 2
    		rp_hat(i)=0.0
    		kp_hat(i)=0.0
    		do j=0, 2
      			rp_hat(i)=rp_hat(i)+rotation(i,j)*r_hat(j)
      			kp_hat(i)=kp_hat(i)+rotation(i,j)*k_hat(j)
    		enddo
  	enddo
  	theta=atan2(sqrt(rp_hat(0)*rp_hat(0)+rp_hat(1)*rp_hat(1)),rp_hat(2))
  	phi=atan2(rp_hat(1),rp_hat(0))
  	call reduce(theta,phi)
  	sint=sin(theta)
	cost=cos(theta)
  	sinp=sin(phi)
	cosp=cos(phi)
  	sinz=-(kp_hat(0)*cost*cosp+kp_hat(1)*cost*sinp-kp_hat(2)*sint)
  	cosz=-kp_hat(0)*sinp+kp_hat(1)*cosp
	zeta=atan2(sinz,cosz)
	return
	end
