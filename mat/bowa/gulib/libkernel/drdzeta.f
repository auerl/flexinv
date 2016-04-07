      subroutine drdzeta(r,drdz,ifrelative,th,viso)
c  returns the dr/dz for radius r.
c  xa, xc, xn, xl, eta are A, C, N, L, eta respectively in PREM,
c  Input:  	r     ---  radius
c		ifrelative --- 1=relative perturbation  2= absolute ..
c		th   ---  current running angle along ray path
c  Output:
c		drdz   --- dr/dzeta for apparent velocity
c			   1. ph  2. pv  3. sh  4. sv  5. eta
c		viso  --- either S or P velocity, depending on ray type
c 			  as well as anisotropy information
c
      implicit double precision (a-h,o-z)
      logical isotrp
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/coeff/coef(4,8,20)
      common/premvel/vphp,vpvp,vshp,vsvp,etap,rhop,vpiso,vsiso
      common/rayinfo/theta_bt(20),rad_bt(20),nbot,ithbot,thetamax,
     &		rayseg_th(20),nrayseg_type(20),nseg,delreq
      dimension drdz(1)
      parameter (pi2 = 3.141592653579*0.5d0)
c
c	print*, 'r=', r, '   rnorm=', rnorm
      do j=1, numlyr
	      if(r.gt.xb(j).and.r.le.xt(j)) then
		      iq = j
		      y=r/rnorm
		      ierror = 0
		      goto 999
	      endif
      enddo
      do i=1, 5
	  drdz(i)=0.d0
      enddo
      vpeq=0.d0
      vseq=0.d0
      vsiso=0.d0
      vpiso=0.d0
      viso=0.d0

999   if(ierror.eq.1) then
	      write(*,*) 'error in the raypath in getgreen.f'
	      stop
      endif
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
c
c  the following saves the anisotropic (or isotropic) model values
c
      vph0 = vph
      vpv0 = vpv
      vsh0 = vsh
      vsv0 = vsv
      eta0 = eta
c
c put them in common block
c
      vphp = vph
      vpvp = vpv
      vshp = vsh
      vsvp = vsv
      etap = eta
      rhop = rho
c
      rsq = r*r
      psq = p*p
      rp_ratio = r/p
c
      xc=rho*vpv*vpv	! C (check PREM paper)
      xl=rho*vsv*vsv	! L
      xa=xc		! A = C if isotropic
c
c  below is the isotropic velocity constructed from anisotropic
c  components
c
      if(.not.isotrp) then
        xa=rho*vph*vph	! A
        xn=rho*vsh*vsh	! C
        xf=eta*(xa-2.d0*xl)	! F = eta*(A-2*L)
c        if(imod.ne.0) then
           xkapa=(4.d0*xa+xc+4.d0*xf-4.d0*xn)/9.d0
           xmu=(xa+xc-2.d0*xf+5.d0*xn+6.d0*xl)/15.d0
           xa=xkapa+4.d0*xmu/3.d0
           xf=xkapa-2.d0*xmu/3.d0
           xn=xmu 
           vph=dsqrt(xa/rho)
           vpv=vph
           vsh=dsqrt(xn/rho)
           vsv=vsh
	   vpeq=vph   ! these are equivalent velocities
	   vseq=vsh
           eta=1.d0
c	endif
      endif
c             below finds the apparent velocity vs for
c             anisotropic case, use Woodhouse (1984), tan(i) = r*qdel
c
c        xn=rho*vsh*vsh		! N
c	temp = qtau*xl
c	qdel = p*xn		! numerator of qdel
c	if(dabs(temp).gt.0.00001) then
c		qdelr=qdel/(temp*r)  ! p*N/(r*qtau*L)
c		ang=pi2*0.5-datan(qdelr)
c		ang=datan(qdelr)
c		print*, 'ang =', ang, ' r=', r, '  qdelr =', qdelr
c		if(ang.lt.0.d0) ang=-ang   !postive and negative equal
c		sang = dsin(ang)
c		cang = dcos(ang)
c		vseq = dsqrt((sang*sang*vsh*vsh+cang*cang*vsv*vsv))
c	else
c		vseq = vsh    ! if horizontally traveling, use vsh
c	endif
c
c  this is another possible formulation for equivalent shear
c  velocity
c	vseq = sqrt((vsh0*vsh0 + vsv0*vsv0)*0.5d0)
c	vpeq = sqrt((vph0*vph0 + vpv0*vpv0)*0.5d0)
cc
c
c
c -------finds the ray type----------------------
c
	th1=0.d0
	do i=1, nseg
		th2=rayseg_th(i)
		if(th.gt.th1.and.th.le.th2) then
			ix=nrayseg_type(i)
			goto 10
		endif
		th1=th2
	enddo
10	continue
	if(ix.eq.1) then
c****  SH-wave
           if(isotrp) then
	      vseq=vsv0
	      viso=vsv0
	   else
	      viso=vsh0
	   endif
	else if(dabs(vsv0).lt.0.00001d0) then
c*** P in fluid
c*** reflected SV at the fluid
	      if(ix.eq.3) then
		viso = vpv0
	      else
	      	vpeq=vph
	      	viso=vpv0
	      endif
        else
c*** P/SV in solid 
	   if(isotrp) then
		if(ix.eq.2) then
		   vpeq=vpv0
		   viso=vpv0
	        else
		   vseq=vsv0
		   viso=vsv0
		endif
	   else
		if(ix.eq.2) then
		   viso=vpv0
		else
		   viso=vsv0
		endif
	   endif
        endif

c
c  the following calculates drdz = r*sqrt(r*r-V*V*p*p)/(V*p)
c  The velocity V uses equivalent shear and compressional velocities
c  which are computed 
c
      temp1 = vpeq*vpeq*psq
      
      if(rsq.ge.temp1)then
        drdz(1) = rp_ratio*dsqrt(rsq-temp1)
      else
        drdz(1) = 0.d0
      endif
      if(dabs(vpeq).ge.0.0001) then
	 drdz(1)=drdz(1)/vpeq
      else
	 drdz(1)=0.d0
      endif
      drdz(2)=drdz(1)


      temp2=vseq*vseq*psq
      if(rsq.ge.temp2) then
         drdz(3)=rp_ratio*dsqrt(rsq-temp2)
      else
         drdz(3)=0.d0
      endif
      if(dabs(vseq).ge.0.0001) then
	  drdz(3)=drdz(3)/vseq
      else
	  drdz(3)=0.d0
      endif

      drdz(4)=drdz(3)
      drdz(5)=1.d0

c
c  the following multiplies by a model vector to allow the use
c  of  relative perturbation rather than absolute model perturbation
c

      if(ifrelative.ne.0) then
         drdz(1)=drdz(1)*vph0
         drdz(2)=drdz(2)*vpv0
         drdz(3)=drdz(3)*vsh0
         drdz(4)=drdz(4)*vsv0
         drdz(5)=drdz(5)*eta0
      endif
      
      return
      end
