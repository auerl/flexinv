      subroutine anisoker_ttpath_topo(dep,del,ray,trueray,idersv,xtu,qvec,ierr)
c
c  DESCRIPTION:
c  This subroutine computes travel times,etc and derivatives of travel 
c  time w.r.t. the model parameters the travel time calculation is 
c  based on a major rewrite of tpg50 by guy masters, oct 97.  This is 
c  modified to allow for an anisotropic inversion using B-splines, 
c  J.G., 2002.
c
c  CHANGELOG: 
c  This version does not support diffracted waves J.G.
c
c  INPUT:
c  dep    - source depth in km
c  ray name (will ask of SH or SV for pure S rays)
c  iso    - =0(isotropic layers remain isotropic)/ =1(otherwise)
c  imod   - =0(anisotropic model is used)/ =1(averaged isotropic model)
c  idersv - =1 if derivatives are to be calculated
c  ierr  - =1 if there is an error
c
      implicit double precision (a-h,o-z)
c
c  To build kernel, only need lower precision
c  ray1 contains the true name of the phase and ray is the modified in
c  case it is a SS (PP) precursors or PdS (SdP) converted phases
c
      character*256 ray, trueray
      character*128 outfl
      dimension qvec(500)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)
      common/ray$/pmin,pmax,nps(20,2),npsrc(2),nn,iup,jref,npleg,nsleg
      common/rayinfo/theta_bt(20),rad_bt(20),nbot,ithbot,thetamax,
     &     rayseg_th(20), nrayseg_type(20),nseg,delreq
      common/xnorm/xnorm        
c
c  /savepath/ saves combined raypath info for kernel calculation
c  note for nrayty array, element = 1 for P, element = 2 for S
c
      common/anisopath/theta_path(50000), rad_path(50000),ntheta ! L.A. 3200->50000
c
c  In /kernel/, ifsplit= radial function type (1=continuous 2=split)
c  ifdiff=1 (diffracted wave)
c
      common/kernel/sum(362,20,10),ifsplit,ifdiff
      common/premvel/vphp,vpvp,vshp,vsvp,etap,rhop,vpiso,vsiso
c
c  commons added to deal with SS precursors or converted phases
c
      common/secphases/isectype,rsec,indexsec,lsec
      common/slow1/l220,l400,l520

      integer iradfun
      external f




c...............................................................................
c *** Set the following without asking
c *** Order of GL integration
      idleg=6
      iso=1       ! means we compute isotropic derivs for isotropic layers 
      imod=0      ! imod =0 for anisotropic layer, =1 for equiv. isotropic model

c *** Assume continuous radial functions as default
      ierr = 0
      ifdiff=0	! no diffracted wave allowed
      lsrc=0

      rdep=6371.d0-dep
      delreq = del
      do i=1,numlyr
        if(xb(i).lt.rdep.and.xt(i).ge.rdep) lsrc=i
      enddo

c 
c *** case iopt=1, computes a table with evenly spaced deltas
  120 p1=pmin
      p2=pmax
c *** in d1 and d2 we put the values of d(delta)/d(p) for p=p1 and
c *** p=p2; if they have opposite signs, d(delta)/d(p) must be zero in
c *** between. the value for which this happens is found by a call
c *** to subroutine zero then the user is prompted for a new p range
      itim=2
c
c *** ider=0 here integrates only ttime and t* in deriv (see fqanis.f for details)
c *** routine below returns:
c *** qvec(1)=distance, qvec(2)=dxdp, qvec(3)=time, qvec(4)=t*
c *** xtu = turningray radius
c
      ider=0   ! here ider is set to 0 to get time and t* only, no kernels
      call deriv(p1,qvec,xtu)
      delex1=qvec(1)
      d1=qvec(2)
      call deriv(p2,qvec,xtu)
      delex2=qvec(1)
      d2=qvec(2)
c *** check if d1 and d2 have opposite sign
      if(d1*d2.lt.0.d0) then
         ptarg=0.d0
c *** If opposite, finds root and put in p3
         call zero(p3,p1,p2,f,ptarg)
         write(*,*) '    old pmin = ', pmin, '  old pmax = ', pmax 
         write(*,*) '    d(delta)/d(p) vanishes at p= ', p3
         write(*,"('enter new pmin,pmax (neg values to quit):')")
         read(*,*) pmin,pmax
         if(pmin.lt.0.d0.or.pmax.lt.0.d0) stop
         goto 120
      end if
c *** Monotonic
      dismin=delex1
      dismax=delex2
      if(delex1.gt.delex2) then
         p1=pmax
         p2=pmin
         dismin=delex2
         dismax=delex1
      end if
c *** Use 0.5 degree as a cursion
      dmin=dismin+0.5
      dmax=dismax-0.5  
      if(delreq.lt.dmin.or.delreq.gt.dmax) then
c	   print*, 'out of range, min=',dmin, ' max=',dmax, ' delta=', delreq
	   ierr = 1
	   return
      endif 

c
c *** Set parameters: ider=idersv to make kernel calculation available.
c *** itim:    1 --> computes delta given p
c ***          2 --> computes p and dx/dp 
c ***          3 --> computes both of the above + time and tstar
c

      itim=1
      ider=idersv
      if(ifdiff.eq.0) then
	    call zero(pp,p1,p2,f,delreq)
          if(pp.lt.0.d0) then
            print*, 'stopped ad second ierr'
		ierr = 1
		return
	    endif
      endif
      itim=3

c *** note: ider = idervs = 1, so he calculates derivs of  tts wrt model params
      call deriv(pp,qvec,xtu)
      call prray_kernel(trueray,pp,qvec(1),qvec(5),xtu,iopt)  ! only prints some things
      call raypath_kernel_topo(ray)  ! computes the raypath

      return
      end

