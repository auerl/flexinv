      subroutine addpath_kernel(th,k,i1,i2,kup)
c  This version calculates the raypath and save the info in
c  theta_path, rad_path and ntheta in order to calculate kernels kIn.
c  input:
c   k       = 1 for P leg, 2 for S leg
c   i1,i2 = indices in depth to go between
c   kup     = 1 for downgoing wave, -1 for upgoing (strange!)
c   th is the current distance of raypath
c   OUTPUT:  theta_path = theta increments of the path
c   	     rad_path   = radius corresponding to theta increments
c	     nseg	= the counter for the number of ray leg
c 	     rayseg_th  = theta for each segment of ray
c	     nrayseg_type= ray type along path, useful in calculating
c                         qtau and tder.
c			  1 = SH;  2 = P;  3 = SV
c	     ntheta     = number of segment in the path
c   Note:    ntheta increases due to contributions from each leg
c
      implicit real*8(a-h,o-z)
      common/coeff/coef(4,8,20)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/path$/delt(2,50000,2),nfin(2),kdep !la 1600->50000
      common/anisopath/theta_path(50000),rad_path(50000),ntheta !la 3200->50000
      common/kernel/sum(362,20,10),ifsplit,ifdiff
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/rayinfo/theta_bt(20),rad_bt(20),nbot,ithbot,thetamax,
     &     rayseg_th(20), nrayseg_type(20),nseg,delreq
c      dimension rayset_th(20), nrayseg_type(20)
      dimension rayset_th(20)!lapo 21.1.2010
      i2s = min(i2,nfin(k))
      if (kup.eq.1) then
c downgoing leg
c ---------add first point near source, YG, 1998.-----------------
	if(ntheta.eq.0) then
		ntheta=ntheta+1
		rad_path(ntheta)=rdep
		theta_path(ntheta)=0.000001
	 endif
c ----------------------------------------------------------------
        do l=i1,i2s
          rray=delt(1,l,k)/6371.
          th=th+delt(2,l,k)
	  if(rray.ne.0.d0) then
	  	ntheta = ntheta + 1
	 	rad_path(ntheta) = delt(1,l,k)
	  	theta_path(ntheta) = th
	  endif
        enddo
c ---------------save ray type ------------------------------------
	nseg=nseg+1
	rayseg_th(nseg)=th
	if(k.eq.1) then
		nrayseg_type(nseg)=2
	else
		if(ish.eq.0) then
			nrayseg_type(nseg)=3	! SV wave
		else
			nrayseg_type(nseg)=1	! SH wave
		endif
	endif
	
	
c ------diffracted wave && bottoming points, YG, 98.---------------
	call check_bot(nrayseg_type(nseg),rad_path(ntheta),iflag)
	if(iflag.ne.0) then
		nbot = nbot +1
        	ithbot =  ntheta
		theta_bt(nbot) = th
c		write(*,*) 'nbot=', nbot, ' ithbot=', ithbot,  theta_bt(nbot)*180.d0/3.141592653579
	endif
c------------------------------------------------------------------
      else
c upgoing leg
        i1s=max(i1,2)
        do l=i2s,i1s,-1
          rray=delt(1,l-1,k)/6371.
          th=th+delt(2,l,k)
	  if(rray.ne.0.d0) then
	  	ntheta = ntheta + 1
	 	rad_path(ntheta) = delt(1,l-1,k)
	  	theta_path(ntheta) = th
	  endif
        enddo
c ---------------save ray type ------------------------------------
	nseg=nseg+1
	rayseg_th(nseg)=th
	if(k.eq.1) then
		nrayseg_type(nseg)=2
	else
		if(ish.eq.0) then
			nrayseg_type(nseg)=3	! SV wave
		else
			nrayseg_type(nseg)=1	! SH wave
		endif
	endif
c
c -------saves the missing point at the surface, YG, 1998.----------
c
	if(rad_path(ntheta).gt.xb(1)) then
		th = th+delt(2,i1s,k)
		ntheta = ntheta + 1
		rad_path(ntheta) = 6370.999999
		theta_path(ntheta) = th
cTEST
c	print*,"in addpath"
c	print*,"missing point:",ntheta,rad_path(ntheta),theta_path(ntheta)*57.295779579
c	print*,"*********************************************************"
	endif
c
c ---for diffracted waves, add the flat ray path,-------------------
c    use spacing of 1 degree,  YG, 1998.
c
	if(ifdiff.ne.0) then
		call diffract_path(th)
	endif
c ------------------------------------------------------------------
      endif
      return
      end
