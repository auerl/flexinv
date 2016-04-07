      subroutine addpath_kernel(th,k,i1,i2,kup)
c  This version calculates the raypath and save the info in
c  theta_path, rad_path and ntheta in order to calculate kernels kIn.
c  input:
c   k       = 1 for P leg, 2 for S leg
c   i1,i2 = indices in depth to go between
c   kup     = 1 for downgoing wave, -1 for upgoing (strange!)
c   th is the current distance of raypath
c   OUTPUT:  theta_path = theta increments of the path
c   	     rad_path   = radius corresponding to theta increment
c	     ntheta     = number of segment in the path
c	     raytp      = the type of wave (P or S) for this leg
c	     rayst 	= the starting number of this leg
c	     ileg	= the counter for the number of ray leg
c   Note:    ntheta increases due to contributions from each leg
c
      implicit real*8(a-h,o-z)
      common/path$/delt(2,800,2),nfin(2),kdep
      common/anisopath/theta_path(1600),rad_path(1600),derv_path(1600,10)
     +   ,ntheta, nraytp(20), nrayst(20),ileg
      common/kernel/sum(362,20,10),ifsplit, kerpath(10,800,2)
      real*8 kerpath(10,800,2)
      i2s = min(i2,nfin(k))
      if (kup.eq.1) then
c downgoing leg
        do l=i1,i2s
          rray=delt(1,l,k)/6371.
          th=th+delt(2,l,k)
	  if(rray.ne.0.d0) then
	  	ntheta = ntheta + 1
	 	rad_path(ntheta) = delt(1,l,k)
	  	theta_path(ntheta) = th
		do kk=1, 5
		    derv_path(ntheta, kk) = kerpath(kk,l,k)
		enddo
	  endif
	  if(l.eq.i2s) then
		ileg = ileg + 1
		nrayst(ileg) = ntheta
		nraytp(ileg) = k
	  endif
        enddo
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
		do kk=1, 5
		    derv_path(ntheta, kk) = kerpath(kk,l-1,k)
		enddo
	  endif
	  if(l.eq.i1s) then
		ileg = ileg + 1
		nrayst(ileg) = ntheta
		nraytp(ileg) = k
	  endif
        enddo
      endif
      return
      end
