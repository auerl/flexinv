      subroutine raypath_kernel_topo(ray)
c Reads the ray name, breaks the ray down into its constituent parts,
c and calls subroutine addpath_kernel to compute the ray path. No checking
c of validity of rayname is performed as this has been done in raynam.
c The final raypath is saved in /anisopath/ in order to calculate 
c kernels kIn in Dziewonski, 1984, or spline kernels.
c This version deals with secondary phases as well, YG 2002.
c
      implicit real*8(a-h,o-z)
      common/coeff/coef(4,8,20)
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/path$/delt(2,50000,2),nfin(2),kdep !la 1600->50000
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/anisopath/theta_path(50000),rad_path(50000),ntheta  !la 3200->50000
      common/kernel/sum(362,20,10),ifsplit,ifdiff
      common/rayinfo/theta_bt(20),rad_bt(20),nbot,ithbot,thetamax,
     &     rayseg_th(20),nrayseg_type(20),nseg,delreq
c----------------------------------------------------------------------
c...commons added to deal with SS precursors or converted phases
c
      common/secphases/isectype,rsec,indexsec,lsec
c----------------------------------------------------------------------
      integer flnum
      character*(*) ray      
      data if1/1/
c initialize for /anisopath/ parameters
	
      do i=1, 3200
	 theta_path(i) = 0.d0
	 rad_path(i) = 0.d0
      enddo
c  ----track ray path and bottoming points, YG, 98.------------------
      ntheta = 0
      ithbot=0   ! temporary layer for bottoming point, will
c		   be changed in addtoplayer.f
      nbot = 0   ! number of bottoming points
      nseg = 0
c--------------------------------------------------------------------
c*** track ray legs
      th=0.d0
      j1=0
      do 10 i=1,70
c ' PSIJKpsic' is the order in ichdec
c note j is 2 for P, 3 for S, 4 for I.  1 is reserved for a space
c k=2 for S, J or s types, otherwise k=1 (P waves)
c
      j=ichdec(ray(i:i),k)
      goto(10,11,11,12,12,13,14,14,100,100,10),j
c
c**** P or S ***
c
   11 goto(10,111,112,10,10,114,113,113,10,100),j1
      call addpath_kernel(th,k,kdep,klay(loc-1),1)
      goto 100
c---------------------------------------------------------------
c**** P (Up and Down going wave respectively) ***
c
  111 if(isectype.eq.0) then
c... major phase
	call addpath_kernel(th,1,1,klay(loc-1),-1)
      	call addpath_kernel(th,k,1,klay(loc-1),1)
      else
	if(isectype.eq.1) then
c... PP precursor
	   call addpath_kernel(th,1,klay(lsec-1)+1,klay(loc-1),-1)
      	   call addpath_kernel(th,k,klay(lsec-1)+1,klay(loc-1),1)		
	else if(isectype.eq.2) then
c... PS conversion
      	   call addpath_kernel(th,k,1,klay(lsec-1),-1)		
	endif
      endif
      goto 100
c---------------------------------------------------------------
c***  S (Up and Down going wave) ***
c
  112 if(isectype.eq.0) then
	call addpath_kernel(th,2,1,klay(loc-1),-1)
      else
c... SS precursors or SP conversions
	call addpath_kernel(th,2,klay(lsec-1)+1,klay(loc-1),-1)
      endif
  113 if(isectype.eq.0) then
	call addpath_kernel(th,k,1,klay(loc-1),1)
      else
	if(isectype.eq.1) then
c... SS precursors
	   call addpath_kernel(th,k,klay(lsec-1)+1,klay(loc-1),1)
	else if(isectype.eq.2) then
c... SP conversions
	   call addpath_kernel(th,k,1,klay(lsec-1),-1)
	endif
      endif
      goto 100
c---------------------------------------------------------------
  114 call addpath_kernel(th,1,klay(loc-1)+1,klay(lic-1),-1)
      goto 100
c
C**** I or J ****
c
   12 call addpath_kernel(th,k,klay(lic-1)+1,klay(numlyr),1)
      call addpath_kernel(th,k,klay(lic-1)+1,klay(numlyr),-1)
      goto 100
c
C**** K ****
c
   13 goto(20,132,132,100,100,131,20,20,100,20),j1
  131 call addpath_kernel(th,1,klay(loc-1)+1,klay(lic-1),-1)
  132 call addpath_kernel(th,1,klay(loc-1)+1,klay(lic-1),1)
      goto 100
C**** p or s ****
   14 call addpath_kernel(th,k,1,kdep,-1)
  100 j1=j
      k1=k
   10 continue
C****  receiver leg, not needed for PdS or SdP converted waves ***
      if(isectype.ne.2) call addpath_kernel(th,k1,1,klay(loc-1),-1)
c---------------------------------------------------------------
   20 thetamax=th
      return                  
      end
