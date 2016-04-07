      program build_ttime_splines_topo
c
c DESCRIPTION:
c Calculate the travel time A matrix for a model with B-spline parameterization. 
c This version allows the joint inversion of perturbations to Vsh, Vsv
c Vph, Vpv relative to a given reference Earth model.  The code is written
c for non-varing resolution only, Yu J. Gu. 2002.
c
c CHANGELOG:
c does not compute the ATA matrix Y.G.
c allows to output ray sensitivity/paths L.A.
c allows for variable layer thickness L.A.
c allows to use a 2x2 degree parameterization L.A.
c added some comments and made everything more readable L.A.
c removed some obsolete stuff (tesselation, radial splines) L.A.
c layer information is now read from a file L.A.
c is supposed to be able to take major arc data L.A.
c

      implicit double precision (a-h, o-z)

c
c basic parameters, naming is obsolete
c

      parameter(maxfreq=8)
      parameter(maxnode=maxfreq*maxfreq*10+2)
      parameter(maxrknot=50)
      parameter(maxparm=maxnode*maxrknot)
      data pi /3.141592653579/
      data rad/57.2957795130824/
      data reprad /57.295779579/    


c
c commons for the path and kernels      
c
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)
      common/ray$/pmin,pmax,nps(20,2),npsrc(2),nn,iup,jref,npleg,nsleg
        common/rayinfo/theta_bt(20),rad_bt(20),nbot,ithbot,thetamax,
     &     rayseg_th(20), nrayseg_type(20),nseg,delreq
      common/xnorm/xnorm        
      dimension rowps(4,50000),indro(50000)
c rowsh(50000),rowsv(50000)

c
c /savepath/ saves combined raypath info for kernel calculation
c note for nrayty array, element = 1 for P, element = 2 for S
c
      common/anisopath/theta_path(50000), rad_path(50000),ntheta
      common/parm/rho,vpv,vsv,vph,vsh,eta,s(5),r
      common/parm1/xa,xc,xn,xl,xf,tdif(5)
      common/isot/isotrp
      common/partvel/pathder(5,50000),pathvel(5,50000)

c
c in /kernel/, ifsplit= radial function type (1=continuous 2=split)
c ifdiff=1 (diffracted wave)
c
      common/kernel/sum(362,20,10),ifsplit,ifdiff
      common/premvel/vphp,vpvp,vshp,vsvp,etap,rhop,vpiso,vsiso

c
c commons added for spline integrations only, /spathinfo/ contains
c velocity along path, 
c
      common/spathinfo/velpath(4,50000),kradpath(50000),devpath(5,50000)
      common/amatrix/arow(maxparm),indarow(maxparm),nonzero_a

c
c this dominates the inversion options, e.g., which velocity to invert
c
      common/invopt/numinvvar,invparlist(5)

c
c commons added to deal with SS precursors or converted phases
c
      common/secphases/isectype,rsec,indexsec,lsec

c
c the dimension of 30000 should be perserved to be consistent 
c with calc_tessa(). ray1 is modified ray name depending on 
c whether it is a secondary phase like SdS or Pds, etc.
c
      parameter(maxspnode=30000)
      parameter(maximpact=20)
      dimension bslat(maxspnode),bslon(maxspnode)
      character*256 resfile,amatname,qkfile,stfile,dtfile
      character*256 ray, ray1, raybuf(2)
      character eqname*8,stname*4
      dimension rknot(maxrknot)
      dimension thpath(50000),radpath(50000)
      dimension qvec(500)
      dimension isplayer(maxrknot,maxrknot)
      dimension arytrans(maximpact), aryrefl(maximpact)
      integer   numinvvar,iadapt,ishbuf(2)
      real*4    elat,elon,edep,slat,slon,stelv,res,datum,
     &		ttime,elipc,crustcor,evcor,delta,coric,pred

      real*4    row0 !rowsh0,rowsv0,
      integer*4 indro0

c
c commons for voxel parameterization
c
      parameter(rearth=6371.d0,rcmb=3480.d0) ! careful, this is hardwired
      parameter(nlatzomax=180,nlaym=50) ! used to be 29
      parameter(maxblocks=31808,numpar=2) ! corresponding to 1.25 degree
      parameter(maxpar=maxblocks*nlaym*numpar)
      parameter(eqthick=0) ! 0: read from file, 1: equal inc, 2: 20l, 3: 26 l 
      parameter(iswit=1)
      dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
      dimension rbnd(0:nlaym)
      real*4 rescale(nlaym)
      real*4 laypts(30)
      real*4 angle_hitcount(4,maxpar) 
      character*80 rad_fil
      character*120 outdir

c
c some testing parameters
c
      dimension qplud(5)
      parameter(writeact=1)
      parameter(delim=123456)
      parameter(plotraysens=1)

c
c a couple of initializations
c
      ifsplit = 0
      ifdiff  = 0    ! assume non-diffracted waves
      ifwate  = 0    ! no weighting based on igrade necessary for major phases

c
c read info on voxel parameterization
c
      print*,"Enter size of pixels:"
      read*,eq_incr 
      print*,eq_incr
      print*,"Enter if adaptive (1) or regular (0) grid"
      read*,iadapt

      if(iadapt.eq.0)then
         compgrid=eq_incr
      elseif(iadapt.eq.1)then
         compgrid=5.0 ! reference grid with which fine grid shall be compatible
      endif

      nlatzones=180./eq_incr
      if(nlatzones.gt.nlatzomax)stop "Pixels too small"
      print*,"Enter folder to store matrices:"
      read(*,*) outdir
      print*, outdir
      print*,"Enter layer file path:"
      read(*,*),rad_fil
      print*,rad_fil


c
c define nsqrs,nsqtot and layers; compgrid is the refgrid
c

        call param(eq_incr,nsqrs,nsqtot,nlatzones,n1layer,iswit,compgrid,nlatzomax,iadapt)

c
c define radial parametrization
c
      rbnd(0)=rearth
      if(eqthick.eq.0)then
         ! Read radial parametrization from layers.in (default)
         print*,rad_fil(1:lnblnk(rad_fil))
         open(333,file=rad_fil(1:lnblnk(rad_fil)),status='old')
         read(333,*),nlay
         if(nlay.gt.nlaym)then
            print *,'Nr. of layers is out of bounds ',nlay,nlaym
            stop
         endif
         do il=1,nlay
            read(333,*),layer
            laypts(il)=layer
         enddo
         close(333)
         print*,"Layers (bottoms): "
         do iyn=1,nlay
             rbnd(iyn)=rearth-dble(laypts(iyn))
             print*,rbnd(iyn)
         enddo                 
      elseif(eqthick.eq.1)then
         ! Equally thick layers from crust to cmb
         yncr=(rearth-rcmb)/dble(nlay)      
         print*,"Layers (bottoms): "
         do iyn=1,nlay
            rbnd(iyn)=rearth-(yncr*dble(iyn))
            print*,rbnd(iyn)
         enddo
      endif

c
c some verbose
c
       nvx=n1layer*nlay
       print*,nvx," voxels in total!"

c
c flush dws
c
       do i=1,nvx*2
          do j=1,4
             angle_hitcount(j,i)=0.
          enddo
       enddo


c
c open the binary files in which matrix will be stored
c


      open(111,file=trim(outdir)//'/tmp/a.vx_mat',access='direct',recl=4,form='unformatted')
      open(112,file=trim(outdir)//'/tmp/a.vx_ind',access='direct',recl=4,form='unformatted')
      open(113,file=trim(outdir)//'/tmp/a.vx_pnt')
      open(114,file=trim(outdir)//'tmp/a.vx_vec')
      if(writeact.eq.1)then
         open(808,file=trim(outdir)//'/tmp/quakes_actu.dat')
         open(909,file=trim(outdir)//'/tmp/receiv_actu.dat')
      endif

      write(*,"('Number of velocity vectors to invert [max=4]:')")
      read(*,*) numinvvar
      if(numinvvar.gt.5.or.numinvvar.lt.1) stop 'Unreasonable number of velocities!'
      write(*,"('***List: 1=Vph,  2=Vpv,  3=Vsh,  4=Vsv')")
      do i=1, numinvvar
      	 write(*,"('enter selection ', i3, ' :')") i
	 read(*,*) ii
	 if(ii.gt.4.or.ii.lt.1) stop 'Option not found!!!'
	 do j=1, i-1
            if(invparlist(j).eq.ii) stop 'This has been selected already!!'
	 enddo
	 invparlist(i)=ii
      enddo
	
      do i=1, numinvvar
         print*, 'parameter to invert:  i=', i, '  option =', invparlist(i)
      enddo 

c
c read reference model
c
      call readmod


c
c we must set up source layer before calling raynam
c To trick the raynam routine, use a hypothetical dep=1.0
c This is just a hook and does have the danger to bypass some of 
c the functionality of the depth checking of the raynam routine.
c

      lsrc=0
      dep=1.0 	      
      rdep=6371.d0-dep
      do i=1,numlyr
        if(xb(i).lt.rdep.and.xt(i).ge.rdep) lsrc=i
      enddo
      write(*,"('Compute A matrix with 400 and 670 topography [1=yes, 0=no]?')")
      read(*,*) iftopo
      write(*,"('Is this an absolute time data set [1=yes, 0=no]?')")
      read(*,*) ifabsolute
      numray=1
      if(ifabsolute.ne.1) numray=2

c
c get ray types and set up the global variables "ray" and "ish" in
c raybuf() and ishbuf() respectively to allow looping over differential
c times.
c

      do iray=1,numray
      	 ray=' '
	 if(numray.eq.1) then
      	 	write(*,"('enter the ray name:')")
	 else
		if(iray.eq.1) then
      	 		write(*,"('enter the slower traveling ray name:')")
      	 		write(*,"('e.g., for SS-S time, enter SS here..')")
		else
      	 		write(*,"('enter the faster traveling ray name:')")
      	 		write(*,"('e.g., for SS-S time, enter S here..')")
		endif	
	 endif

     	 read*,ray ! lapo 17.2.2010
	 raybuf(iray)=ray
	 ifsecond=ifsecondary(ray,ray1)

	 if(ifsecond.ne.0) ifwate=1
         call raynam_topo(ray1,ierr0)
      	 if(ierr0.ne.0) then
           write(*, "('error in ray name')")
           stop
      	 end if
      	 ish=0
         if(ifanis.ne.0) then

c
c Check to see if this is pure S, then ask if SH or SV
c
               if(npleg.eq.0) then ! this is S SS SSS SSSS
                  write(*,"('We deal with a ray not having a p leg')")
                  write(*,"('choose wave type (0 = SV, 1 = SH)')")
                  write(*,"('note: generally should be SH')")
                  read(*,*) ish
                  if(ish.ne.1) ish=1
               else ! this is SKS SKKS 
                  write(*,"('We deal with a ray having a p leg')")                  
                  write(*,"('wave type is set to (0 = SV)')")
                  read(*,*) ish ! this is just a "trick" from yu, L.A. 
                  ish=0
               end if
          end if
	  ishbuf(iray)=ish
      enddo

c
c Read desired distance rage, the residual file and major or minor arc flag
c
      write(*,"('desired distance range[if not sure, enter in 0--180]:')")
      read(*,*) dist0, dist1     
      write(*,"('distance range :', f8.2,f8.2)") dist0, dist1
      write(*,"('enter residual file:')")
      read(*,"(a)") resfile 
      write(*,"('enter 0 or 1 for minor or major arc:')")
      read(*,*) arc 

c
c Some increment variables
c
      nrec=1
      irec=0
      isucc=0
      ifail=0
      nerrvoxint=0
      nrsenspath=1   ! only plot 100 raypaths per subset
c
c Read data in Lapos format
c
      dtfile=resfile(1:lnblnk(resfile))//".dtprem.bin"
      qkfile=resfile(1:lnblnk(resfile))//".quakes.bin"
      stfile=resfile(1:lnblnk(resfile))//".receiv.bin"
      print*,"opening",dtfile,qkfile,stfile

      open(921,file=dtfile,access='direct',recl=4*4,status='old',form='unformatted')
      open(422,file=qkfile,access='direct',recl=4*12,status='old',form='unformatted')
      open(522,file=stfile,access='direct',recl=4*8,status='old',form='unformatted')


c
c ^^^^^^^^^^^^^^^ MAJOR "GOTO LOOP" ^^^^^^^^^^^^^^^
c 

10    read(921,rec=nrec,err=100)datum
      read(422,rec=nrec)elon,elat,edep
      read(522,rec=nrec)slon,slat

c
c to test certain ranges of the dataset
c      if (nrec.lt.101168) then
c         goto 90
c      end if
c      print*,nrec

c
c escape sequences, stop criteria, 2nd part was added by ludwig
c
      if(slat+90..lt.0.1)then
         slat=-89.
      endif
      if(edep.gt.660.)then
         print*,"Earthquake located too deep, discarding datum: ", nrec
         goto 90
      elseif(edep.lt.0)then
         print*,"Earthquake depth below 0, discarding datum: ", nrec
         goto 90
      endif      
      if(abs(res).gt.20.0) then
	  write(*,*) 'nrec=', nrec, 'measurement unreasonable!! Skip!'
	  goto 90
      endif

      delat=elat
      delon=elon
      dslat=slat
      dslon=slon
      dedep=edep

c
c compute the epicentral distance
c
      call conv2geocen(delat,delon,sth,sphi)
      call conv2geocen(dslat,dslon,rth,rphi)
      del=dist(sth,sphi,rth,rphi)*reprad

c
c convert epdist in case of major arc 
c
      if (arc.eq.1) then 
         del=360-del
      endif

c
c escape sequence, when del outside def range 
c
      if (del.lt.dist0.or.del.ge.dist1) then
         print*,dist0,dist1,del
         print*, "del in radians outside defined range!"
         goto 90 
      endif

c
c this block finds the source layer
c
      lsrc=0
      rdep=6371.d0-dedep
      do i=1,numlyr
        if(xb(i).lt.rdep.and.xt(i).ge.rdep) lsrc=i
      enddo
c
c initialize a row of the A matrix
c
      do i=1, numatd
	 arow(i)=0.d0
      enddo

c
c this block loops over the ray types
c
      do iray=1, numray

         ray=raybuf(iray)
	 ish=ishbuf(iray)
	 ifsecond=ifsecondary(ray, ray1)

c
c identify raytype
c
             call raynam_topo(ray1,ierr0)
             if(ierr0.ne.0) then
                write(*,"('problem with raynam routine for ', i6)") nrec
                goto 90
             endif

c
c compute the travel time, ray path and derivatives of travel
c time with velocity.  The "1" before ierr here chooses to use anisotropic 
c velocity. qvec(1)=delta, qvec(2)=dX/dp, qvec(3)=time, qvec(4)=tstar
c        
      	     call anisoker_ttpath_topo(dedep,del,ray1,ray,1,xtu,qvec,ierr)
             if(ierr.ne.0) then
c                print*,"error in anisoker"
                goto 90
             endif
        
c
c output ray sensitivity and raypath
c
             
             if(plotraysens.eq.1) then
                  if (nrsenspath.le.100) then
                     open(666,file=trim(outdir)//"/tmp/sens.txt")
                     do tiefe=1.,6371.,1.
                        call calcqvec(tiefe,qplud)
                        write(666,"(f20.7,1x,5(f14.7,1x))"),qplud(1),qplud(2),qplud(3),qplud(4),qplud(5)   
                     enddo
                     write(666,"(f20.7,1x,5(f14.7,1x))"),delim,delim,delim,delim,delim
                     nth=ntheta
                     open(777,file=trim(outdir)//'/tmp/ray.txt')
                     do j=1,nth
                        thpath(j)=theta_path(j)
                        radpath(j)=rad_path(j)
                        write(777,"(f20.7,1x,5(f14.7,1x))"),thpath(j)*180/3.1415926,radpath(j)
                     enddo
                     write(777,"(f20.7,1x,5(f14.7,1x))"),delim,delim
                     nrsenspath=nrsenspath+1
                  endif
             endif

c
c project kernels on voxels and computes kernels (not done in anisoker)
c but in anisoker, a bit unlogical. 
c



      	 call voxel_int(delat,delon,dedep,dslat,dslon,bslat,bslon,ngpt,hh,
     &		 rknot,numknot,numsplit,isplayer,xtu,iray,nsqrs,nsqtot,rbnd,
     &		 nlatzones,nlay,nlatzomax,nlaym,eq_incr,n1layer,rowps,
     &           indro,nnz,ierrvoxint,nvx,angle_hitcount)
                   if(ierrvoxint.ne.0)then
                      nerrvoxint=nerrvoxint+1     
                      print*,"datum:",nrec,nerrvoxint,"discarded data for errors in voxel_int"
                      goto 90
                   endif


c
c write out actually used data in case requested
c
         if(writeact.eq.1)then
            write(808,*),elon,elat,edep
            write(909,*),slon,slat
         endif
 
c 
c update matrix, this part was hardcoded to vsh/vsv, now allows to select an
c arbitrary number of inversion parameters b/w 1 and 4 to invert for, LA 05/14
c  
          if(nnz.gt.50000) then
             print*,"ERROR nnz is larger than matrix dim"
             stop
          endif


          ! Loop over all voxels affected by current ray
	  do inz=1,nnz

            ! convert to short 4 byte int
            indro0=indro(inz)

            ! loop over selected inv params
            do iniv=1,numinvvar
               ipar=invparlist(iniv)              
               row0=sngl(rowps(ipar,inz))
               irec=irec+1
               write(111,rec=irec)row0
               write(112,rec=irec)indro0+nvx*(iniv-1)
            enddo
                           
            if(indro(inz).gt.nvx)then
               print*,"index out of bounds, core parameterized?"
               print*,indro(inz)
            endif

c
c display out a row of the matrix, for diagnostics
c	write(*,"(2(1x,i6),2(1x,f15.7))")irec,indro(irec),sngl(rowsh(irec)),sngl(rowsv(irec))
c	write(*,*)irec,indro(irec),sngl(rowsh(irec)),sngl(rowsv(irec))
          enddo


c
c write out measurements and irec 
c
      write(113,*)irec
      write(114,*)datum


c
c some things to deal with topgraphy
c Ludwig has never used these parts of the
c code and doesn't know whether they are working
c

	if(iftopo.ne.0) then
c surface topography
	 	pp=p
c 400 discontinuity
	 	rdisp=5971.d0
	 	call find_impactpt(rdisp,arytrans,maximpact,
     &			ntrans,aryrefl,maximpact,nrefl)
	 	call addsurfspline(delat,delon,dslat,dslon,
     &		 bslat,bslon,ngpt,hh,arytrans,
     &		 ntrans,aryrefl,nrefl,pp,ray1,iray,rdisp,num3d)
c 670 discontinuity
	 	rdisp=5701.d0
	 	call find_impactpt(rdisp,arytrans,maximpact,
     &			ntrans,aryrefl,maximpact,nrefl)
	 	call addsurfspline(delat,delon,dslat,dslon,
     &		 bslat,bslon,ngpt,hh,arytrans,
     &		 ntrans,aryrefl,nrefl,pp,ray1,iray,rdisp,num3d+ngpt)
	endif

      enddo ! End of loop over raytype

      isucc = isucc + 1	

c
c uncomment lines in savearow_topo to get printout of matrix row
c     call savearow_topo(ioamat,delat,delon,dslat,dslon,dedep,finalres,qvec(3),
c     &	numatd,ig,ifwate)
c

      goto 12
90    ifail=ifail+1                                                       
12    nrec=nrec+1

      verbose=verbose+1
      if (verbose.ge.1000) then
         print*,"Processing record nr: ",nrec
         verbose=0
      endif

      goto 10
100   continue

      open(345,file=trim(outdir)//'/tmp/a.vx_dws')
      do i=1,nvx*2

c
c Note: Eventually I decided to use the hitcount instead of the dws (derivative weight sum) to define
c my paramterization. Surface wave hitcounts are multiplied with a depth dependent function to accoun
c for their smaller sensitivity at depth
c

         write(345,'(4F17.4)') angle_hitcount(1,i),angle_hitcount(2,i),angle_hitcount(3,i),angle_hitcount(4,i)
      end do
      close(345)


c.......................................................................

      close(ioamat)
      write(*,*) 'Rows computed        : ',isucc
      write(*,*) 'Failures             : ',ifail
      call system("date")
      stop
      end
