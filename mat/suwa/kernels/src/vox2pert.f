c ***
c *** Reads voxel domain .pcn file and expands them in .pert
c *** Files on a 2x2 deg reference grid. pert files are used
c *** by prtlyrmod to perturb a Crust 2.0 + Refmod 1D model
c *** at a certain location. This allows to iterate!
c ***

	parameter(pi=3.1415926536)
	parameter(nmax_rad_spl=29)
	parameter(num_blocks_max=12436)
      parameter(nlatzomax=180.)
	parameter(iswit=1,refgrid=2.)
      parameter(nlatper=180./refgrid)
      parameter(nlonper=2*nlatper)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax)
	dimension vpcoeff(nmax_rad_spl*num_blocks_max)
	dimension vscoeff(nmax_rad_spl*num_blocks_max)
	character solutionp*80,solutions*80,miter*3
      character pertfile*80,type*1,answi*1
      character*80 outdir
	type='F' ! This type corresponds to "read layer info from file"

c *** Read parameters from stdin
99	print*,'PS model (1), SHSV model (2) or SH=PH/SV=PV model (3)?'
	read*,iso_choose
	if(iso_choose.ne.1.and.iso_choose.ne.2.and.iso_choose.ne.3)goto99
	print*,'Impose isotropic transition zone?'
	read*,answi
	if(answi.eq.'y'.or.answi.eq.'Y')then
	   print*,'How many isotropic splines?'
	   read*,iisos
	else
	   iisos=0
	endif
      print*,'Which eq_incr horizontal pixel size'
      read*,eq_incr
      print*,'Nominal number of blocks'
      read*,num_blocks
      print*,'Number of layers'
      read*,num_rad_spl
      print*,'Number of iteration'
      read*,iterationnr
      write(miter,'(i3.3)')iterationnr
      print*,'Directory to output .pert files'
      read(*,*),outdir
      do k=1,80
          if(outdir(k:k).eq." ") goto 706
      enddo
706   kchodir=k-1
      compgrid=eq_incr
      nlatzones=180./eq_incr
      ! Model .pcn files
	if(iso_choose.eq.1)then
          print*,'P model file'
          read(*,*),solutionp
          print*,'S model file'
          read(*,*),solutions
       elseif(iso_choose.eq.2.or.iso_choose.eq.3)then
          print*,'SH model file'
          read(*,*),solutionp
          print*,'SV model file'
          read(*,*),solutions
	endif
      do k=1,80
          if(solutionp(k:k).eq." ") goto 707
      enddo
707   kchpfil=k-1
      do k=1,80
          if(solutions(k:k).eq." ") goto 708
      enddo
708   kchsfil=k-1


c *** Define horizontal parametrization
      call param(eq_incr,nsqrs,nsqtot,nlatzones,numto,iswit,compgrid,nlatzomax)

      print*,solutionp(1:kchpfil),solutions(1:kchsfil)

c      pause
c *** Note that model coefficients are multiplied by 100.
	open(3,file=solutionp(1:kchpfil))
	do i=1,num_blocks*num_rad_spl
	   read(3,*),dummmy,vpcoeff(i)
c         print*,vpcoeff(i)
	enddo
	close(3)
	open(3,file=solutions(1:kchsfil))
	do i=num_blocks*num_rad_spl+1,2*num_blocks*num_rad_spl
	   read(3,*),dummmy,vscoeff(i-num_blocks*num_rad_spl)
c         print*,vscoeff(i-num_blocks*num_rad_spl)
c	   vscoeff(i-num_blocks*num_rad_spl)=
c     &	vscoeff(i-num_blocks*num_rad_spl)*100.
	   numbers=i-(num_blocks*num_rad_spl)
	   if(numbers.gt.((num_rad_spl-iisos)*num_blocks))then
	   vscoeff(i-num_blocks*num_rad_spl)=vpcoeff(i-num_blocks*num_rad_spl)
	   endif
	enddo
	close(3)

      print*,"came here"
	dummy=0.

c *** Loop over reference grid
	do ilat=1,nlatper
        do ilon=1,nlonper

            ! Convert to colat and mult by 100 (accuracy)
            INPLA=(90.-(ILAT*REFGRID))+(REFGRID/2.)
            INPLA=INPLA*100.
            INPLO=(ILON*REFGRID)-(REFGRID/2.)
            INPLO=INPLO*100.
            
            ! Find block id where lat and lon is
            iblo=superisqre(inpla,inplo,nsqrs,nsqtot,nlatzones,
     &	num_blocks,eq_incr)

            if(mod(iblo,100).eq.0)print*,'First',iblo,' blocks done'


            ! Write out pert files
            write(pertfile,"('/KER_',a3,'/L_',i3.3,i3.3,'_',a3,'.pert')")miter,ilat,ilon,miter
            open(4,file=outdir(1:kchodir)//pertfile)
            write(4,*)type,8,nmax_rad_spl
            do ispl=1,nmax_rad_spl
               index=iblo+(num_blocks*(ispl-1))
               if(ispl.gt.num_rad_spl)then
                  write(4,'(8e12.4)')
     &	      dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
               elseif(ispl.le.num_rad_spl)then
                  if(iso_choose.eq.2)then
	                write(4,'(8e12.4)')dummy,dummy,dummy,
     &	          vpcoeff(index),vscoeff(index),dummy,dummy,dummy
                  elseif(iso_choose.eq.1)then
                      write(4,'(8e12.4)')dummy,
     &	          vpcoeff(index),vpcoeff(index),vscoeff(index),vscoeff(index),
     &	          dummy,dummy,dummy
                  ! Note that this means, I perturb my vpv vph with my vsv model!         
	            elseif(iso_choose.eq.3)then
                      write(4,'(8e12.4)')dummy,
     &	          vpcoeff(index),vscoeff(index),vpcoeff(index),vscoeff(index),
     &	          dummy,dummy,dummy
	            endif
               endif
            enddo
	  close(4)
        enddo
	enddo
	end

c ////////////////////////////////////////////////////////////////////
c *** finds the number of the square where (lat, lon) is *************
	function superisqre(lat,lon,nsqrs,nsqtot,nlatzones,n,eq_incr)
	dimension nsqrs(nlatzones),nsqtot(nlatzones)
	incr=eq_incr*100.
	lazone=(9000-lat)/incr+1
	if(lazone.gt.nlatzones)lazone=nlatzones
	llon=lon
	if(llon.lt.0)llon=36000+llon
	superisqre=(llon*nsqrs(lazone))/36000+1
	superisqre=superisqre+nsqtot(lazone)
	if(superisqre.gt.n)superisqre=n
	return
	end



c ////////////////////////////////////////////////////////////////////
c ****** parameterization (regular grid only) *** OUTDATED VERSION!!!!
        subroutine param(eq_incr,nsqrs,nsqtot,nlatzones,numto,iswit,
     &compgrid,nlatzomax)
c---find vectors nsqrs and nsqtot that define a block parameterization
c---eq_incr,iswit,refgrid are input, the rest is output
        integer nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	parameter(pi=3.1415926536)

        numto=0
	colat=-eq_incr/2.



	do k=1,nlatzones
	   colat=colat+eq_incr
	   theta=(colat/180.)*pi
c--for this latitudinal zone, compute number of blocks (nsqrs)
	   deltalon=eq_incr/(sin(theta))
	   nsqrs(k)=(360./deltalon)+1


c--if requested, correct nsqrs(k) so the grid is compatible to reference grid
	   if(iswit.eq.1)then
              if(360./nsqrs(k).ge.compgrid)then
 100             if((mod(360./nsqrs(k),compgrid).ne.0).or.
     &              (mod(nsqrs(k),2).ne.0))then
                    nsqrs(k)=nsqrs(k)+1
                    goto100
                 else
                 endif
              elseif(360./nsqrs(k).lt.compgrid)then
 101             if((mod(compgrid,360./nsqrs(k)).ne.0).or.
     &              (mod(nsqrs(k),2).ne.0))then
                    nsqrs(k)=nsqrs(k)-1
                    goto101
                 else
                 endif
              endif
         else ! modified by lud
              if(mod(nsqrs(k),2).ne.0)nsqrs(k)=nsqrs(k)-1
	   endif
         if(MOD(NSQRS(K),2).ne.0)then
              stop "nsqrs has to be even"
         endif


c---------------------------
	   nsqtot(k)=numto
	   numto=numto+nsqrs(k)
	enddo
	nsqtot(nlatzones+1)=numto


	return
	end
