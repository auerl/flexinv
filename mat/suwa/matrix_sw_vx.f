c
c
c

      program matrix_sw_vx
c
c
c DESCRIPTION:
c computes regular grid surface wave matrices. this version was heavily 
c amended, rewritten and reorganized by ludwig auer bw 2012 and 2014
c
c CHANGELOG: 
c . version slightly amended for ETL database (lb July 2009) 
c . a comment below says this holds only for isotropic starting model?
c . version reads ker_actu from input file and read in the folder 
c . compgrid added as a hook to get a 3 deg grid upon a 2 deg ref grid
c . which is a situation the code was not amended for, before
c . information to make everything more flexible
c . added the possibility to output hitcount files
c . reads format/unit of data, phase delay (sec) or phase vel
c . allows for adaptive parameterization
c . allows for up to 4 inversion parameters (la may, 2014)
c . removed remnant stuff related to computing ATA


c
c *** some hardcoded parameters
c
      parameter(num_blocks_max=31808)
      parameter(iswit=1,refgrid=2.,adaptgrid=5.)    ! set iswit to 1 for grid to be with compgrid
      parameter(nlatzomax=180.)                     ! nr of latitudinal zones, 180 = 1 degree parametrization
      parameter(m=2000000)                          ! maximal number of datums to read /just for array dims
      parameter(max_ker=30)                         ! corresponds to 100 km equally thick layers
      parameter(max_par=8)                          ! maximal nr of inversion parameters
      parameter(max_per=16)                         ! maximal number of periods
      parameter(max_lat=180./refgrid)               ! maximal latitude for kernel grid 
      parameter(max_lon=2*max_lat)                  ! maximal nr of longitudes for kernel grid
      parameter(maxover=6,max_over=6)               ! availabl number of overtones in dataset
      parameter(ntopiso=100)                        ! index of layer at and below which isotropy is imposed
      parameter(maxinvpar=4)                        ! Maximum number of inversion parameters 4=vsh,vsv,vph,vpv
      parameter(itotnmax=num_blocks_max*maxinvpar*max_ker)  ! dimension of ata and atd arrays
      parameter(pi=3.1415926536)                    ! Pi

c
c *** allocate array memory
c
      dimension pargrid(max_lat,max_lon,max_ker,max_par,max_per,0:maxover)
      dimension simple_pargr(max_lat,max_lon,max_ker,maxinvpar,max_per,0:maxover)
      dimension iperiod(max_per,2,0:maxover),c0(max_per,2,0:maxover)
      dimension pvelgrid(max_lat,max_lon,max_per,0:maxover)
      dimension xi_in(max_ker,max_lat,max_lon) ! ratio deltavsv/deltavsh
      dimension true_per(max_per,2,0:maxover)
      dimension narcs(max_per,2,0:maxover),nperiods(2,0:maxover),nover(2)
      dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1) ! bug removed 12.04.13 la
      dimension cor3d(num_blocks_max)

      character*3 miter,perio,overt3
      character*1 typ(2),wave_type,chiarc
      character txtfile*120,perio4*4,txtfileconv*120
      character*120 dat_dir,wrk_dir,obs_fil,ker_dir
      character*10 ident
     
      dimension invparlist(5)
      integer numinvvar
      integer dataunit
c
c *** adaptive parameterization
c
      integer n1layer0,n1layer,nlatzones0
      dimension nsqrs0(nlatzomax),nsqtot0(nlatzomax+1)
      integer,dimension(4,itotnmax) :: angle_hitcount
      real*8, dimension(4,itotnmax) :: angle_dws
      integer,  dimension(itotnmax) :: hitcount
      real*8,   dimension(itotnmax) :: dws
      integer*8		            :: i

c
c *** types of data
c
      typ(1)='R' ! Rayleigh
      typ(2)='L' ! Love

c
c *** read information from input file
c
      print*,"Enter data directory"
      read(*,*),dat_dir
      print*,"Enter working directory"
      read(*,*),wrk_dir
      print*,"Enter kernel directory"
      read(*,*),ker_dir
      print*,"Enter obs modes list file"
      read(*,*),obs_fil
      print*,"RL (1), L (2), or R (3)"
      read*,irorl
c
c *** added by ludwig to make selection of inversion parameters
c *** more general, now can use b/w 1 and 4 inv parameters
c *** la, may 2014
c
      write(*,"('Number of velocity vectors to invert [max=4]:')")
      read(*,*) numinvvar
      if(numinvvar.gt.5.or.numinvvar.lt.1) stop 'Unreasonable number of velocities!'
      write(*,"('***List: 1=Vph,  2=Vpv,  3=Vsh,  4=Vsv')")
      do i=1,numinvvar
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

      print*,"Enter data unit"
      read(*,*),dataunit
      print*,"Enter horizontal pixel size eq_incr"
      read*,eq_incr
      print*,"Enter if adaptive (1) or regular (0) grid"
      read*,iadapt
      print*,"Enter number of layers ker_actu"
      read*,ker_actu
      print*,"Enter nominal number of blocks num_blocks"
      read*,num_blocks
      print*,"Enter number of iteration"
      read*,iterationnr
      write(miter,'(i3.3)')iterationnr
      print*,"Enter identifier for sw"
      read*,ident
      write(miter,'(i3.3)')iterationnr

c
c *** which datatypes to read
c
      if(irorl.eq.1)then ! Love and Rayleigh
         nrl1=1
         nrl2=2
      elseif(irorl.eq.2)then ! Love only
         nrl1=2
         nrl2=2
      elseif(irorl.eq.3)then ! Rayleigh only
         nrl1=1
         nrl2=1
      else
         print*,"Error, irorl has to be 1, 2 or 3!"
         stop
      endif

c
c *** regular parameterization
c
      nlatzones=180./eq_incr 
      itotn=num_blocks*numinvvar*ker_actu
c      ndimata=itotn*(itotn+1)/2
      compgrid=eq_incr

c
c *** adaptive parameterization
c
      hitcount=0
      angle_hitcount=0
      angle_cat=0
      angle_dws=0
      dws=0

c
c *** in case adaptive param is used
c
      if(iadapt.eq.1)then
         compgrid=adaptgrid
      endif

c
c *** horizontal parametrization
c
	print*,'Number of latitudinal zones: ',nlatzones
	if(mod(nlatzones,2).ne.0)then
	   print*,'Error! nlatzones should not be odd!'
	   stop
	endif
c        
c *** define parameterization grid
c
      call param(eq_incr,nsqrs,nsqtot,nlatzones,numto,iswit,compgrid,nlatzomax,iadapt)
      print*,numto," horizontal pixels defined"
      if(numto.ne.num_blocks)stop "Error! Wrong value of num_blocks parameter!"

c
c *** check available data
c
	do ity=nrl1,nrl2 ! loop over R, L
         ! need to do a
	   open(75,file=obs_fil(1:lnblnk(obs_fil))//typ(ity)//".txt")
	   iover0=-1
440	   read(75,*,end=442)wave_type,iover,intt
	   if(iover.gt.iover0)then
	      if(iover.gt.0)nperiods(ity,iover0)=i
	      i=1
	      iover0=iover
	   elseif(iover.eq.iover0)then
	      i=i+1
	   else
	      stop "must sort file with data info"
	   endif
	   iperiod(i,ity,iover)=intt
	   narcs(i,ity,iover)=1  ! no greater-arc data available
	   goto440
442	   continue
	   nover(ity)=iover
	   nperiods(ity,iover0)=i
	enddo ! end of loop over R,L
	write(*,*)nover(1),"R Ts:",(nperiods(1,iover),iover=0,maxover)
	write(*,*)nover(2),"L Ts:",(nperiods(2,iover),iover=0,maxover)

c
c *** loop to check the exact periods for which we have measurements
c
	print*,"Checking periods of available data" 
	do ity=1,2
         do iover=0,nover(ity)
            write(overt3,"(i3.3)")iover
            do k=1,nperiods(ity,iover)
               write(perio,"(i3.3)")iperiod(k,ity,iover)
               write(perio4,"(i4.4)")iperiod(k,ity,iover)
               do iarc=1,narcs(k,ity,iover)
                  write(chiarc,'(i1)')iarc
                  txtfile=dat_dir(1:lnblnk(dat_dir))//"/summary/summ."//typ(ity)//"."//overt3//"."//perio//".txt"
                  open(2,file=txtfile,status="old")
                  read(2,*)wave_type,true_per(k,ity,iover),c0(k,ity,iover)
                  if(wave_type.ne.typ(ity))stop "wrong data file"
                  close(2)
               enddo
            enddo
         enddo
	enddo
	jtot=0


c
c *** loop over L, R
c

c
c *** simple_pargr array is a "simplified" version of the full kernel array pargr
c *** which is reduced to the selected inversion parameters
c

	do ity=nrl1,nrl2
	  call load_kernels(pargrid,max_lat,max_lon,max_ker,max_par,! calcuate kernels for all overtones
     &	          max_per,miter,typ(ity),nperiods,pvelgrid,true_per,c0,
     &                dat_dir,wrk_dir,obs_fil,ker_dir)
	  ntola=180./refgrid
	  ntolo=360./refgrid
	  do ilat=1,ntola
           do ilon=1,ntolo
              do iover=0,max(nover(ity),maxover)
                 do iper=1,nperiods(ity,iover)
                    do iker=1,ker_actu
                       ! In case of a parameterization in vsh, vsv, vpv and vph i can 
                       ! freely combine these parameters. In case of a Xi-Voigt parameterization
                       ! this is done already at an earlier stage (the kernel making stage).                       
                       ! in case of a xi voigt parameterization invparlist should just 
                       ! obtain the right (4 and 5) parameter indices from setup.py

                       do iniv=1,numinvvar
                          ! increase by one since the parameter "index" is different in the 
                          ! body wave and the surface wave version, in the surface wave version
                          ipar=invparlist(iniv)+1
                          simple_pargr(ilat,ilon,iker,iniv,iper,iover)=
     &	                  pargrid(ilat,ilon,iker,ipar,iper,iover)
                       enddo

c
c                        ! Write out phase velocity maps for diagnosis
c	                    if(iper.eq.1)write(199,"(3(i4,1x),2(e16.8,2x))")iker,ilon,ilat,
c                       &simple_pargr(ilat,ilon,iker,1,iper,iover),simple_pargr(ilat,ilon,iker,2,iper,iover)
c	                    if(iper.eq.1.and.iker.eq.1)
c                       &write(299,*)ilon,ilat,pvelgrid(ilat,ilon,iper,iover)                        
c
                     enddo
                  enddo
               enddo
            enddo
         enddo
c
c *** finally, build the matrix and store it in sparse matrix CRS format
c *** and save hitcount files associated to each matrix
c
	  do nmod=0,max_over ! Calculate matrices for some overtones
            do iper=1,nperiods(ity,nmod)
               do iarc=1,narcs(iper,ity,nmod)
                  jold=jtot
                  datamean=0.
                  write(overt3,"(i3.3)")nmod
                  write(perio,"(i3.3)")iperiod(iper,ity,nmod)
                  ! Once per mode
	               open(77,file=wrk_dir(1:lnblnk(wrk_dir))//"/sw.xxx."//ident(1:lnblnk(ident))//"."
     &               //typ(ity)//"."//overt3//"."//perio//".d",access='direct',recl=4,form='unformatted')
	               open(99,file=wrk_dir(1:lnblnk(wrk_dir))//"/sw.indx."//ident(1:lnblnk(ident))//"."
     &               //typ(ity)//"."//overt3//"."//perio//".d",access='direct',recl=4,form='unformatted')
	               open(66,file=wrk_dir(1:lnblnk(wrk_dir))//"/sw.rhs."//ident(1:lnblnk(ident))//"."
     &               //typ(ity)//"."//overt3//"."//perio//".d")
	               open(88,file=wrk_dir(1:lnblnk(wrk_dir))//"/sw.poin."//ident(1:lnblnk(ident))//"."
     &               //typ(ity)//"."//overt3//"."//perio//".d")

                  ! Initialize nr of measurement
                  nnn=0
                  ! Build A/ATA matrix


                  call build_matrix(hitcount,numto,nsqrs,
     &	          nsqtot,nlatzones,eq_incr,m,pi,refgrid,
     &	          iperiod(iper,ity,nmod),typ(ity),pvelgrid,
     &	          max_lat,max_lon,max_ker,max_par,max_per,
     &            iper,simple_pargr,jtot,itotn,
     &            ker_actu,iaorata,nnn,cor3d,num_blocks,
     &            datamean,iarc,nmod,maxover,itotnmax,
     &            num_blocks_max,nlatzomax,
     &            dat_dir,wrk_dir,obs_fil,ker_dir,
     &            angle_hitcount,dws,angle_dws,
     &            dataunit,numinvvar,invparlist,
     &            maxinvpar)

                  ! Some verbose
                  print*,datamean,datamean/(Jtot-jold)

                  ! Close files
                  close(66)
                  close(77)
                  close(88)
                  close(99)
c
c *** i just use the hitcount and not the dws, from now on
c

                  open(unit=30,file=wrk_dir(1:lnblnk(wrk_dir))//'dws/dws.'//
     &ident(1:lnblnk(ident))//'.'//typ(ity)//'.'//overt3//'.'//perio//'.d')

                  do i=1,itotn
                     write(30,'(4I15)') angle_hitcount(1:4,i)
                     ! flush hitcounts
                     hitcount(i)=0.
                     angle_hitcount(1:4,i)=0.
                     dws(i)=0.
                     angle_dws(1:4,i)=0
                  enddo
                  close(30)

                  ! Some verbose
                  print*,'period ',iper,' : ',Jtot-jold,' data included'
                  print*,typ(ity),iperiod(iper,ity,nmod),Jtot-jold,(datamean/(Jtot-jold))

               enddo
            enddo
         enddo
	enddo !end of loop over L,R


        end program matrix_sw_vx






c
c =====================================================================
c  subroutine build_matrix
c  subroutine that actually augments the kernel matrix by combining
c  measurements and preloaded kernels
c
	subroutine build_matrix(hitcount,n1lay,nsqrs,
     &	nsqtot,nlatzones,eq_incr,m,pi,refgrid,iperiod,type,pvelgrid,
     &	max_lat,max_lon,max_ker,max_par,max_per,iper,simple_pargr,
     &	j,itotn,ker_actu,iaorata,nnn,cor3d,num_blocks,datamean,
     &      iarc,iover,maxover,itotnmax,num_blocks_max,nlatzomax,
     &      dat_dir,wrk_dir,obs_fil,ker_dir,angle_hitcount,dws,angle_dws,
     &      dataunit,numinvvar,invparlist,maxinvpar)

	parameter(rad=0.0174533) ! rad=pi/180.
        parameter(ndim=10000)
        parameter(conv_to_phdel=0) ! Small mode to be used to convert from
                                   ! Phase velocity to phase delay in s

	dimension pvelgrid(max_lat,max_lon,max_per,0:maxover)
	dimension simple_pargr(max_lat,max_lon,max_ker,maxinvpar,max_per,0:maxover)
	dimension cor3d(num_blocks_max)
	dimension jsqrg(ndim),ddelg(ndim),jsqrm(ndim),ddelm(ndim)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax)
        dimension invparlist(5) 

        integer dataunit,numinvvar,n1lay
	character*1 wave_type,type,chiarc
	character perio4*4,perio3*3,overt3*3
        character*120 dat_dir,wrk_dir,obs_fil
	character*120 txtfile,txtfileconv*120

        ! stuff related to adaptive parameterization
        integer*8 nnn
        real*8 eplat,eplon,stlat,stlon
        real*8 xlat_new,xlat_old,xlon_old,xlon_new
        real*8 deltax,az_old,az_new
        integer angle_cat,nperlay
        integer,dimension(4,itotnmax) :: angle_hitcount
        integer,dimension(itotnmax) :: hitcount
        real*8,dimension(itotnmax) :: dws
        real*8,dimension(4,itotnmax) :: angle_dws
        real*8,dimension(ndim) :: xlonm,xlong,xlatm,xlatg   





	print*,'Starting routine buildmatrix. Augmenting the matrix with ',type,iperiod,iover,' data'

        nnn=0.

        ! LA, May 2014: when i increased to a maximum of 4 inversion parameters, a weird bug
        ! occured, n1lay is set to zero for some reason, the following is just a "hook". 
        ! Ned to further investigate that. Seems to work at least...

        nperlay=n1lay ! la, may 2014: 

        print*,n1lay,"------------------------"
c        stop


        ! concenate names of data file
	write(chiarc,'(i1)')iarc
	write(perio3,"(i3.3)")iperiod ! this routine is called only once per mode so iperiod can be scalar
	write(perio4,"(i4.4)")iperiod
	write(overt3,"(i3.3)")iover



	txtfile=dat_dir(1:lnblnk(dat_dir))//"/summary/summ."//type//"."//overt3//"."//perio3//".txt"



        ! open datafile and read header information
	open(2,file=txtfile)
	read(2,*)wave_type,period,c0 ! read header information
	if(type.ne.wave_type)stop "mismatch toroidal vs spheroidal mode"
	print*,wave_type,period,' (',iper,') ',c0
	read(2,*) ! read ### line dummy

	print*,"Start reading phase delay data"

        if(conv_to_phdel.eq.1)then
           txtfileconv="/tmp/summary/summ."//type//"."//overt3//"."//perio3//".txt"
           open(43,file=txtfileconv)
           write(43,1003)wave_type,period,c0
           write(43,'(a46)') '# columns: eqlat eqlon stlat stlon dt error'
           print*,"CONVERTING TO PHASE DELAY IN SECONDS"
        endif




c ********************************************************
c *** Major GOTO loop over each line in a summary file ***
c ********************************************************
23	j=j+1
	read(2,*,end=999)eplat,eplon,stlat,stlon,dphi_k,error

        ! determine great circle path
	call ocav(nsqrs,nsqtot,nlatzones,eq_incr,n,
     &	eplat,eplon,stlat,stlon,ientm,jsqrm,ddelm,ient,
     &	jsqrg,ddelg,pi,delta,xlatm,xlonm,xlatg,xlong)

        ! don't know what happens here
	if(jsqrg(ient).eq.jsqrg(1))then
	   ddelg(1)=ddelg(1)+ddelg(ient)
	   ient=ient-1
	endif

        ! unit conversion operations (is sometimes to be amended)
        if(dataunit.eq.1)then        
           delphi=dphi_k ! have already phase delay in sec, eq etl data
        elseif(dataunit.eq.2)then 
           t0=(delta*6371.)/c0 
           delphi=-dphi_k*t0 ! get absolute phase delay in sec, eg vtk, hvh data
        endif
        
        ! if requested, write out converted data
        if(conv_to_phdel.eq.1)then
           write(43,1002)eplat,eplon,stlat,stlon,delphi,error
           if(mod(j,1000).eq.0) print*,j,' data converted'
           goto 23
1002       format(4(f11.4,1x),1x,f11.6,1x,f11.7)
1003       format(8x,a1,7x,f8.4,10x,f8.5)
        endif

        ! this is probably obsolete, but im afread to delete it
	corr_to_delph=0.

        ! Generalized to R3, L3, R4, L4, etc.:
	iarx=iarc

c
c **** Major raypath loop
c **** Major loop over arcs of raypath
c

	do iturns=1,max(iarc-1,1)
         if(iarx.eq.1)then
            iloop1=1
            iloop2=ientm
         elseif(iarx.eq.2)then
            iloop1=ientm+1
            iloop2=ient
         elseif(iarx.gt.2)then ! Orbits >2 not yet implemented
            print*,'Orbits >2 not implemented, error!'
            stop
            iloop=ient
            iarx=iarx-2
         endif

         ! Augment the A matrix for all points over the raypath
         ! Major DO loop over all points along the raypath
         do i=iloop1,iloop2

            ! reference velity and kernel of current pixel
            call coordsuper(jsqrg(i),blocla,bloclo,nsqrs,nlatzones,eq_incr)

            ! find azimuth of raypath, written by julia
            xlat_new=xlatg(i)
            xlon_new=xlong(i)
            if (xlon_new.lt.0.0_8) xlon_new = xlon_new + 360.0_8              
            if (i/=iloop1) then
               xlat_old=xlatg(i-1)
               xlon_old=xlong(i-1)
            end if
            if (xlon_old.lt.0.0_8) xlon_old = xlon_old + 360.0_8               
            call delazs(xlat_old,xlon_old,xlat_new,xlon_new,deltax,az_old,az_new)
            if (az_new.gt.360.0_8) az_new = az_new - 360.0_8
            if (az_new.gt.180.0_8) az_new = az_new - 180.0_8

            ! converts blocla to kernel domain index ilat
            ilat=(90.-blocla)/refgrid+1  
            ilon=bloclo/refgrid+1        ! same as above for longitude

            ! Stopping criteria
            if((ilat.gt.max_lat).or.(ilon.gt.max_lon))then
               print*,"Problem at location",ilat,ilon,max_lat,max_lon
               stop
            endif

            ! Increment the correction to delphi for this block
            corr_to_delph=corr_to_delph+(ddelg(i)*6371.*pvelgrid(ilat,ilon,iper,iover)) ! pvelgrid is d_slowness(km/s)
            cor3d(jsqrg(i))=-pvelgrid(ilat,ilon,iper,iover)*c0*100. ! So cor3d is percent d_velocity

            ! Augment the A matrix for current point along raypath
            do i1=1,numinvvar ! loop over vsv and vsh .. and so on
               do i2=1,ker_actu ! loop over layers
                  rad_ker=simple_pargr(ilat,ilon,i2,i1,iper,iover)
                  if(abs(rad_ker*(ddelg(i)*6371.)).gt.(1.e-5))then ! Warning: to be changed if units are changed, might also be (1.e-8)
                     
                     write(77,rec=nnn+1)rad_ker*(ddelg(i)*6371.)
                     write(99,rec=nnn+1)jsqrg(i)+(nperlay*ker_actu*(i1-1))+(nperlay*(i2-1))                       
                     
                     ! define hitcounts, dws and angle dws
                     jj=jsqrg(i)+(nperlay*ker_actu*(i1-1))+(nperlay*(i2-1))
                     hitcount(jj)=hitcount(jj)+1     ! add to hitcount value      
                     dws(jj)=dws(jj) + abs( rad_ker*(ddelg(i)*6371.) )
                     if(az_new>=0.and.az_new<45) then
                        angle_cat=1
                        angle_hitcount(angle_cat,jj)=angle_hitcount(angle_cat,jj)+1
                        angle_dws(angle_cat,jj)=angle_dws(angle_cat,jj)+ abs( rad_ker*(ddelg(i)*6371.) )
                     elseif(az_new>=45.and.az_new<90) then
                        angle_cat=2
                        angle_hitcount(angle_cat,jj)=angle_hitcount(angle_cat,jj)+1
                        angle_dws(angle_cat,jj)=angle_dws(angle_cat,jj)+ abs( rad_ker*(ddelg(i)*6371.) )
                     elseif(az_new>=90.and.az_new<135) then
                        angle_cat=3
                        angle_hitcount(angle_cat,jj)=angle_hitcount(angle_cat,jj)+1
                        angle_dws(angle_cat,jj)=angle_dws(angle_cat,jj)+ abs( rad_ker*(ddelg(i)*6371.) )
                     elseif(az_new>=135.and.az_new<180) then
                        angle_cat=4
                        angle_hitcount(angle_cat,jj)=angle_hitcount(angle_cat,jj)+1
                        angle_dws(angle_cat,jj)=angle_dws(angle_cat,jj)+ abs( rad_ker*(ddelg(i)*6371.) )
                     else	
                        print*,"angle definition didn't work:"
                        print*,"ray data:", eplat, eplon, stlat, stlon !, nsqrs, nsqtot, nlatzones, numto
                        print*,xlat_old, xlon_old, xlat_new, xlon_new, az_new, az_old
                     endif
                     
                     ! Increment nr of matrix entry
                     nnn=nnn+1 
                  endif
               enddo
            enddo
         enddo ! end of major loop over all points along raypath
      enddo ! end of major loop over arcs of ray
      
      ! correct the current datum for 3-D starting model (pvelgrid):
      delphi=delphi-corr_to_delph
      
      ! Store the d vector and the pointer on disc
      write(88,*)nnn
      write(66,*)delphi
      
      ! for diagnostics
      datamean=datamean+delphi
      
      ! verbose, nr of data
      if(mod(j,1000).eq.0)then
         print*,j,' data read'
      endif
      goto 23
c ************** end of major GOTO loop *****************
      ! some more verbose

999   print*,'Number of paths used:',j-1
      j=j-1
      nnntot=nnn
      if(conv_to_phdel.eq.1)close(43)           
      return
      end subroutine build_matrix ! end of subroutine buildmatrix


c
c =====================================================================
c subroutine isqre
c finds the number of the square where (xlat,xlon) is
c
	function isqre(xlat,xlon,nsqrs,nsqtot,nlatzones,n,eq_incr)
	parameter(nlatzomax=180,nlatzhmax=720)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax+1)
	dimension nsqrsh(nlatzhmax),nsqtoth(nlatzhmax+1)
        lazone=(90.-xlat)/eq_incr+1
	if(lazone.gt.nlatzones)lazone=nlatzones
	isqre=(xlon/360.)*nsqrs(lazone)+1
	isqre=isqre+nsqtot(lazone)
	return
        end function isqre ! end of subroutine isqre






c
c =====================================================================
c  SUBROUTINE RANGE
c  FINDS THE COORDINATE RANGE OF SQUARE NUMBER 'NSQ'
c
	SUBROUTINE RANGE(NSQ,XLAMIN,XLAMAX,XLAMID,XLOMIN,XLOMAX,
     &	XLOMID,nsqrs,nsqtot,nlatzones,n,eq_incr)
        PARAMETER(NLATZOMAX=180.)
	DIMENSION NSQRS(NLATZONES),NSQTOT(NLATZONES)
	LAZONE=2
	DO WHILE (NSQ.GT.NSQTOT(LAZONE))
	   LAZONE=LAZONE+1
	ENDDO
	LAZONE=LAZONE-1
	NNSQ=NSQ-NSQTOT(LAZONE)
	XLAMIN=90.-LAZONE*eq_incr
	XLAMAX=XLAMIN+eq_incr
	XLAMID=XLAMIN+(eq_incr*0.5)
	GRSIZE=360./NSQRS(LAZONE)
	XLOMAX=NNSQ*GRSIZE
	XLOMIN=XLOMAX-GRSIZE
	XLOMID=XLOMAX-(GRSIZE*0.5)
	RETURN
        END SUBROUTINE RANGE ! end of subroutine RANGE





c
c =====================================================================
c  subroutine coordsuper
c  given cell index on the surface, finds lon and lat of its center
c
	subroutine coordsuper(nbloc,blocla,bloclo,nsqrs,nlatzones,eq_incr)
        parameter(nlatzomax=180.)
	dimension nsqrs(nlatzomax)
	ntot=0
        ! Loop(s) over all the blocks
	do 500 ila=1,nlatzones
           ! Increment latitude
	   rlati=90.-(eq_incr*(ila-1))
           ! Calculate increment in longitude for this band
	   rinlo=(360./nsqrs(ila))
	   do 400 isq=1,nsqrs(ila)
	      rlong=(360./nsqrs(ila))*(isq-1)
	      ntot=ntot+1
	      if(ntot.eq.nbloc)then
	         bloclo=rlong+(rinlo/2.)
	         blocla=rlati-(eq_incr/2.)
	         goto 600
	      endif
400	   continue
500	continue
600	return
	end subroutine coordsuper ! end of subroutine coordsuper





c
c =====================================================================
c  subroutine rsoinc
c
      SUBROUTINE RSOINC(A,N,IDX)
      DIMENSION A(1),IDX(1) 
      IF (N.EQ.1) GO TO 65
      IF (N.LE.0) GO TO 60
      DO 1 I = 1,N
      IDX(I) = I
    1 CONTINUE
      N2 = N/2
      N21 = N2 + 2
      ICT=1 
      I=2 
   11 N1=N21-I
      NN=N
      IK=N1 
   15 C=A(IK) 
      IC=IDX(IK)
  100 JK=2*IK 
      IF (JK.GT.NN) GO TO 140 
      IF (JK.EQ.NN) GO TO 120 
       IF (A(JK+1).LE.A(JK)) GO TO 120
      JK=JK+1 
  120 IF (A(JK).LE. C) GO TO 140
      A(IK)=A(JK) 
      IDX(IK)=IDX(JK) 
      IK=JK 
      GO TO 100 
  140 A(IK)=C 
      IDX(IK)=IC
      GO TO (3,45) ,ICT 
    3 IF (I.GE.N2) GO TO 35 
      I=I+1 
      GO TO 11
   35 ICT=2 
      NP2=N+2 
      I=2 
   37 N1=NP2-I
      NN=N1 
      IK=1
      GO TO 15
  45  CONTINUE
      T = A(1)
      A(1) = A(N1)
      A(N1) = T 
      IT = IDX(1) 
      IDX(1) = IDX(N1)
      IDX(N1) = IT
      IF (I.GE.N) GO TO 55
      I=I+1 
      GO TO 37
   55 RETURN
   60 WRITE(16,500)
  500 FORMAT('ERROR RETURN FROM SORTD1 - N LESS THAN OR EQUAL TO 1')
      STOP
   65 IDX(1)=1
      RETURN
      END SUBROUTINE RSOINC





c
c =====================================================================
c  subroutine ocav
c  computes great circle path of measurment 
c
	subroutine ocav(nsqrs,nsqtot,nlatzones,eq_incr,n,
     &	eplat,eplon,stlat,stlon,ientm,jsqrm,ddelm,
     &	ient,jsqrg,ddelg,pi,del,xlatm,xlonm,xlatg,xlong)
        parameter (ndim=10000)
        parameter(nlatzomax=180)
	dimension geo_colat(ndim),geo_long(ndim)
	dimension pht(ndim),idx(ndim),jsqrg(ndim),ddelg(ndim)
	dimension jsqrm(ndim),ddelm(ndim)
	dimension nsqrs(nlatzomax),nsqtot(nlatzomax)
        real*8,dimension(ndim) :: xlonm,xlong, xlatm,xlatg
        real*8 eplat, eplon, stlat, stlon     


	PI2=PI*2.
	FORPI=4.*PI
	RADIAN=180./PI
	nth1=nlatzones+1
	dth=pi/(nlatzones*1.)

        TH1=(90.-EPLAT)/RADIAN
        PH1=EPLON/RADIAN
        TH2=(90.-STLAT)/RADIAN
        PH2=STLON/RADIAN
        STH1=SIN(TH1)
        STH2=SIN(TH2)
        CTH1=COS(TH1)
        CTH2=COS(TH2)
        SPH1=SIN(PH1)
        SPH2=SIN(PH2)
        CPH1=COS(PH1)
        CPH2=COS(PH2)
        CPH21= CPH1*CPH2+SPH1*SPH2
        SPH21=SPH2*CPH1-SPH1*CPH2
        CDEL=STH1*STH2*CPH21+CTH1*CTH2
        CCAPTH=STH1*STH2*SPH21/SQRT(1.-CDEL*CDEL)
        SCAPTH=SQRT(1.-CCAPTH*CCAPTH)
        capth=atan2(scapth,ccapth)
        SCAPPH=CTH1*STH2*CPH2-CTH2*STH1*CPH1
        CCAPPH=STH1*CTH2*SPH1-STH2*CTH1*SPH2
        CAPPH=ATAN2(SCAPPH,CCAPPH)
        SCAPPH=SIN(CAPPH)
        CCAPPH=COS(CAPPH)
c
        DEL=ATAN2(SQRT(1.-CDEL*CDEL),CDEL)
        CPHSP=CCAPTH*STH1*(CPH1*CCAPPH+SPH1*SCAPPH)-SCAPTH*CTH1
        SPHSP=STH1*(SPH1*CCAPPH-CPH1*SCAPPH)
        PHSP=ATAN2(SPHSP,CPHSP)
c
        thet=capth
        if(capth.gt.0.5*pi) thet=pi-capth
        thmin=0.5*pi-thet
        thmax=0.5*pi+thet
        lat_zone=0
        IENT=0
        IF(SCAPTH.EQ.0) GOTO 10
c
        DO 20 I=2,NTH1
           lat_zone=lat_zone+1
           TH=FLOAT(I-1)*DTH
           CTH=COS(TH)
           CPHT=-CTH/SCAPTH
           CPHT2=CPHT*CPHT
           if(th.lt.thmin.or.th-dth.gt.thmax) then
              go to 20
           endif
           IF(CPHT2.GT.1.) GOTO 21
           IENT=1+IENT
           PHT(IENT)=ATAN2(SQRT(1.-CPHT2),CPHT)
           spht=sin(pht(ient))
           geo_colat(ient)=th
           geo_long(ient)=atan2(spht,-cth*ccapth/scapth)
           xel=AMOD(PHT(ient)-phsp+FORPI,PI2)
           IENT=1+IENT
           PHT(IENT)=-PHT(IENT-1)
           geo_colat(ient)=th
           geo_long(ient)=atan2(-spht,-cth*ccapth/scapth)
           if(geo_long(ient).lt.0) geo_long(ient)=pi2+geo_long(ient)
           xel=AMOD(PHT(ient)-phsp+FORPI,PI2)
  21       numlong=nsqrs(lat_zone)
           dphi=pi2/float(numlong)
           DO 40 j=1,numlong
              PH=FLOAT(j-1)*DPHi
              angr=ph-capph
              thlo=atan(-ccapth/(scapth*cos(angr)))
              if(thlo.lt.0.) thlo=pi+thlo
              if(thlo.gt.th-dth.and.thlo.lt.thmin) go to 40
              if(thlo.lt.th+dth.and.thlo.gt.thmax) go to 40
              if(thlo.gt.th.or.thlo.lt.th-dth) go to 40
              SPH=SIN(PH)
              CPH=COS(PH)
              IENT=IENT+1
              PHT(IENT)=ATAN2(CCAPTH*(SPH*CCAPPH-CPH*SCAPPH),
     1        CPH*CCAPPH+SPH*SCAPPH)
              IF(PHT(IENT).GT.PI) PHT(IENT)=PHT(IENT)-PI2
              geo_colat(ient)=thlo
              geo_long(ient)=ph
              xel=AMOD(PHT(ient)-phsp+FORPI,PI2)
   40      CONTINUE
   20   CONTINUE
   10 CONTINUE
c 
      DO 60 I=1,IENT
   60 PHT(I)=AMOD(PHT(I)-PHSP+FORPI,PI2)
      IENT=IENT+1
      PHT(IENT)=0.
      IENT=IENT+1
      PHT(IENT)=DEL
c
      CALL RSOINC(PHT,IENT,IDX)
c
      PHT(IENT+1)=PI2
      ientm=0
      ientg=0
      cdelta=0.
      DO 50 I=1,IENT
      I1=1+I
      PHTT=.5*(PHT(I)+PHT(I1))+PHSP
      CPHT=COS(PHTT)
      SPHT=SIN(PHTT)
      CTH=-CPHT*SCAPTH
      TH=ATAN2(SQRT(1.-CTH*CTH),CTH)
      CPH=CPHT*CCAPTH*CCAPPH-SPHT*SCAPPH
      SPH=CPHT*CCAPTH*SCAPPH+SPHT*CCAPPH
      IF(SPH.EQ.0..AND.CPH.EQ.0.) GOTO 9873
      PH=ATAN2(SPH,CPH)
      if(ph.lt.0.) ph=ph+pi2
      if(ph.gt.pi2) ph=ph-pi2
      GOTO 9874
 9873 PH=0.
 9874 CONTINUE
      XLAT=90.-TH*RADIAN
      ilat=xlat*100.+0.5
      XLON=PH*RADIAN
      ilon=xlon*100.+0.5
      jsqre=isqre(xlat,xlon,nsqrs,nsqtot,nlatzones,n,eq_incr) ! lapo 16/5/2006
      RD=(PHT(I1)-PHT(I))
      IF(PHT(I1).LE.DEL) then
      ientm=ientm+1
      jsqrm(ientm)=jsqre
      ddelm(ientm)=rd
      xlatm(ientm)=xlat
      xlonm(ientm)=xlon
c
c some type of correction
c
      if(ientm.gt.1.and.jsqrm(ientm).eq.jsqrm(ientm-1)) then
      ientm=ientm-1
      ddelm(ientm)=ddelm(ientm)+rd
      endif
c
      cdelta=cdelta+rd
      call range(jsqre,XLAMIN,XLAMAX,blocla,
     &	XLOMIN,XLOMAX,bloclo,nsqrs,nsqtot,nlatzones,n,eq_incr)
      endif
c
      ientg=ientg+1
      jsqrg(ientg)=jsqre
      ddelg(ientg)=rd
      xlong(ientg)=xlon
      xlatg(ientg)=xlat
      if(ientg.gt.1.and.jsqrg(ientg).eq.jsqrg(ientg-1)) then
      ientg=ientg-1
      ddelg(ientg)=rd+ddelg(ientg)
      endif
   50 continue
      RETURN
      END SUBROUTINE OCAV! end of subroutine ocav






c
c =====================================================================
c  Subroutine load_kernels
c
c  delivers integrated kernels, stored in the array pargrid, in the
c  order: rho, vph, vpv, vsh, vsv, eta, Q_mu, Q_kappa (4-th index
c  of pargrid, from 1 to max_par=8)
c
      subroutine load_kernels(pargrid,max_lat,max_lon,max_ker,
     &   	max_par,max_per,miter,type,nperiods,pvelgrid,true_per,c0,
     &          dat_dir,wrk_dir,obs_fil,ker_dir)

      parameter (maxover=6) ! highest overtone to be computed 
      parameter (twopi=6.2831853072)
      parameter (maxker=100)
      parameter (maxpar=8)
      parameter (maxmod=100)
      
      dimension pargrid(max_lat,max_lon,max_ker,max_par,max_per,0:maxover)
      dimension true_per(max_per,2,0:maxover),c0(max_per,2,0:maxover),nperiods(2,0:maxover)
      dimension pvelgrid(max_lat,max_lon,max_per,0:maxover)
      dimension parms(maxmod,maxker,maxpar,0:maxover)
      dimension parmt(maxmod,maxker,maxpar,0:maxover)
      dimension omegs(maxmod,0:maxover)
      dimension omegt(maxmod,0:maxover)
      dimension omegain(maxmod)
      dimension gvels(maxmod,0:maxover)
      dimension gvelt(maxmod,0:maxover)
      dimension pvels(maxmod,0:maxover)
      dimension pvelt(maxmod,0:maxover)
      dimension nmods(0:maxover),nmodt(0:maxover)
      dimension array(maxmod)
      dimension arraycube(3,maxmod)
      dimension work(3,maxmod)

      character*120 string
      character*3 miter ! iteration index or index of starting model
      character*120 kernelfile
      character*1 tmod
      character*1 type
      character*120 dat_dir,wrk_dir,obs_fil,ker_dir

      print*,'Buffering the ',type,' kernels'
      if(type.ne.'L'.and.type.ne.'R')stop 'check Love or Rayleigh'
c
c *** loop over all reference grid locations
c
      do ilat=1,max_lat
         do ilon=1,max_lon
            ! Concenate kernelfile name
            write(kernelfile,"('/L_',i3.3,i3.3,'_',a3,'_kernels')")ilat,ilon,miter
            ! Open kernelfile
            open(1,file=ker_dir(1:lnblnk(ker_dir))//kernelfile,iostat=ios,status="old")
	    do nmod=0,maxover
               nmods(nmod)=0
               nmodt(nmod)=0
            enddo
            if(ios.eq.0) then ! ? whats ios ?
               do while (ios.eq.0) ! we keep reading as long as there are lines in the kernel file
                  read(1,"(a)",iostat=ios) string ! first read header for each mode
                  if(ios.eq.0) then
                     read(string,"(1x,i2,1x,a1,i4,2e15.7,i3,i3)",iostat=ios)
     &               nmod,tmod,lmod,ommod,gvmod,nker,npar ! nmod is overtone number
                  endif
                  if(ios.eq.0) then
                     if(tmod.eq.'T') then
                        lora=1
                        nmodt(nmod)=nmodt(nmod)+1
                        omegt(nmodt,nmod)=ommod
                        gvelt(nmodt,nmod)=gvmod
                        pvelt(nmodt,nmod)=ommod*6371./(float(lmod)+0.5)
                        factor=(pvelt(nmodt(nmod),nmod)**2)/(gvmod*ommod)
                     else if(tmod.eq.'S') then
                        lora=2
                        nmods(nmod)=nmods(nmod)+1
                        omegs(nmods,nmod)=ommod
                        gvels(nmods,nmod)=gvmod
                        pvels(nmods,nmod)=ommod*6371./(float(lmod)+0.5)
                        factor=(pvels(nmods(nmod),nmod)**2)/(gvmod*ommod)
                     else
                        stop 'check mode type'
                     endif
                     ! increment arrays parms and parmt
                     do iker=1,nker !loop over integrated kernels
                        if(lora.eq.1) then
                           read(1,"(3x,8e12.4)") (parmt(nmodt(nmod),iker,ipar,nmod),ipar=1,npar)
                           do ipar=1,npar
                              parmt(nmodt(nmod),iker,ipar,nmod)=parmt(nmodt(nmod),iker,ipar,nmod)*factor
                           enddo
                        else if(lora.eq.2) then
                           read(1,"(3x,8e12.4)") (parms(nmods(nmod),iker,ipar,nmod),ipar=1,npar)
                           do ipar=1,npar
                              parms(nmods(nmod),iker,ipar,nmod)=parms(nmods(nmod),iker,ipar,nmod)*factor
                           enddo
                        endif
                     enddo !end of loop over integrated kernels
                  endif ! end of if over ios
               enddo ! end of while loop (ios.eq.0) that reads kernel file for this location
               close(1)

               if(type.eq.'L') then
                  ! loop over all L modes for which we have kernels
                  do nmod=0,maxover! loop over all branches
                     ! spline all parameters and store in gridded array
                     do iker=1,nker
                        do ipar=1,npar
                           do imod=1,nmodt(nmod) !copy relevant entries of parmt and omegt to working vectors
                              array(imod)=parmt(imod,iker,ipar,nmod)
                              omegain(imod)=omegt(imod,nmod)
                           enddo
                           call rspln(1,nmodt(nmod),omegain,array,arraycube,work)
                           do iper=1,nperiods(2,nmod)
                              omega=twopi/true_per(iper,2,nmod)
                              ! modes are like nodes between which we must interpolate
                              value=rsple(1,nmodt(nmod),omegain,array,arraycube,omega) 
                              pargrid(ilat,ilon,iker,ipar,iper,nmod)=value
                           enddo
                        enddo
                     enddo
                     ! spline phase velocity perturbation and store in gridded array (no need to loop over ker and par)
                     do imod=1,nmodt(nmod) ! copy relevant entries of parmt and omegt to working vectors
                        array(imod)=pvelt(imod,nmod)
                        omegain(imod)=omegt(imod,nmod)
                     enddo
                     call rspln(1,nmodt,omegain,array,arraycube,work)
                     do iper=1,nperiods(2,nmod)
                        omega=twopi/true_per(iper,2,nmod)
                        value=rsple(1,nmodt,omegain,array,arraycube,omega)
                        valueref=c0(iper,2,nmod)! here i used to call premgephvelo
                        ! pvelgrid: absolute difference in slowness between our 3d reference and prem models:
                        pvelgrid(ilat,ilon,iper,nmod)=-(value-valueref)/(valueref*valueref)
                        ! correct pargrid: we want slowness not velocity kernels:
                        do iker=1,nker
                           do ipar=1,npar
                              pargrid(ilat,ilon,iker,ipar,iper,nmod)=-pargrid(ilat,ilon,iker,ipar,iper,nmod)/(value*value)
                           enddo
                        enddo
                     enddo ! end of loop over periods
                  enddo ! end of loop over branches
                  ! loop over all L modes for which we have kernels
               else if(type.eq.'R') then
                  do nmod=0,maxover! loop over all branches
                     do iker=1,nker
                        do ipar=1,npar
                           do imod=1,nmods(nmod)
                              array(imod)=parms(imod,iker,ipar,nmod)
                              omegain(imod)=omegs(imod,nmod) ! BUG FIXED LB 13.5.09
                           enddo
                           call rspln(1,nmods(nmod),omegain,array,arraycube,work)
                           do iper=1,nperiods(1,nmod)
                              omega=twopi/true_per(iper,1,nmod)
                              value=rsple(1,nmods,omegain,array,arraycube,omega)
                              pargrid(ilat,ilon,iker,ipar,iper,nmod)=value
                           enddo
                        enddo
                     enddo
                     do imod=1,nmods(nmod)
                        array(imod)=pvels(imod,nmod)
                        omegain(imod)=omegs(imod,nmod)
                     enddo
                     call rspln(1,nmods(nmod),omegain,array,arraycube,work) !bug fixed lb 13.5.09
                     do iper=1,nperiods(1,nmod)
                        omega=twopi/true_per(iper,1,nmod)
                        value=rsple(1,nmods(nmod),omegain,array,arraycube,omega) !bug fixed lb 13.5.09
                        valueref=c0(iper,1,nmod)
                        ! pvelgrid: absolute difference in slowness between 3d starting model and prem
                        pvelgrid(ilat,ilon,iper,nmod)=-(value-valueref)/(valueref*valueref)
                        ! correct pargrid: we want slowness not velocity kernels:
                        do iker=1,nker
                           do ipar=1,npar
                              pargrid(ilat,ilon,iker,ipar,iper,nmod)=-pargrid(ilat,ilon,iker,ipar,iper,nmod)/(value*value)
                           enddo
                        enddo
                     enddo
                  enddo !enf of loop over branches
               endif !end of if L/R
            endif !if ios.eq.0
         enddo !latitude
      enddo !longitude
      return
      end subroutine load_kernels! end of subroutine load_kernels


c
c =====================================================================
c  Subroutine delazs
c
      subroutine delazs(eplat,eplong,stlat,stlong,delta,azep,azst)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      data hpi,twopi,rad,reprad/1.5707963268
     1   ,6.2831853,.017453293,57.2957795/
      arcos(x)=atan2(sqrt(1.-x*x),x)
      el=eplat*rad
      el=hpi-el
      stl=stlat*rad
      stl=hpi-stl
      elon=eplong*rad
      slon=stlong*rad
      as=cos(stl)
      bs=sin(stl)
      cs=cos(slon)
      ds=sin(slon)
      a=cos(el)
      b=sin(el)
      c=cos(elon)
      d=sin(elon)
      co=1.0
      cdel=a*as+b*bs*(c*cs+d*ds)
      if(abs(cdel).gt.1.0) cdel=sign(co,cdel)
      delt=arcos(cdel)
      delta=delt*reprad !distance
      sdel=sin(delt)
      caze=(as-a*cdel)/(sdel*b)
      if(abs(caze).gt.1.0) caze=sign(co,caze)
      aze=arcos(caze)
      if(bs.gt.0.0) cazs=(a-as*cdel)/(bs*sdel)
      if(bs.eq.0.0) cazs=sign(co,cazs)
      if(abs(cazs).gt.1.0) cazs=sign(co,cazs)
      azs=arcos(cazs)
      dif=ds*c-cs*d
      if(dif.lt.0.0) aze=twopi-aze
      azep=reprad*aze
      if(dif.gt.0.0) azs=twopi-azs
      azst=reprad*azs
      return
      end subroutine delazs





c
c =====================================================================
c  Parameterization
c

        subroutine param(eq_incr,nsqrs,nsqtot,nlatzones,numto,iswit,
     &                   compgrid,nlatzomax,iadapt)
c
c *** this version is parameterized exactly as julias version
c *** and allows to distinguish between adaptive and regular
c *** parameterization, L.A. Aug. 2012
c
        integer nlatzones0
        integer nsqrs(nlatzomax),nsqtot(nlatzomax+1)
        integer nsqrs0(nlatzomax),nsqtot0(nlatzomax+1)
	parameter(pi=3.1415926536)

c
c *** this is ludwigs preferred regular parameterization
c
        if(iadapt.eq.0)then
           print*,"Using Ludwigs prefered regular grid parameterization"
           numto=0
           colat=-eq_incr/2.
           do k=1,nlatzones
              colat=colat+eq_incr
              theta=(colat/180.)*pi
              ! for this latitudinal zone, compute number of blocks (nsqrs)
              deltalon=eq_incr/(sin(theta))
              nsqrs(k)=(360./deltalon)+1
              ! if requested, correct nsqrs(k) so the grid is compatible to reference grid
              if(iswit.eq.1)then
                 if(360./nsqrs(k).ge.compgrid)then
 100             if((mod(360./nsqrs(k),compgrid).ne.0).or.
     &              (mod(nsqrs(k),2).ne.0))then
                    nsqrs(k)=nsqrs(k)+1
                    goto 100
                 else
                 endif
              elseif(360./nsqrs(k).lt.compgrid)then
 101             if((mod(compgrid,360./nsqrs(k)).ne.0).or.
     &              (mod(nsqrs(k),2).ne.0))then
                    nsqrs(k)=nsqrs(k)-1
                    goto 101
                 else
                 endif
              endif
           else ! modified by lud
              if(mod(nsqrs(k),2).ne.0)nsqrs(k)=nsqrs(k)-1
	   endif
         if(MOD(NSQRS(K),2).ne.0)then
              stop "nsqrs has to be even"
         endif
         ! print regular grid parameterization
         print*,k,nsqrs(k)
         nsqtot(k)=numto
         numto=numto+nsqrs(k)
         enddo
	 nsqtot(nlatzones+1)=numto         
c
c *** this is julias adaptive grid, which is only slighly differs from
c *** my version of param. In fact, there is not really a reason why 
c *** i distinguish b/w the both... historical reasons
c
         elseif(iadapt.eq.1)then
           print*,"Using julias prefered adaptive grid parameterization"
            n1layer0=0
            nlatzones0=180./compgrid
            colat=-compgrid/2.

            do k=1,nlatzones0
               colat=colat+compgrid
               theta=(colat/180.)*pi
               ! for this latitudinal zone, compute number of blocks (nsqrs)
               deltalon=compgrid/(sin(theta))
               nsqrs0(k)=(360./deltalon)+1
               if(mod(nsqrs0(k),2).ne.0)nsqrs0(k)=nsqrs0(k)-1
               ! if requested, correct nsqrs(k) so the grid is compatible to reference grid
               if(iswit.eq.1)then
                  if(360./nsqrs0(k).ge.compgrid)then
 102              if(mod(360./nsqrs0(k),compgrid).ne.0)then
                     nsqrs0(k)=nsqrs0(k)+1
                     goto 102
                  else
                  endif
               elseif(360./nsqrs0(k).lt.compgrid)then
 103              if(mod(compgrid,360./nsqrs0(k)).ne.0)then
                     nsqrs0(k)=nsqrs0(k)-1
                     goto 103
                  else
                  endif
               endif
	       endif
               
               if(mod(nsqrs0(k),2).ne.0)stop "nsqrs has to be even"
               ! print rough parameterization acc. julia
               print*,k,nsqrs0(k)
               nsqtot0(k)=n1layer0
               n1layer0=n1layer0+nsqrs0(k)
            enddo

            ! define fine grid
            fact=compgrid/eq_incr
            n1layer=0
            print*,"ratio of coarser to finest grid is",fact
            do k=1,nlatzones
               k0=((k-1)/fact)+1
               nsqrs(k)=nsqrs0(k0)*fact
               nsqtot(k)=n1layer
               n1layer=n1layer+nsqrs(k)
               print*, "k=",k,"k0=",k0," nsqrs=", nsqrs(k),nsqrs0(k0)," lon=", 360./nsqrs(k)
               print*, "latitudinal zone =",k, "lon=", 360./nsqrs(k), "#pixels=", nsqrs(k) 
            enddo
            nsqtot(nlatzones+1)=n1layer
            print*,'Number of pixel with finest parameterization: n1layer=',n1layer
            numto=n1layer            
         endif
       
         return
         end subroutine param
