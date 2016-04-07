!====================================================================
!
! Least-squares inversion, original: 
! originally from L. Boschi, 2/2009; refactored and ported to 
! Fortran 90 port by L. Auer in 2012; variable grid roughness damping
!
!====================================================================

	program flexinv90

        use flexinv_module ! contains all global variables
        implicit none

!======= LOCAL VARIABLES ============================================        

        ! LSQR arrays
        real, allocatable, dimension(:) :: se !nmax
        real, allocatable, dimension(:) :: aty
	real, allocatable, dimension(:) :: v
        real, allocatable, dimension(:) :: w
        real, allocatable, dimension(:) :: x
	real, dimension(300)  :: npnt

        ! LSQR variables
        real*4      :: damp0,atol0,btol0
        real*8      :: anorm, acond, rnorm
        real*8      :: arnorm, xnorm
        integer*8   :: ilim0,nout0,clim0,itn

        ! Variance reduction
        real*8      :: rnum,dnum,vtot,varred
        real*4      :: xsolrms,subvr,iw,tinit
        real*4      :: roughness,rowprod,inval
        real, allocatable, dimension(:) :: wei

        ! Various variables
        integer*8   :: n1layer,ndata,nperlay
        integer*8   :: nplt,nsets,ivre
        integer*8   :: dmy1,dmy2,dmy3
        integer*8   :: istop,info
        
        ! Local running indices
        integer*8   :: jj,kk,ll,zz
        integer*8   :: pp,ii,ii0    ! most important local indices
        integer*8   :: i,l,k,p,u,j  ! running indices       
        integer*8   :: inp,lev 

        ! hardcoded LSQR parameters
        damp0   = 0.d0   ! additional norm damping factor, not needed
        atol0   = 0.0001 ! 0.000001
        btol0   = 0.0001 ! 0.000001 
        clim0   = 0.d0   ! con limit
        ilim0   = 10000  ! iteration limit
        nout0   = 1      ! dont know
        iata    = 0      ! lsqr on the rectangular system or cholesky on normal equations

!
! NOTE: THE ATA/CHOLESKY MODE DOESN'T WORK FOR SOME REASONS 
!       PRESUMABLY DESTROYED IT WHEN CONVERTING TO F90. 
!       DEACTIVATING IT FOR THE MOMENT
!

        namexxx = ""       

!======= READ FROM STDIN ============================================        
        print*,""
        print*,""
        print*,""
        print*,""
        print*," _|_|_|_|  _|        _|_|_|_|  _|      _|  _|_|_|  _|      _|  _|      _|  "
        print*," _|        _|        _|          _|  _|      _|    _|_|    _|  _|      _|  "
        print*," _|_|_|    _|        _|_|_|        _|        _|    _|  _|  _|  _|      _|  "
        print*," _|        _|        _|          _|  _|      _|    _|    _|_|    _|  _|    "
        print*," _|        _|_|_|_|  _|_|_|_|  _|      _|  _|_|_|  _|      _|      _|      "
        print*,""
        print*,""
        print*,""
        print*,"Welcome to Flexinv version alpha 0.1 (c) L.Auer, L.Boschi, J.Schaefer 2013"
        
        print*,""
        print*,""
        print*,"--------------------- READ PARAMETERS FROM STDIN ------------------------"
        print*,""
        print*,""
        
        write(*,'(A,$)'),"Enter dimension of right-hand-side vector?"
        read*,m
        write(*,*),m
        write(*,'(A,$)'),"Enter dimension of sparse matrix arrays?"
        read*,nonz
        write(*,*),nonz

!     -------------------------------------------------------------
!       read inversion parameters from stdin

	write(*,'(A,$)'),"Enter number of inversion parameters (1-4):"
	read*,npar
        write(*,'(3f8.3)'),npar

	write(*,'(A,$)'),"Enter equatorial cell size in degree:"
	read*,eqincr
        write(*,'(3f8.3)'),eqincr

        write(*,'(A,$)'),"Compatible with crustal model of same gridsize (0=no,1=yes): "
	read*,iswit
        write(*,*),iswit

        write(*,'(A,$)'),"Enter number of layers: "
	read*,nlay
        write(*,*),nlay

        write(*,'(A)'),"Working directory (in which results are stored):"
        read*,wdir
        write(*,*),trim(wdir)

        write(*,'(A,$)'),"Adaptive grid (1) or non (0): "
        read*,iadapt
        write(*,*),iadapt

        if (iadapt.ge.1) then
           write(*,'(A)'),"File containing adaptive grid info:"
           read*,gridinfo
           write(*,*),trim(gridinfo)
        end if

!======= DEFINE/READ PARAMETERIZATION ===============================         

        ! some basics variables
        nlatzones = 180 / eqincr

        ! regular grid parameterization 
        if (iadapt.eq.0) then ! grid is regular
           refgrid = eqincr  ! refgrid in the not adaptive case!
           write(*,'(A)'),"Calling subroutine param to define regular grid:"
           call param(eqincr,nsqrs,nsqtot,&
                      nlatzones,n1layer,iswit,&
                      refgrid,nlatzomax,iadapt)
           n0 = n1layer ! number of pars per layer
           n  = n0 * nlay * npar ! number of inversion parameters
           n1 = n0 * nlay ! number of pars per phys. par
           write(*,'(A)'),"... successful!"

        ! adaptive grid parameterization
        else if (iadapt.ge.1) then ! grid is adaptive
           write(*,'(A,$)'),"Reading adaptive grid information from file: "
           refgrid = 5.0 ! reference grid is always 5.0 degree
           n  = 0
           n0 = 0
           open(620,file=trim(gridinfo)//".lay",status="old")
           do k = 1,nlay
              read(620,*), adpx1layer(k)
              if (k.eq.1) then 
                 n0 = adpx1layer(k)
              endif
              n = n + adpx1layer(k)
           end do
           n1 = n
           n  = n * npar
           write(*,*),"... successful!"
        end if        

        ! initialize some more things
        ii       = 1 ! set rhs index
	jj       = 1 ! set elm index
	nsets    = 1 ! number of datasets
	tinit    = secnds(0.0)

!
! NOTE: ATA/CHOLESKY MODE DEFUNCTIONAL
!       natamax  = n*(n+1)   ! calculate number of val in upper triangular of ata
!       natamax  = natamax/2                            
!

        nperlay  = n/(2*nlay) ! number of columns/pixels per layer
        nplt     = nperlay*nlay  ! here sensitivity of 2nd parameter starts

!     -------------------------------------------------------------
!       Check if something is wrong

	if (nlay.gt.nlaym) then
	   print*,"nlay = ",nlay," is larger than nlaym = ",nlaym
	   stop "ERROR: too many layers, stopping ..."
	end if
	if (n0.gt.n0max) then 
           print*,"n0 = ",n0," is larger than n0max = ",n0max
           stop "ERROR: too many pixels, stopping ..."
        end if
	if (n.gt.nmax) then
           print*,"n = ",n," is larger than nmax = ",nmax
           stop "ERROR: too many voxels, stopping ..."
        end if

	write(*,'(A,$)'),'Total number of parameters: '
        write(*,*),n
        print*,""


!======= ALLOCATE MEMORY ============================================        


        allocate(se(nmax),aty(nmax),v(nmax),w(nmax),x(nmax))        

        if (iata.eq.0) then
           allocate(val(nonz))   ! large m x n val array
           allocate(ind(nonz))   ! large m x n ind array
           allocate(val0(ntemp)) ! temporary val array
           allocate(ind0(ntemp)) ! temporary ind array
           allocate(rhs(m))      ! rhs data vector
           allocate(rhs0b(m))    ! temp rhs vector
           allocate(pnt(0:m))    ! point verctor
           allocate(wei(0:m))    ! weight verctor

!
! NOTE: ATA/CHOLESKY MODE DEFUNCTIONAL
!
!        elseif(iata.eq.1) then
!           allocate(val(ntemp))  ! just temp array
!           allocate(ind(ntemp))  ! just temp array
!           allocate(val0(ntemp))
!           allocate(ind0(ntemp))
!           allocate(rhs0b(m))
!           allocate(pnt(0:m))
!           allocate(wei(0:m))    
!           allocate(rhs0(m))
!           allocate(ata(natamax)) ! ata matrix
!           allocate(atd(n))       ! atd vector
!           allocate(rhs(m))
!           allocate(pnt0(0:m))
!           do k=1,natamax
!              ata(k)=0.
!              if (mod(k,int(natamax/100)).eq.0) then
!                 print*,"Initialising ata,",int(100*k/natamax)+1,"%"
!              end if
!           end do
!           do k=1,n
!              atd(k)=0.
!              if (mod(k,int(n/100)).eq.0) then
!                 print*,"Initialising atd,",int(100*k/n)+1,"%"
!              end if
!           end do
        else
           stop "ERROR: iata chosen incorrectly, stopping ..."
        endif




!======= READ SUBMATRICES ===========================================

        print*,""
        print*,""
        print*,"-------------------------- READ SUBMATRICES -----------------------------"
        print*,""
        print*,""


        ! MAJOR LOOP OVER SUBMATRICES
        pnt(0) = 0 ! necessary for readmatrix
        do while (1.eq.1)

           npnt(nsets)=ii
 
           print*,"Path to VAL array of sub-matrix:"
           read(*,*)namexxx
           print*,trim(namexxx)

           if(namexxx.eq."finished") then
              goto 222
           end if
           
           print*,"Path to IND array of sub-matrix:"
           read*,nameind
           print*,trim(nameind)
           print*,"Path to PNT vector of sub-matrix:"
           read*,namepoi
           print*,trim(namepoi)
           print*,"Path to RHS vector of sub-matrix:"
           read*,namerhs
           print*,trim(namerhs)          
           write(*,'(A,$)'),"Relative weight of submatrix:"
           read*,relwei
           write(*,*),relwei
           write(*,'(A,$)'),"Threshold delay times for downweighting (pos,neg):"
           read*,cutoff,cutoffn
           write(*,*),cutoff,cutoffn
           write(*,'(A,$)'),"Number of observations in this subset:"
           read*,ndata
           write(*,*),ndata

                    
           if (iata.eq.0) then

              open(33,file=namepoi,status='old')
              open(77,file=namerhs,status='old')
              print*,"ii,jj,n: ",ii,jj,n

              ! read pointer and rhs
              ii0=ii

              do ii=ii0,ii0+ndata-1  ! loop over pointer val
                 read(33,*) zz ! read pointer value up to which val of first measurement of subset is going
                 pnt(ii) = zz + pnt(ii0-1) ! this is the pointer value up to which this subset is going
                 read(77,*) rhs(ii) ! read rhs value for the first measurement in this subset

                 ! assign additional weight to penalize large anomalies (outliers)
                 adtime  = abs(rhs(ii))
                 wei(ii) = 1.

                 if(adtime.gt.cutoff)then
                    wei(ii) = exp(cutoff-adtime)
                 end if
                 if(adtime.lt.cutoffn)then
                    wei(ii) = exp(-cutoffn+adtime)
                 end if

                 rhs(ii) = rhs(ii) * relwei * wei(ii)

              enddo
              close(33)
              close(77)

              open(1,file=namexxx,status='old',access='direct',form='unformatted',recl=4*zz)
              open(4,file=nameind,status='old',access='direct',form='unformatted',recl=4*zz)
              read(4,rec=1) (ind0(kk), kk=1,zz) ! index in row = index of corresponding pixel
              read(1,rec=1) (val0(kk), kk=1,zz) ! val: raypath through this pixel
              print*,"... sub-matrices read, applying weight of ",relwei

              ! apply weighting and store weighted cummulative sensitivitie             
              pp = ii0
              do kk=jj,jj+zz-1
                 if (kk.gt.pnt(pp)) then
                    pp = pp + 1                    
                 end if
                 ind(kk)=ind0(kk-jj+1)
                 ind0(kk-jj+1)=0
                 if(ind(kk).gt.n)then
                    print *,"ERROR: undefined voxel index", n, ind(kk)
                    stop
                 endif
                 val(kk)=val0(kk-jj+1) * relwei * wei(pp)
                 val0(kk-jj+1)=0.
              end do
              jj = jj + zz
              close(4)
              close(1)

! NOTE: I HAVE DEACTIVATE THE ATA MODE
!
!           elseif (iata.eq.1) then
!
!              open(33,file=namepoi,status='old')
!              open(77,file=namerhs,status='old')
!              print*,"ii,jj,n: ",ii,jj,n
!              ! read pointer and rhs
!              ii0=ii
!              pnt0(0)=0.
!              do ii=1,ndata ! loop over pointer val
!                 read(33,*)zz 
!                 pnt0(ii)=zz
!                 read(77,*)rhs0(ii) ! read rhs value for the first measurement in this subset
!                 ! assign additional weight to penalize large anomalies (outliers)
!
!                 ! I don't exclude outliers in the cholesky mode
!                 ! adtime=abs(rhs0(ii))
!                 ! wei2=1.
!                 ! if(adtime.gt.cutoff)then
!                 !    wei2=exp(cutoff-adtime)
!                 ! endif
!                 ! if(adtime.lt.cutoffn)then
!                 !    wei2=exp(-cutoffn+adtime)
!                 ! endif
!
!                 rhs0(ii)=rhs0(ii)*relwei !*wei2
!              enddo
!              close(33)
!              close(77)
!              
!              open(1,file=namexxx,status='old',access='direct',form='unformatted',recl=4*zz)
!              open(4,file=nameind,status='old',access='direct',form='unformatted',recl=4*zz)
!              read(4,rec=1) (ind0(kk), kk=1,zz) ! index in row = index of corresponding pixel
!              read(1,rec=1) (val0(kk), kk=1,zz) ! val: raypath through this pixel
!              print*,"... sub-matrices read, amending ATA ",relwei
!
!              do ii=1,ndata
!                 irec=0
!
!                 if(mod(ii,int(ndata/100)).eq.0) then
!                    print*,"contribution_ata,",int(100*ii/ndata)+1,"%"
!                 end if
!
!
!                 ! I don't exclude outliers, Can be ignored
!                 ! adtime=abs(rhs(ii))
!                 ! wei2=1.
!                 ! if(adtime.gt.cutoff)then
!                 !    wei2=exp(cutoff-adtime)
!                 ! endif
!                 ! if(adtime.lt.cutoffn)then
!                 !    wei2=exp(-cutoffn+adtime)
!                 ! endif
!
!                 do kk=pnt0(ii-1)+1,pnt0(ii)                                  
!                    irec=irec+1
!                    val(irec)=val0(kk)*relwei !*wei2
!                    val0(kk)=0.
!                    ind(irec)=ind0(kk)
!                    if(ind(irec).eq.0) then
!                       print*,"ERROR"
!                       stop
!                    end if
!                    if(ind(irec).gt.n)then
!                       print *,"ERROR: undefined voxel index", n, ind(irec)
!                       stop
!                    endif
!                    ind0(kk)=0
!                 end do
!
!                 call contribution_ata2(val(1:irec),ind(1:irec),ata,rhs0(ii),&
!                                       atd,n,natamax,irec)
!              end do
!              close(4)
!              close(1)
!
           else
              stop "ERROR: iata chosen incorrectly, stopping ..."
           endif
            
           nsets = nsets + 1 ! # of datasets

           
           print*,""
           print*,"* * * * * * * * Statistic * * * * * * * "
           print*,'Data so far: ',ii,"out of max",m
           print*,'Nonz so far: ',jj-1,"out ot max",nonz
           print*,"~ entries per datum: ",float(jj-1)/float(ii)
           print*,"* * * * * * * * * * * * * * * * * * * * "
           print*,""


	end do ! MAJOR LOOP OVER SUBMATRICES
  
        ! initialize, want to start at 0 for damping in the iata=1 case
        ii=1
        jj=1
   
222     continue

        print*,"No more matrices to read!"
        print*,"Matrices read in",secnds(tinit)," seconds!"



!======= DO A COUPLE OF THINGS ======================================
        
	nrhs      = ii - 1 ! global nrhs = local ii
	nelm      = jj - 1 ! global nelm = local jj      
	nrhs0     = nrhs
        nelm0     = nelm 
	do u = 1,nrhs0
	   rhs0b(u) = rhs(u)
	enddo
        do u = nrhs0+1,nrhs 
           rhs(u) = 0.d0
        enddo
        do u = nelm0+1,nelm
           val(u) = 0.d0
        enddo
        nelm      = nelm0
        nrhs      = nrhs0
	do u = 1,nrhs0
	   rhs(u) = rhs0b(u) ! t gets modified by lsqr
	enddo        


!======= READ DAMPING SCHEME ======================================

        write(*,'(A,$)'),"Roughness damping scheme: "
        read*,dmpsched_r
        write(*,*),trim(dmpsched_r)
        write(*,'(A,$)'),"Difference damping scheme: "
        read*,dmpsched_d
        write(*,*),trim(dmpsched_d)
        write(*,'(A,$)'),"Additional weight factor, Roughness damping: "
	read*,wrdamp
        write(*,*),wrdamp
        write(*,'(A,$)'),"Additional weight factor, Similarity damping: "
	read*,wddamp
        write(*,*),wddamp
        write(*,'(A,$)'),"Additional weight factor, Norm damping: "
	read*,wndamp
        write(*,*),wndamp

        write(*,'(A,$)'),"Prescribe anisotropy, yes (1) or no (0)?"
        read*,preani
        write(*,*),preani
        if (preani.gt.0) then
           read*,preani_model
           write(*,*),preani_model
        elseif (preani.eq.-1.) then
           read*,preani_model
           write(*,*),preani_model
        end if

        write(*,'(A,$)'),"Name of solution model file: "
        read*,outfile
        write(*,*),trim(outfile)

        if (wndamp.gt.0) then
           if (iadapt.eq.1) then
              iadaptnormdmp=1
              inormdmp=0
           elseif (iadapt.eq.0) then
              iadaptnormdmp=0
              inormdmp=1
           end if
        end if
        
        print*,""
        print*,""
        print*,"------------------ AMENDING A MATRIX WITH DAMPING EQUATIONS ------------- "
        print*,""
        print*,""


!======= APPLY ROUGHNESS DAMPING ====================================

        ! We always need to apply roughness damping, so...
        print*,'Applying variable/regular grid roughness damping'
        print*,"Nelm and nrhs before damping:", nelm,nrhs

        ! Read damping schedule
        open(7878,file=trim(dmpsched_r),status="old")
        read(7878,*)
        read(7878,*)
        read(7878,*)
        read(7878,*)
        read(7878,*),af1,af2,af3,af4,afv
        read(7878,*)
        do l=1,nlay
           read(7878,*),wgradh(l,1),wgradh(l,2),wgradh(l,3),wgradh(l,4),wgradv(l)
           wgradh(l,1) = wgradh(l,1) * wrdamp * af1 ! multiply by global damp factor
           wgradh(l,2) = wgradh(l,2) * wrdamp * af2 ! multiply by global damp factor
           wgradh(l,3) = wgradh(l,3) * wrdamp * af3 ! multiply by global damp factor
           wgradh(l,4) = wgradh(l,4) * wrdamp * af4 ! multiply by global damp factor
           wgradv(l)   = wgradv(l) * afv
           print*,"Read rdamp facs"
           print*,l,wgradh(l,1),wgradh(l,2),wgradh(l,3),wgradh(l,4),wgradv(l)
        end do
        close(7878)
       
        ! call subroutine
        call damp_roughness
                                         
        ! Some statistics 
        print*," "
        print*,"* * * * * * * * Statistic * * * * * * * "
        print*,'Data so far: ',nrhs,"out of max",m
        print*,'Nonz so far: ',nelm,"out ot max",nonz
        print*,"* * * * * * * * * * * * * * * * * * * * "
        print*," "


!======= APPLY NORM DAMPING =========================================
        ! the regular grid case
	if (inormdmp.eq.1) then
	   print*,'Applying additional norm damping in the mid-mantle, regular grid'
           ! weights for add ndamp
           ndampvec0(1)=(wndamp**2)*((eqincr*111.11)**2)
	   ndampvec0(2)=(wndamp**2)*((eqincr*111.11)**2)
 	   do p=1,npar
	      dampsize=ndampvec0(p)
	      print*,'damp parameter',p,';',dampsize
	      inp=n1*(p-1)
	      call damp_norm(dampsize,nelm,nrhs,val,ind,pnt,rhs,m,inp,nonz,n1)
	   end do
           print*,""
           print*,"* * * * * * * * Statistic * * * * * * * "
           print*,'Data so far: ',nrhs,"out of max",m
           print*,'Nonz so far: ',nelm,"out ot max",nonz
           print*,"* * * * * * * * * * * * * * * * * * * * "
           print*,""
	end if
        ! the same for the adaptive case
        if (iadaptnormdmp.eq.1) then
           print*,'Applying norm damping, variable grid'
           area(1)=(5.0*111.1)
           area(2)=(2.5*111.1)
           area(3)=(1.25*111.1)
           area(4)=(0.625*111.1)              
           open(610,file=trim(gridinfo),status="old")
           adpxtotal=0
           do k=1,nlay
              adpxtotal=adpxtotal+adpx1layer(k)
              print*,'Total number of adaptive parameters',adpxtotal
           end do
           do j=1,adpxtotal
              read(610,*) dmy1, lev, dmy2, dmy3
              ndampvec1(j)=wndamp*area(lev)
           end do
           close(610)
           do p=1,npar
	      inp=n1*(p-1)
              call damp_norm_adapt(ndampvec1,nelm,nrhs,val,ind,&
                                   pnt,rhs,m,inp,nonz,n1)
           end do
           print*,""
           print*,"* * * * * * * * Statistic * * * * * * * "
           print*,'Data so far: ',nrhs,"out of max",m
           print*,'Nonz so far: ',nelm,"out ot max",nonz
           print*,"* * * * * * * * * * * * * * * * * * * * "
           print*,""
        end if


!======= APPLY ANISO DAMPING ========================================
        if (wddamp.gt.0.) then

           ! Read difference damping schedule
           open(7879,file=trim(dmpsched_d),status="old")
           read(7879,*)
           read(7879,*)
           read(7879,*)
           read(7879,*)
           read(7879,*),ad1,ad2
           read(7879,*)
           do l=1,nlay
              read(7879,*),wddampvec(l,1),wddampvec(l,2)
              wddampvec(l,1) = wddampvec(l,1) * wddamp * ad1 ! multiply by global damp factor
              wddampvec(l,2) = wddampvec(l,2) * wddamp * ad2 ! multiply by global damp factor
              print*,"Read ddamp facs"
              print*,l,wddampvec(l,1),wddampvec(l,2)
           end do
           close(7879)

           if (npar.gt.1) then
              print*,'Damp anisotropy with weight',wddamp

              call damp_difference 

              print*,""
              print*,"* * * * * * * * Statistic * * * * * * * "
              print*,'Data so far: ',nrhs,"out of max",m
              print*,'Nonz so far: ',nelm,"out ot max",nonz
              print*,"* * * * * * * * * * * * * * * * * * * * "
              print*,""
           else
              print*,"Difference damping between only makes sense"
              print*,"in the case of 2,3 or 4 physical parameters"
           end if
        end if


!======= PRESCRIBE ANISOTROPY ======================================
        if (preani.gt.0.) then
           print*,'Prescribe anisotropy model with weight',preani
           call prescribe_aniso(preani,preani_model,nelm,nrhs,val,ind,pnt,rhs,m,n1,nonz)
           print*,""
           print*,"* * * * * * * * Statistic * * * * * * * "
           print*,'Data so far: ',nrhs,"out of max",m
           print*,'Nonz so far: ',nelm,"out ot max",nonz
           print*,"* * * * * * * * * * * * * * * * * * * * "
           print*,""
        end if     
        if (preani.eq.-1.) then
           ! Skip inversion and just test variance reduction for 
           ! a specified input model (of vs and ani)
           open(747,file=trim(preani_model))
           print*,"Testing variance reduction"
           do k=1,n1*2
              print*,"Parameter: ",k
              read(747,*) dmy1,inval
              x(k)=inval
           enddo
           close(747)
           goto 555 ! jump to variance red comp
        end if     



!======= SOLVE MATRIX ===============================================
        print*,""
        print*,""
        print*,"------------------------- SOLVING LINEAR SYSTEM ------------------------- "
        print*,""
        print*,""

        tinit = secnds(0.0)
        if (iata.eq.0) then

           print*,""
           print*,""
           print*,"LSQR on rectangular system of size n x m"
           print*,""
           print*,""

	   call lsqr( m, n, damp0, 1, n, iw, aty, rhs, v, w, x, se,&
                      atol0, btol0, clim0, ilim0, nout0, istop, itn,& 
                      anorm, acond, rnorm, arnorm, xnorm, ind,&
                      val, pnt, nonz, nrhs)

           print*,""
           print*,""
           print*,"Matrix solved! Lsqr inversion ran in",secnds(tinit)," seconds"
           print*,""
           print*,""

!
! NOTE: ATA/CHOLESKY MODE IS DEACTIVATED FOR THE MOMENT
!
!        elseif (iata.eq.1) then
!
!           print*,""
!           print*,""
!           print*,"CHOLESKY on normal equations"
!           print*,""
!           print*,""
!
!           call sppsv( 'U', n, 1, ata, atd, n, info )
!
!           do k=1,n
!              x(k)=atd(k)
!           end do
!
!           print*,""
!           print*,""
!           print*,"Matrix solved! Cholesky inversion ran in",secnds(tinit)," seconds"
!           print*,""
!           print*,""
!           
        else

           stop "ERROR: iata chosen incorrectly, stopping ..."           

        end if


!     -------------------------------------------------------------
!       compute cummulative variance reduction
555     continue
	rnum=0.
	dnum=0.
	do l=1,nrhs0
	   dnum=dnum+(rhs0b(l)*rhs0b(l))
	   vtot=0.
	   do p=pnt(l-1)+1,pnt(l)
	      vtot=vtot+(val(p)*x(ind(p)) )
	   enddo
	   rnum=rnum+(rhs0b(l)-vtot)**2
	enddo
	varred=1.-(rnum/dnum)
	print*,"Cumulative variance reduction",varred

        ! store cumulative variance reduction and log-likelihood
 	open(72,file="fits."//trim(outfile)//".dat")
	write(72,*)wrdamp,wgradv(1),wndamp,itn,varred
	close(72)

        ! log-likelihood function
	open(72,file="logl."//trim(outfile)//".dat")
	write(72,*)wrdamp,wgradv(1),wndamp,-0.5*float(nrhs0)*log(rnum),itn,n
	close(72)


        ! variance reduction for each dataset used in this inversion
	print*,'Computing vr for',nsets-1,' subsets of data'
	open(99,file="varr."//trim(outfile)//".dat")
	do ivre=1,nsets-1
	   rnum=0.
	   dnum=0.
	   xsolrms=0.
	   do l=npnt(ivre),npnt(ivre+1)-1
	      dnum=dnum+(rhs0b(l)*rhs0b(l))
	      vtot=0.
	      do p=pnt(l-1)+1,pnt(l)
	         vtot=vtot+(val(p)*x(ind(p)))
	      enddo
	      xsolrms = xsolrms + vtot**2
	      rnum=rnum+(rhs0b(l)-vtot)**2
	   enddo
	   subvr=1.-(rnum/dnum)
	   xsolrms = sqrt(xsolrms/(npnt(ivre+1)-npnt(ivre)-1))! rms of solution
	   print*,'set',ivre,' vr ',subvr,' rms ',xsolrms,' np ',npnt(ivre+1)-npnt(ivre)-1
	   write(99,*)subvr
	enddo


!     -------------------------------------------------------------
!       Compute normalized model roughness

        print*,"Computing normalized model roughness, need to recompute damping!"
	xnorm=0.
        do i=1,n
	   xnorm=xnorm+x(i)*(x(i)*100)
        enddo
        nrhs=0
        nelm=0
        do ll=1,nlay
           wgradh(ll,1) = 1.0
           wgradh(ll,2) = 1.0
           wgradv(ll)   = 1.0
        end do

        call damp_roughness

        ! dot product with model vector
        roughness=0.
        do i=1,nrhs
           rowprod=0. !roughness=0.,rowprod
           do j=pnt(i-1)+1,pnt(i)
              rowprod=rowprod+(val(j)*x(ind(j)))
           end do
           roughness=roughness+(rowprod*rowprod)
        end do
        open(35,file="rough."//trim(outfile)//".dat")
        write(35,*),roughness,roughness/xnorm
        close(35)
        print*,"... done"


!     -------------------------------------------------------------
!       Write solution to disc

        print*,"Writing solution to: ",outfile
        open(40,file="sol."//trim(outfile)//".dat")
        do k=1,n
           write(40,*)k,x(k)
        enddo
        close(40)


        ! -------------------------------------------> WE ARE DONE!
        print*,""
        print*,""
        print*," ____ ____ ____ ____ "
        print*,"||D |||o |||n |||e ||"
        print*,"||__|||__|||__|||__||"
        print*,"|/__\|/__\|/__\|/__\|"
        print*,""
        print*,""

        stop "... check your model!" ! end everything when encountering "666"


	end program flexinv90
!=======================================================================

