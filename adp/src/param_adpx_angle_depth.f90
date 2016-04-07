!--------------------------------------------------------------
! from L. Auer 2012 and J. Schaefer 2009
! Major rewrite by Ludwig to allow for depth dependent
! variable (adaptive) grid parameterization
!


subroutine param_adpx

  use project_adpx_module
  implicit none

  real                  :: xlamin,xlamax,xlomin,xlomax
  real                  :: colat,theta,deltalon,size               
  real                  :: xx1,xx2
  
  character*3           :: i3   
  character*1           :: dmy3

  integer	        :: k,k0,ixx,jxx,layers
  integer               :: lev0
  integer               :: nlatzones0                 	! latitudinal zones in reference grid
  integer		:: ipxnew, ipxnewtmp			! total number of new pixel
  integer               :: ifa,ipar                  	! actual pixel size / smallest pixel size
  integer               :: ipx,ipx0,ipxin0    
  integer               :: layer,layer0                	! index for loop over layers
  integer               :: ipx1layer,tmp,summe

  integer,dimension(4**nlev)  :: iwithin,iwithin_tmp			! indexes of the smallest pixel contained in a pixel
  integer,dimension(4)	      :: iwithin_next
  integer,dimension(4)	      :: htctmp,htctmp_layernumber	! temporary hitcount for all smallest pixel in one pixel in first layer and layer for which parameterization is defined
  integer,dimension(nlatzomax):: nsqrstmp,nsqtottmp
  integer,dimension(4)	      :: wedgecnt

  integer,dimension(:,:),allocatable	:: htctmp_each
  integer,dimension(:,:),allocatable	:: htctmp_each_layernumber,angle_hitcount

  
  allocate(htctmp_each(4,4**(nlev-1)))
  allocate(htctmp_each_layernumber(4,4**(nlev-1)))
  allocate(angle_hitcount(4,n1layer(nlev)))
  
  nlatzones0=180/refgrid 
  print*, "Number of latitudinal zones, nlatzones = ", nlatzones0

! -------------------------------------------------------------
! define rough grid, compatible with reference grid
 
  lev=1
  colat=-refgrid/2.
  n1layer(lev)=0
  
  do k=1,nlatzones0
      colat=colat+refgrid
      theta=(colat/180.)*pi

      ! for this latitudinal zone, compute number of blocks (nsqrs)
      deltalon=refgrid/(sin(theta))
      nsqrs(k,lev)=(360./deltalon)+1
      if(mod(nsqrs(k,lev),2).ne.0)nsqrs(k,lev)=nsqrs(k,lev)-1

      ! if requested, correct nsqrs(k,lev) so the grid is compatible to reference grid
      if(iswit.eq.1)then
           if(360./nsqrs(k,lev).ge.refgrid)then                   ! blocks > reference grid
100             if(mod(360./nsqrs(k,lev),refgrid).ne.0)then    
                    nsqrs(k,lev)=nsqrs(k,lev)+1                   ! add blocks -> blocks become smaler
                    goto 100
                 else
                 endif
           elseif(360./nsqrs(k,lev).lt.refgrid)then               ! blocks < reference grid
101             if(mod(refgrid,360./nsqrs(k,lev)).ne.0)then
                    nsqrs(k,lev)=nsqrs(k,lev)-1                   ! subtract blocks -> blocks become bigger
                    goto 101
                 else
                 endif
           endif
         endif

         if(mod(nsqrs(k,lev),2).ne.0) stop "ERROR: param_adpx_angle_depth: nsqrs has to be even, stopping ..."         

         nsqtot(k,lev)=n1layer(lev)
         n1layer(lev)=n1layer(lev)+nsqrs(k,lev)

      enddo
      
      nsqtot(nlatzones0+1,lev)=n1layer(lev)
      print*,'Pixels at resolution level 1: ',n1layer(lev) 

! -------------------------------------------------------------
! define finer grids compatible with rough one  

      do lev=2,nlev
         ifa=2**(lev-1)
         size=refgrid/ifa
         nlatzones=180./size
         colat=-size/2.
         n1layer(lev)=0
         do k=1,nlatzones
            k0=((k-1)/2)+1
            nsqrs(k,lev)=nsqrs(k0,lev-1)*2
            nsqtot(k,lev)=n1layer(lev)
            n1layer(lev)=n1layer(lev)+nsqrs(k,lev)
         enddo
         nsqtot(nlatzones+1,lev)=n1layer(lev)
         print*,"Pixels at resolution level",lev,n1layer(lev)
         write(82,*),n1layer(lev)," pixels at resolution level",lev
  enddo

! -------------------------------------------------------------
! read hitcount values

  open(27,file=trim(namehits))
  allocate(htc(4,n1layer(nlev)*npar*nlayi))
  allocate(htcad(4,nadmax))
  ! initialize hitcount matrix
  htcad=0  
  do i=1,n1layer(nlev)*npar*nlayi
     if(mod(i,30000).eq.0) print*, 'Read',i,'hitcount values'
     read(27,*)htc(1:4,i)
  enddo
  print*, i, " hitcount values read"
  close(27)

! -------------------------------------------------------------
! define adaptive grid
!
! We start at the lowest level. If the hitcount is below the threshold 
! the pixel can't be refined and the pixel just gets a new pixel number.
! Then the next level is treated in the same way, but just pixel which 
! do not have a new pixelnumber yet. At the end the pixel in the highest
! level, which do not have a new number yet also get a new number.
!


!  if (writegrid.eq.1) then
     open(100,file=trim(namegrid))
     open(200,file=trim(namegrid)//".lay")
!  end if
  n1layer_adpx=0

! -------------------------------------------------------------
! Loop over layers
  do layer0=1,nlayi
     if (adaptmode.eq.1) then
        print*,"Using a full 3d variable grid"
        layer=layer0
        if(layer0.eq.1) then ! never use the first layer
           layer=2
        endif
     elseif (adaptmode.ge.2) then
        print*,"Defining grid according to layer",adaptmode
        layer=adaptmode
     else
        stop 'ERROR: param_adpx_angle_depth: adaptivity mode chosen incorrectly, stopping ...'
     end if

     print*,"Layer ",layer0
     inew0=0
     ipxnew=0
     ipxnewtmp=0
     htctmp(:)=0
     htctmp_each(:,:)=0
     htctmp_layernumber(:)=0
     htctmp_each_layernumber(:,:)=0


     ! * * * * * REGIONALIZATION * * * * * * *
     ! keep big pixel outside Europe or USA     
     if (outsideREG==1) then
        print*,"Preconditioning grid to keep large pixels outside EU"
        lev=1
        ifa=1
        size=refgrid/ifa
        nlatzones=180./size
        do k=1,nlatzones
           nsqrstmp(k)=nsqrs(k,lev)
           nsqtottmp(k)=nsqtot(k,lev)
        enddo
        nsqtottmp(nlatzones+1)=nsqtot(nlatzones+1,lev)
        do ipx=1,n1layer(1)
           call rang(ipx,xlamin,xlamax,xlomin,xlomax,nsqrstmp,nsqtottmp,nlatzones,size,nlatzomax)
           if (xlamin>=70.or.xlamax<=30.or.(xlomax<=335.and.xlomin>=45)) then ! outside the region we focus on
              write(77,*) xlamin, xlamax, xlomax, xlomin       
              call pxwithin(ipx,iwithin) ! check for lowest order pixel inside 
              htctmp(:)=0                ! count together hits of all of the smallest pixel in this pixel                
              do k=1,4**(nlev-lev)       ! loop over number of smallest pixel in pixel of actual level
                 htctmp(:)=htctmp(:)+htc(:,iwithin(k))
              enddo
              ipxnew=ipxnew+1
              do k=1,4**(nlev-lev)       ! all small pixel inside are asigned to be a big pixel
                 inew0(iwithin(k))=ipxnew
              enddo
              htcad(:,ipxnew)=htctmp(:)       ! hitcount of new pixel is hitcount of all small pixel together
              write(100,'(8I8)') ipx,lev,ipxnew,sum(htctmp),htctmp(1:4)
           end if
        end do  
     else if (outsideREG==2) then
        print*,"Preconditioning grid to keep large pixels outside US"
        lev=1
        ifa=1
        size=refgrid/ifa
        nlatzones=180./size
        do k=1,nlatzones
           nsqrstmp(k)=nsqrs(k,lev)
           nsqtottmp(k)=nsqtot(k,lev)
        enddo
        nsqtottmp(nlatzones+1)=nsqtot(nlatzones+1,lev)
        do ipx=1,n1layer(1)
           call rang(ipx,xlamin,xlamax,xlomin,xlomax,nsqrstmp,nsqtottmp,nlatzones,size,nlatzomax)
           if (xlamin>=70.or.xlamax<=10.or.(xlomax<=180.and.xlomin>=315)) then ! outside the region we focus on
              write(77,*) xlamin, xlamax, xlomax, xlomin       
              call pxwithin(ipx,iwithin) ! check for lowest order pixel inside 
              htctmp(:)=0                ! count together hits of all of the smallest pixel in this pixel                
              do k=1,4**(nlev-lev)       ! loop over number of smallest pixel in pixel of actual level
                 htctmp(:)=htctmp(:)+htc(:,iwithin(k))
              enddo
              ipxnew=ipxnew+1
              do k=1,4**(nlev-lev)       ! all small pixel inside are asigned to be a big pixel
                 inew0(iwithin(k))=ipxnew
              enddo
              htcad(:,ipxnew)=htctmp(:)       ! hitcount of new pixel is hitcount of all small pixel together
              write(100,'(8I8)') ipx,lev,ipxnew,sum(htctmp),htctmp(1:4)
           end if
        end do  
     end if
     print*,"... done with Regionalization!"


     ! * * * * REST OF GLOBE * * * *
     do lev=1,nlev-1	! rough -> fine  !--------do lev
        do ipx=1,n1layer(lev)	                !-------do ipx

           ! in iwithin the number of the finest pixels within the actual one is stored
           call pxwithin_j(nsqrs,nsqtot,nlatzomax,nlev,lev,lev+1,ipx,iwithin_next)   !next level pixel inside this
           call pxwithin(ipx,iwithin) !smallest pixel inside this

           if (inew0(iwithin(1)).eq.0) then	! only check for hitcount if not already new index assigned for lower level	    
              ! there are only pixels of the smallest size inside the actual one        
              if ((lev+1)==nlev) then		!--------if
                 htctmp(:)=0
                 htctmp_each(:,:)=0
                 htctmp_layernumber(:)=0
                 htctmp_each_layernumber(:,:)=0
                 do k=1,4 !hitcount for each wedge all together and each px
                    htctmp(:)=htctmp(:)+htc(:,iwithin(k))
                    htctmp_layernumber(:)=htctmp_layernumber(:)+htc(:,iwithin(k)+n1layer(nlev)*(layer-1)) 
                    htctmp_each(:,k)=htc(:,iwithin(k))
                    htctmp_each_layernumber(:,k)=htc(:,iwithin(k)+n1layer(nlev)*(layer-1))                      
                 end do
              else
                 htctmp(:)=0
                 htctmp_each(:,:)=0
                 htctmp_layernumber(:)=0
                 do k=1,4**(nlev-lev)
                    htctmp(:)=htctmp(:)+htc(:,iwithin(k))  
                    htctmp_layernumber(:)=htctmp_layernumber(:)+htc(:,iwithin(k)+n1layer(nlev)*(layer-1)) ! hitcount in layer for which parameterization is defined
                 end do
                 ! define hitcount for subpixel               
                 htctmp_each(:,:)=0
                 htctmp_each_layernumber(:,:)=0
                 do k=1,4 ! smaller pixel
                    call pxwithin_j(nsqrs,nsqtot,nlatzomax,nlev,lev+1,nlev,iwithin_next(k),iwithin_tmp) !iwithin contains numbers of tiny pixel
                    do i=1,4**(nlev-lev-1)    !sum up hitcount of tiny pixel
                       htctmp_each(:,k)=htctmp_each(:,k)+htc(:,iwithin_tmp(i))
                       htctmp_each_layernumber(:,k)=htctmp_each_layernumber(:,k)+htc(:,iwithin_tmp(i)+n1layer(nlev)*(layer-1))
                    end do
                 end do
              end if                             !------if
              
              tmp=0
              do k=1,4
                 tmp=htctmp_each(k,1)+htctmp_each(k,2)+htctmp_each(k,3)+htctmp_each(k,4)
                 if (htctmp(k)/=tmp) stop "ERROR: param_adpx_angle_depth: undefined problem, stopping ..."
              end do
             
              ! check if all of the small pixel fulfill the threshold criterion,
              ! otherwise assign the pixel to this level
              wedgecnt(:)=0
              do j=1,4 ! check all px
                 do i=1,4 !check how many wedges are covered properly according to threshold value
                    if(htctmp_each_layernumber(i,j)>=ithres(lev)) wedgecnt(j)=wedgecnt(j)+1
                 end do
              end do
              do i=1,4	
                 if (wedgecnt(i)>4) stop "problem wedgecnt"
              end do
              do j=1,4 
                 if (wedgecnt(j)<3) then !if less than 3 wedges are covered no split up possible, assign px to actual level 
                    ipxnew=ipxnew+1
                    do k=1,4**(nlev-lev) !assign 1 big pixel 
                       inew0(iwithin(k))=ipxnew
                    end do
                    htcad(1:4,ipxnew)=htctmp(1:4)
                    write(100,'(8I8)') ipx,lev,(n1layer_adpx_total+ipxnew),sum(htctmp),htctmp(1:4)                                        
                    goto 33 
                 end if
              end do
33         end if 		!-------endif px not assigned
        end do			!-------enddo px
        ipxnewtmp=ipxnew
     end do			!-------enddo lev	 
     
   
     ! assign unassigned pixels to level nlev (highest resolution)
     do ipx=1,n1layer(nlev)
        htctmp(1:4)=0
        if (inew0(ipx).eq.0) then
           htctmp(1:4)=htc(1:4,ipx)
           ipxnew=ipxnew+1
           inew0(ipx)=ipxnew
           htcad(:,ipxnew)=htc(1:4,ipx)
           summe=htc(1,ipx)+htc(2,ipx)+htc(3,ipx)+htc(4,ipx)
           write(100,'(8I8)') ipx,nlev,(n1layer_adpx_total+ipxnew),sum(htctmp),htctmp(1:4)        
        end if
     end do
     
     ! Some verbose
     print*,ipxnew-ipxnewtmp," inv coeffs at lev ",nlev
     print*,ipxnew," inv coeffs total "
     n1layer_adpx(layer0)=ipxnew     
     
!     if (writegrid.eq.1) then
        write(200,'(8I8)') ipxnew
!     end if
     write(*,*) 'Parameters: ',n1layer_adpx(layer0),ipxnew
     write(*,*) '* * * * * * * * * * * * * * * * * '


     ! "project"/store current grid on its actual layer
     do ipx1layer=1,n1layer(nlev)
        ipx=ipx1layer+n1layer(nlev)*(layer0-1)                 ! pixelnumber in finest grid
        inew(ipx)=inew0(ipx1layer)+n1layer_adpx_total !*(layer-1)      ! pixelnumber in adaptive grid depending of number in finest grid
        htcad(:,inew(ipx))=htcad(:,inew(ipx))+htc(1:4,ipx) !!!! DONT KNOW IF THIS IS CORRECT
     end do

     n1layer_adpx_total=n1layer_adpx_total+ipxnew
      
  end do ! ----> end of loop over layers  $

  ! "project"/store current grid to the 2nd parameter
  do ipar=1,(npar-1)
     do ipx1layer=1,n1layer(nlev)*nlayi
        ipx=ipx1layer+n1layer(nlev)*nlayi*ipar                 ! pixelnumber in finest grid
        inew(ipx)=inew(ipx1layer)+n1layer_adpx_total*ipar      ! pixelnumber in adaptive grid depending of number in finest grid
        htcad(:,inew(ipx))=htcad(:,inew(ipx))+htc(1:4,ipx)
     end do
  end do
  write(*,*) 'Grid copied for second parameter!'  
  close(100)
  close(200)

  npar_adpx=n1layer_adpx_total*npar ! total number of parameter   
  npar_reg=n1layer(nlev)*npar*nlayi

  
  ! check if too many parameter/pixel assumed 

  print*,""
  print*,""
  print*,'Total pars',npar_adpx,'max pars',nadmax
  print*,""
  print*,""
   
  if (npar_adpx.gt.nadmax) stop "ERROR: param_adpx_angle_depth: too many parameters, stopping ..." 
             


! -------------------------------------------------------------
! save adaptive grid info to file  

  open(99,file=trim(nameadpx))
  open(101,file=trim(namehtcadgrid))
  
  do ipx=1,n1layer(nlev)*nlayi*npar !npx_adpx !n1layer(nlev)
     write(99,*)ipx,inew(ipx)
     tmp=htcad(1,inew(ipx))+htcad(2,inew(ipx))+htcad(3,inew(ipx))+htcad(4,inew(ipx))
     write(101,*)ipx,tmp,htcad(1:4,inew(ipx))
  end do

  deallocate (angle_hitcount,htctmp_each)
  
  close(99)
  close(101)
  
  open(102,file=trim(namenumberadpx))
  write(102,*) n1layer_adpx, npx_adpx
  close(102)
  
  ! print*, "open ", namegrid
  open(100,file=trim(namegrid),status='old')
  open(103,file=trim(namegrid)//".gmt")
  open(104,file=trim(namegrid)//".sh")
  
  lev0=0
  do i=1,n1layer_adpx_total
     read(100,*) ipx,lev,ipxnew,htctmp     
     if(lev.ne.lev0)then
        ifa=2**(lev-1)
        size=refgrid/ifa
        nlatzones=180./size
        do k=1,nlatzones
           nsqrstmp(k)=nsqrs(k,lev)
           nsqtottmp(k)=nsqtot(k,lev)
        enddo
        nsqtottmp(nlatzones+1)=nsqtot(nlatzones+1,lev)
        lev0=lev
     end if

     call rang(ipx,xlamin,xlamax,xlomin,xlomax,nsqrstmp,nsqtottmp,nlatzones,size,nlatzomax)
     
     write(103,"(a1)")">"
     write(103,*)xlomin,xlamin
     write(103,*)xlomin,xlamax
     write(103,*)xlomax,xlamax
     write(103,*)xlomax,xlamin
    
     write(104,"(I7,4(1X,F10.3))")ipxnew,xlamin,xlamax,xlomin,xlomax

  end do
   
  close(100)
  close(103)
  close(104)

  if (writegrid.eq.1) then
     ! Check if reading is possible and store layer by layer gmt file
     open(430,file=trim(namegrid)//".lay")
     open(440,file=trim(namegrid)//".gmt")
     do ixx=1,nlayi
        read(430,*),layers
        write(i3,"(i3.3)")ixx
        ! separate grid in layers and export them to output folder
        open(450,file='./grid.layer.'//i3//'.gmt')
        do jxx=1,layers
           read(440,*)  ,dmy3
           write(450,*) ,dmy3
           read(440,*)  ,xx1,xx2
           write(450,*) ,xx1,xx2
           read(440,*)  ,xx1,xx2
           write(450,*) ,xx1,xx2
           read(440,*)  ,xx1,xx2
           write(450,*) ,xx1,xx2
           read(440,*)  ,xx1,xx2
           write(450,*) ,xx1,xx2
        end do
        close(450)
     end do
     close(430)
     close(440)
  end if

  write(*,*) 'Grid info files saved!'

  goto 12 
9999 continue 
  if (info.ne.0) write(*,*) 'Controlled termination on error.'

12 continue
   info=0
   return
 end subroutine param_adpx


!------------------------------------------------------------
 subroutine rang(nsq,xlamin,xlamax,xlomin,xlomax,nsqrs,nsqtottmp,nlatzones,size,nlatzomax)
   ! finds the coordinate range of square number 'nsq'
   dimension nsqrs(nlatzomax),nsqtottmp(nlatzomax+1)
   lazone=2
   do while (nsq.gt.nsqtottmp(lazone))
      lazone=lazone+1
   enddo
   lazone=lazone-1
   nnsq=nsq-nsqtottmp(lazone)
   xlamin=90.-lazone*size
   xlamax=xlamin+size
   grsize=360./nsqrs(lazone)
   xlomax=nnsq*grsize
   xlomin=xlomax-grsize
   return
 end subroutine rang
  

!------------------------------------------------------------
subroutine pxwithin(ipx,iwithin)
  ! finds indexes of level-nlev pixels within level-n pixel ipx
  use project_adpx_module
  implicit none
  integer :: lzin
  integer, dimension(4**nlev)  :: iwithin
  integer                      :: ipx, ipx0, ipxin0,k,ifa,n
  ifa=2**(nlev-lev)				! = actual pixel size / smallest pixel size
  k=1
  do while(ipx.gt.nsqtot(k,lev))	      ! search for number of pixel till actual latitudinal zone in actual pixel level
     k=k+1
  enddo
  ipx0=ipx-nsqtot(k-1,lev)-1		      ! number of pixel before actual pixel counted from beginning of this latitudinal zone
  ipxin0=ipx0*ifa					! same but for finest pixel
  n=1
  do i=1,ifa					! determines pixel within this pixel
     lzin=(k-2)*ifa+i
     do j=1,ifa
        iwithin(n)=nsqtot(lzin,nlev)+ipxin0+j
        n=n+1
     enddo
  enddo
  return	
end subroutine pxwithin


!------------------------------------------------------------
subroutine pxwithin_j(nsqrs,nsqtot,nlatzomax,nlev,lev,levin,ipx,iwithin)
  ! finds indexes of level-nlev pixels within level-n pixel ipx
  dimension nsqrs(nlatzomax,nlev),nsqtot(nlatzomax+1,nlev)
  dimension iwithin(4**nlev)
  ifa=2**(levin-lev)                  ! = actual pixel size / smaller pixel size (from level for which indicees should be defined)
  k=1
  do while(ipx.gt.nsqtot(k,lev))      ! search for number of pixel till actual latitudinal zone in actual pixel level
     k=k+1
  enddo
  ipx0=ipx-nsqtot(k-1,lev)-1		! number of pixel before actual pixel counted from beginning of this latitudinal zone
  ipxin0=ipx0*ifa				! same but for smaller pixel
  n=1							! determines pixel within this pixel
  do i=1,ifa					! sample latidudinal zones		
     lzin=(k-2)*ifa+i
     do j=1,ifa				! sample all small pixel in actual latidudinal zone
        iwithin(n)=nsqtot(lzin,levin)+ipxin0+j
        n=n+1
     enddo
  enddo
  return
end subroutine pxwithin_j
