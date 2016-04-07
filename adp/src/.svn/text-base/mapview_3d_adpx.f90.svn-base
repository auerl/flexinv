SUBROUTINE mapview_3d_adpx

 !====================================================================!
 !
 !    J. Schaefer 07/2009 from Lapo's mapview_3d.f90, some little changes
 !
 ! from a 3d mantle model stored in lapo's format, to a set of map-files (one per
 ! model layer) in a format compatible with GMT command psxy
 !
 ! uses parameterization defined by param_adpx
 ! use from module  : namegrid, nlay, rbot, rtop, name, average, n1layer_adpx,nlev,n1layer(nlev), imode
 !
 !====================================================================!

  USE inv_adpx_module
  USE plotprep_adpx_module
  IMPLICIT NONE
  REAL,DIMENSION(nlay)        :: ave      ! vector with average values for all layers
  INTEGER                     :: rincr    ! layer distance
  CHARACTER(LEN=80)           :: string   ! just for reading the colorpalette 
  REAL,DIMENSION(100)      	:: v        ! boundary values for colors
  INTEGER,DIMENSION(100,3)    :: irgb     ! colorpalette 
  INTEGER,DIMENSION(3)        :: rgb 
  INTEGER                     :: nint, nend,nl
  CHARACTER(LEN=11)           :: filename3
  REAL                        :: zdepth   ! average depth of actual layer
  INTEGER                     :: ila      ! loop over lateral zones
  INTEGER                     :: rlong, rinlo,rlati                     ! number of pixel in this latitudinal zone and size
  REAL                        :: xlomin, xlomax, xlamin, xlamax   ! coordinates of the pixel
  INTEGER                     :: k
  REAL,DIMENSION(n1layer_adpx):: ccc_layer      ! value for pixel of the layer
  REAL				:: ccc		! one value out of the vector
  REAL                        :: a,b      ! color values
  INTEGER                     :: isq, ipxnew
  CHARACTER(LEN=2)   		:: char 
  REAL                        :: zmin, zmax 	! min/max depth of the actual layer
  INTEGER,DIMENSION(8)		:: point    
  
  !====================================================================!
  ! check if all necessary information available in the module
  write(82,*) "informations used from the module:"
  write(82,*) "# different size pixels:          	 ", nlev
  write(82,*) "grid information: 				 ", namegrid
  write(82,*) "# layers: 					 ", nlayi
  write(82,*) "top and bottom layers: 			 ", rbot, rtop
  write(82,*) "prepare printing for file: 		 ", name
  write(82,*) "remove average y/n: 				 ", average
  write(82,*) "# pixel/ layer in adaptive grid: 	 ", n1layer_adpx
  write(82,*) "# of pixels/layer for different pixel sizes:  ", n1layer(1:nlev)
  write(82,*) "imode:						 ",imode

 !====================================================================!
 ! preparations
 
  open(104,file=trim(namegrid)//".sh")
  open(100,file=trim(namegrid))
  open(101,file="colorpalette.dat")
  open(107,file=trim(layers))
! 
!  ! layer spacing
!   rincr=(rtop-rbot)/nlay
!   print*, "rincr=",rincr


 !====================================================================!
 ! determine averages of all layers

 ! set average for all layers to 0
  do i=1,nlay
     ave(i)=0.
  enddo
  
  open(1,file=trim(name),status='old')     
  print*, "open file ", trim(name)   
  if(average.eq.'y'.or.average.eq.'Y')then
    print*, "calculate average"
    call layav(nlay,ave,n1layer_adpx,rincr,nlev,n1layer(nlev))
    print*, "...done"
  endif

   
  
 !====================================================================!
 ! open files for writing

  if(imode.eq.3)then
     open(unit=3,file='valxyz.dat')
     print *,'writing to valxyz.dat'
  else if(imode.eq.4)then
     open(unit=3,file='centxy.dat')
     print *,'writing to centxy.dat'
  endif

  if(imode.eq.4)then 
     nend = 1
  else
     nend = nlay
  endif



 !====================================================================!
 ! prepare vtk-file
  
  if(imode.eq.5)then          
     open(111,file=trim(name)//'.vtk')
     write(111,*)"# vtk DataFile Version 3.0"
     write(111,*)"vtk output"
     write(111,*)"ASCII"
     write(111,*)"DATASET UNSTRUCTURED_GRID"
     write(111,*)"POINTS ",npx_adpx*8, " float"

     
     
     do nl=1,nend

        ! read layerdepth 
        read(107,*) k, zmax, zmin
        if (k/=nl) stop "problem with layerfile"
        print*,"depth range of layer",nl, ":", zmax, zmin
        
        do isq=1,n1layer_adpx
             read(104,"(I7,4(1X,F10.3))")ipxnew,xlamin,xlamax,xlomin,xlomax      
             write(111,*) xlomin, xlamin,  (-1) * zmin
             write(111,*) xlomax, xlamin,  (-1) * zmin
             write(111,*) xlomin, xlamax,  (-1) * zmin          
             write(111,*) xlomax, xlamax,  (-1) * zmin
             write(111,*) xlomin, xlamin,  (-1) * zmax
             write(111,*) xlomax, xlamin,  (-1) * zmax
             write(111,*) xlomin, xlamax,  (-1) * zmax          
             write(111,*) xlomax, xlamax,  (-1) * zmax
        enddo
        rewind(104)
     enddo   
     rewind(107)
     
     write(111,*)"CELLS", npx_adpx, 9*npx_adpx
     
     point(1)=0
     point(2)=1
     point(3)=2
     point(4)=3
     point(5)=4
     point(6)=5
     point(7)=6
     point(8)=7

     do nl=1,nend
     do isq=1,n1layer_adpx       
        write(111,'(A1,8(1X,I7))') "8",point(1),point(2),point(3),point(4),point(5),point(6),point(7),point(8)  
        do j=1,8
           point(j)=point(j)+8
        enddo   
     enddo
     enddo
     
     write(111,*)"CELL_TYPES", npx_adpx   
     do i=1,npx_adpx
         write(111,*) "11"    
     enddo
     write(111,*)"CELL_DATA" , npx_adpx   
     write(111,*)"SCALARS my_field float"
     write(111,*)"LOOKUP_TABLE default"

   endif
   




 !====================================================================!
 ! loop over layers
 
  open(107,file=trim(layers))
  do nl=1,nend
     print*,"working on layer",nl, "number of pixel in layer:", n1layer_adpx
     
     
    !====================================================================!
    ! read layerdepth 

    read(107,*) k, zmax, zmin
    if (k/=nl) stop "problem with layerfile"
    print*,"depth range of layer",nl, ":", zmax, zmin
 !    zdepth = rincr*(float(nl)-.5)
    zdepth = (zmin + zmax) / 2 
      
      
     if(imode.lt.3)then
        write(filename3,'(a6,i2.2,a2)')'layer_',nl
        open(unit=3,file=filename3)
     endif
     

     
     ! read in all values for that layer ============================!        
     if(imode.ne.4)then
     do isq=1,n1layer_adpx
        read(1,*)k,ccc_layer(isq)
         write(77,*) k,ccc_layer(isq)
         ccc_layer(isq)=ccc_layer(isq)-ave(nl)
!             write(77,*) k,ccc_layer(isq)
!        ccc_layer(isq)=(ccc_layer(isq)-ave(nl))
!        check=k-(nl-1)*k
!        if 
!        write(77,*) isq, k-(nl-1)*n1layer_adpx- n1layer_adpx*nlayi, k
 !       if( k - (nl-1)*n1layer_adpx - n1layer_adpx*nlayi /= isq) stop "PROBLEM!"
!        write(77,*) isq, k-(nl-1)*k-n1layer_adpx*nlayi , k 
!        if( k-(nl-1)*n1layer_adpx-n1layer_adpx*nlayi /= isq) stop "PROBLEM!"
     enddo
     endif
     
    
    !====================================================================!
    ! determine which colorpalette is usefull
    ! read in discrete gmt color palette table 
    
     if(adopt_color=='n' .or. adopt_color=='N'.and. nl>1) then
	  !did already read colorpalette
	  
        if (nl==1 .and. imode == 1) then
           print*, "using colorpalette", colorname  
           open(unit=54,file=colorname,status='old')
           i=1
1          read(54,"(a80)",end=2)string
           if((string(1:1).eq."#").or.(string(1:1).eq."B").or.(string(1:1).eq."F").or.(string(1:1).eq."N")) goto 1
           read(string,*)v(i),(irgb(i,k),k=1,3),v(i+1),(irgb(i+1,k),k=1,3)
          ! write(*,*)i,v(i),(irgb(i,k),k=1,3),v(i+1)
           i=i+1
           goto 1
2          continue 
	     nint=i-1	
	     close(54)             
        endif	  
	  
	  
     
!        
!      elseif(adopt_color=='y' .or. adopt_color=='Y') then
!            print*, "Data ranges for this layer from ",minval(ccc_layer) , "to", maxval(ccc_layer)
!            write(char,'(I2.2)') int(maxval(abs(ccc_layer)))+1
!            write(colorname2,*) trim(colorname)//char//".cpt" 
! 
! 
!         if(imode.eq.1)then
!         print*, "using colorpalette", colorname2
!         write(101,"('+/-', a2,'%')")char
!     
!        !   30 format(1X,'Zeilennummer =',I7,','2X,'Wert =',F12.6)
! 
!            open(unit=54,file=colorname2,status='old')
!            i=1
! 3          read(54,"(a80)",end=4)string
!            if((string(1:1).eq."#").or.(string(1:1).eq."B").or.(string(1:1).eq."F").or.(string(1:1).eq."N")) goto 3
!            read(string,*)v(i),(irgb(i,k),k=1,3),v(i+1),(irgb(i+1,k),k=1,3)
!           ! write(*,*)i,v(i),(irgb(i,k),k=1,3),v(i+1)
!            i=i+1
!            goto 3
! 4          continue 
! 	     nint=i-1	
! 	     close(54)
!         endif
        
      elseif(adopt_color=='y' .or. adopt_color=='Y') then
!           print*, "Data ranges for this layer from ",minval(ccc_layer) , "to", maxval(ccc_layer)


           if (nl<=6) write(char,'(I2.2)') 10
           if (nl>6) write(char,'(I2.2)') 5          
           write(colorname2,*) trim(colorname)//char//".cpt" 


        if(imode.eq.1)then
        print*, "using colorpalette", colorname2
        write(101,"('+/-', a2,'%')")char
    
       !   30 format(1X,'Zeilennummer =',I7,','2X,'Wert =',F12.6)

           open(unit=54,file=colorname2,status='old')
           i=1
3          read(54,"(a80)",end=4)string
           if((string(1:1).eq."#").or.(string(1:1).eq."B").or.(string(1:1).eq."F").or.(string(1:1).eq."N")) goto 3
           read(string,*)v(i),(irgb(i,k),k=1,3),v(i+1),(irgb(i+1,k),k=1,3)
          ! write(*,*)i,v(i),(irgb(i,k),k=1,3),v(i+1)
           i=i+1
           goto 3
4          continue 
	     nint=i-1	
	     close(54)
        endif       
        
        
        
        
        
     endif
    
       
     do isq=1,n1layer_adpx
           ccc=ccc_layer(isq) 
           if (nl==1) write(88,*) isq, ccc
           read(104,"(I7,4(1X,F10.3))")ipxnew,xlamin,xlamax,xlomin,xlomax

!            if(imode.ne.4)then ! read value and remove mean
!               read(1,*)k,ccc
!               ccc=(ccc-ave(nl))
!            endif
           
 !          if (nl<3) WRITE(77,*) nl, isq, ipxnew, ccc
           
 ! imode 1 ===========================================================!             
 ! gmt values           
           if(imode.eq.1)then
              ! determine color
            if(isq==1)  print*, v(1),v(nint+1),j !test
              do j=1,nint
                 if((ccc.ge.v(j)).and.(ccc.le.v(j+1)))then
   if (isq<10.and.nl==1) print*, "1",ccc   !test
                  do k=1,3   
                     ! interpolation of color
                     a=(irgb(j,k)-irgb(j+1,k))/(v(j)-v(j+1))
                     b=irgb(j,k)-(a*v(j))                    
                     rgb(k)=a*ccc+b
			   ! no interpolation of colours
			   ! rgb(k)=irgb(j,k)
                  enddo
                   
                 endif
              enddo
              
              if(ccc.ge.v(nint+1))then
    if (isq<10.and.nl==1) print*, "2",ccc   !test             
                do k=1,3
                  rgb(k)=irgb(nint+1,k)
                enddo
              endif
              
              if(ccc.le.v(1))then
    if (isq<10.and.nl==1) print*, "3",ccc   !test                  
                do k=1,3
                  rgb(k)=irgb(1,k)
                enddo
                
              endif
              
        if (isq<10.and.nl==1) print*, isq, rgb(1:3)          
              
!       print*,ccc,rgb(1),rgb(2),rgb(3)
40            write(3,421)'> -G',rgb(1),'/',rgb(2),'/',rgb(3)
              write(3,*)xlomin,xlamin
              write(3,*)xlomin,xlamax
              write(3,*)xlomax,xlamax
              write(3,*)xlomax,xlamin
              
 ! imode 2 ===========================================================!             
 ! x,y,val             
           else if(imode.eq.2)then
              write(3,*)(xlomin+xlomax)/2.,(xlamin+xlamax)/2.,ccc
              
 ! imode 3 ===========================================================!             
 ! x,y,z val             
           else if(imode.eq.3)then
              write(3,*)(xlomin+xlomax)/2.,(xlamin+xlamax)/2.,zdepth,ccc
              
 ! imode 4 ===========================================================!             
 ! x,y
           else if(imode.eq.4)then
              write(3,*)(xlomin+xlomax)/2.,(xlamin+xlamax)/2.             
              
 ! imode 5 ===========================================================!             
 ! vtk info
           else if(imode.eq.5)then           
              
               write(111,*) ccc           
              
              
              
           endif
         
         enddo


      if(imode.lt.3)then
          close(3)
      endif
      rewind(104)
 

    enddo
!***************************    
    
    if(imode.ge.3)then
       close(3)
    endif
    if(imode.ne.4)then
       close(1)
    endif
421 format(a4,i3.3,a1,i3.3,a1,i3.3)





  go to 11 
9999 CONTINUE ! jump here upon error
  IF (info .NE. 0) WRITE(*,*)'Controlled termination on error.'

11 continue

  info=0
END SUBROUTINE mapview_3d_adpx


!**********************************************************************!
!     SUBROUTINES
!**********************************************************************!

	subroutine layav(nlay,ave,n1layer,rincr,nlev,n1layerfine)
	dimension ave(nlay)
      integer ibla, ipxnew,lev, weight,nlev,n1layerfine,rincr
	open(17,file="vprofile.txt")
	open(18,file="dvprofile.txt")
      print*, "++++++++++++++"

!	rincr=(rtop-rbot)/nlay
  print*, rincr !, rtop, rbot
	do l=1,nlay
	   tot=0.
	   tot2=0.
         icheck2=0
	   do i=1,n1layer
	      read(1,*)ii,value
            read(100,"(4I15)")ibla,lev,ipxnew,ibla
            if(ipxnew/=i)stop "layav has a problem!" 
            weight= 2**(2*(nlev-lev))
	      tot=tot+value*weight
	      tot2=tot2+(value**2)*weight
!	      icheck=(n1layer*(l-1))+i
!	      if(ii.ne.icheck)stop "check check"
            icheck2=icheck2+weight
	   enddo
         rewind(100)
!	   std = sqrt ((n1layer * tot2 - tot**2) / ((n1layer*(n1layer-1))));
         if(n1layerfine.ne.icheck2)stop "layav has a problem with weighting!"  
	   ave(l)=tot/float(n1layerfine)
!       
         print*, rincr, rincr*l
	   print*,'a: ',ave(l),(rincr*(float(l)-.5)),rincr*(l-1),rincr*(l),std
	   write(17,*)">"
	   write(17,*)rincr*(l-1),ave(l)
	   write(17,*)rincr*(l),ave(l)
!	   write(18,*)">"
!	   write(18,*)rincr*float(l-1),std
!	   write(18,*)rincr*float(l),std
	enddo
	close(17)
!	close(18)
	print*,'I computed all the ',nlay,' averages'
	rewind(1)
	return
	end
