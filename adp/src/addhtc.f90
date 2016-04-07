!
! from Ludwig, 01/2013: Program which adds dws/htc files
! and applies depth dependent weighting
!

program addhtc

   implicit none

   character(len=120) :: htcfil
   real,allocatable,dimension(:) :: x11,x22,x33,x44,mm,sme
   real :: x1,x2,x3,x4,dmy,maxim,po,wghtin
   integer :: i,u,k,j,numlay,numpar,o

   print*,""
   print*,""
   print*," * * * * Adding dws/htc files * * * * *"
   print*,""
   print*,""

   write(*,'(A,$)'),"Enter number of Parameters"
   read(*,*),numpar
   write(*,*),numpar
   write(*,'(A,$)'),"Number of layers"
   read(*,*),numlay
   write(*,*),numlay

   allocate(x11(2*numpar*numlay))
   allocate(x22(2*numpar*numlay))
   allocate(x33(2*numpar*numlay))
   allocate(x44(2*numpar*numlay))
!   allocate(mm(4),sme(numlay))

   u=0
   do j=1,numlay
      do i=1,numpar*2
         u=u+1
         x11(u)=0.
         x22(u)=0.
         x33(u)=0.
         x44(u)=0.
      end do
   end do

   do k=1,1000

      u=0
      read(*,*),htcfil ! read name of htc/dws file
      if (htcfil.eq."END") goto 11
      read(*,*),wghtin ! upweight, downweight due to number of observations, normalized b/w 0 and 1

      print*,"Adding ",trim(htcfil)

      open(10,file=trim(htcfil))      

      do o=1,2
      do j=1,numlay
         do i=1,numpar  
            u=u+1
            if (wghtin.eq.0.) then
               po=((real(j)-1.0)/(real(numlay)-1.0)+1.0)*10
               po=1./po
            elseif (wghtin.gt.0.) then
               po=wghtin
            endif
            read(10,*),x1,x2,x3,x4               

            x11(u)=x11(u)+x1*po!*weight
            x22(u)=x22(u)+x2*po!*weight
            x33(u)=x33(u)+x3*po!*weight
            x44(u)=x44(u)+x4*po!*weight
!            sme(j)=sme(j)+(x1+x2+x3+x4)*weight
         end do
      end do
      end do
      close(10)
   end do

11 continue

   print*,"Exporting results ..."

!   mm(1)=maxval(x11(:))
!   mm(2)=maxval(x22(:))
!   mm(3)=maxval(x33(:))
!   mm(4)=maxval(x44(:))
!   maxim=maxval(mm(:))

   u=0
   open(20,file='./dws.all.dat',status='replace')
!   open(30,file='./dws.all.2.dat',status='replace')
!   open(40,file='./dws.all.3.dat',status='replace')
   do o=1,2
   do j=1,numlay
     do i=1,numpar         
         u=u+1
         write(20,*),x11(u),x22(u),x33(u),x44(u)         
!         write(30,*),x11(u)/sme(j),x22(u)/sme(j),x33(u)/sme(j),x44(u)/sme(j)         
!         write(40,*),x11(u)/sme(j)/real(j),x22(u)/sme(j)/real(j),x33(u)/sme(j)/real(j),x44(u)/sme(j)/real(j)         
      end do
   end do
   end do
   close(20)
   close(30)
   close(40)

   print*,"... done!"
   
end program addhtc
