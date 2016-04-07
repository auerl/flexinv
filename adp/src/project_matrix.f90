!--------------------------------------------------------------
! from L. Auer 2012 and J. Schaefer 2009
    

subroutine project_matrix

  use project_adpx_module
  implicit none
  real,dimension(:),allocatable     :: values,valuesa ! values in 1 row of original and adaptive matrix
  real,dimension(m)                 :: t              ! rhs
  integer*4,dimension(:),allocatable:: indx, indxa    ! index of values in 1 row of matrix  
  integer                           :: icol           ! index for loop over rows
  integer                           :: mpoin_row,mpoin_now,mpoin_old  
  integer                           :: i2,k,nrec,n
  integer                           :: rec_written    
  integer                           :: rec_new        
  integer                           :: natamax 
  real                              :: t1,t2,t3,t4,t5,t6
  
! -------------------------------------------------------------
! allocate memory

  call cpu_time(t1)
  natamax=npar_adpx*(npar_adpx+1)
  natamax=natamax/2

  print*,"npar_adpx = ",npar_adpx," natamax = ",natamax

  allocate(indx(npar_reg),values(npar_reg))
  allocate(indxa(npar_adpx),valuesa(npar_adpx))

  if (writeata.eq.1) then
     allocate(ata(natamax),atd(npar_adpx))
  end if

  mpoin_old=0
  mpoin_row=0
  rec_written=0

  call cpu_time(t2)
  print*, "Time for allocating arrays", t2-t1


! -------------------------------------------------------------
! open files for matrix with adaptive grid

  open(1,file=trim(namexxxad),access='direct',recl=4,form='unformatted',status='replace')
  open(2,file=trim(nameindad),access='direct',recl=4,form='unformatted',status='replace')
  open(3,file=trim(namerhsad),status='replace')
  open(4,file=trim(namepoiad),status='replace')

! -------------------------------------------------------------
! read data vector

  open(24,file=trim(namerhs),status='old')
  do i=1,m
        read(24,*,err=154)t(i)          
  end do
  close(24)

! -------------------------------------------------------------
! projecting row by row

  open(23,file=trim(namepoi),status='old')
  k=0
  do icol=1,m 
      ! pointer
      read(23,*,err=153)mpoin_now         ! Nr. values/row, as many as pixel crossed
      mpoin_row=mpoin_now - mpoin_old     ! actual value
      ! read row      
      if (icol>=r1.and.icol<=r2) then
         open(21,file=trim(namexxx),status='old',access='direct',recl=4,form='unformatted')
         open(22,file=trim(nameind),status='old',access='direct',recl=4,form='unformatted')
         ! values, index
         jj=0
         do nrec=mpoin_old+1, mpoin_now
            jj=jj+1
            read(21,rec=nrec,err=156)values(jj) ! values: raypath through this pixel
            values(jj)=values(jj)*relwei
            read(22,rec=nrec,err=155)indx(jj)   ! index in row = index of corresponding pixel
            if(indx(jj).gt.npar_reg)then
               print *,"ERROR: undefined voxel index, stopping ... ", npar_reg, indx(jj)
               stop
            endif
         enddo
         if (jj>n1layer(nlev)) then
            print*, "ERROR: error while reading the matrix, stopping ..."
            stop
         end if         
         close(21)
         close(22)
         if (icol==1)  then
            call cpu_time(t3)
            print*, "Time for reading row", t3-t2            
         end if
         if (rec_written/=k) then
            print*, "ERROR: unknown problem, stopping ..."
            stop
         end if
         
         ! * * * * MAIN PROJECTION LOOP * * * * *
         ! projecting row by row to adaptive grid         
         rec_new=0
         do i=1,jj
            if(values(i).ne.0.)then                            ! after working with value it is set to 0
               k=k+1                                           ! record number to write the adapted value   
               rec_new=rec_new+1                               ! counts matrix entries added in this step
               valuesa(rec_new)=values(i)                      ! value of new a is same as for old a
               indxa(rec_new)=inew(indx(i))
               values(i)=0.
               do i2=1,jj                                      ! are more mat-entries for same pixel? add them
                  if(inew(indx(i)).eq.inew(indx(i2)))then
                     valuesa(rec_new)=valuesa(rec_new)+values(i2)
                     values(i2)=0.                             ! to be sure this value is not used again
                  endif
               enddo
            endif
         end do
         
         if (rec_new/=k-rec_written) then
            print*, "ERROR: something went wrong, stopping ..."
            stop
         end if         
         if (icol==1)  then
            call cpu_time(t4)
            print*, "Time for projecting row", t4-t3
         end if

         ! * * * * * * * * * * * * * * * 
         ! write ata version of new mat
         if(writeata.eq.1)then
            print*,"Starting to write ata matrix!"
            call contribution_ata(valuesa,indxa,ata,t(icol),atd,npar_adpx,natamax,rec_new)
         endif

         call cpu_time(t5)  
         if(writeata.eq.1)then
            if (icol<4)  then
               print*, "Time to write ata", t5-t4
            end if
         end if

         ! * * * * * * * * * * * * * * * 
         ! write row of new matrix
         do i=1,rec_new
            write(2,rec=rec_written+i)indxa(i)           ! index
            write(1,rec=rec_written+i)valuesa(i)         ! values
         enddo
         rec_written=rec_new+rec_written 
         write(4,*)  rec_written                         ! pointer: records/row, added
         if (mod(icol,1000).eq.0) print*,icol," rows projected"
      endif 
      mpoin_old = mpoin_now
   enddo 
   close(23)

   ! Write RHS
   do i=1,m
      write(3,*)t(i) ! should be equal to original file
   end do

   close(1)
   close(2)
   close(3)
   close(4)
  
   ! Save ata matrix
   if(writeata.eq.1)then
      open(31,file=trim(nameata),access='direct',recl=4*natamax,form='unformatted')
      open(32,file=trim(nameatd),access='direct',recl=4*npar_adpx,form='unformatted')     
      call cpu_time(t2)
      write(31,rec=1) (ata(i),i=1,natamax) 
      write(32,rec=1) (atd(i),i=1,npar_adpx) 
      call cpu_time(t3)        
      print*, "ata and atd matrix written in ", t3-t2, "s"
      deallocate(indx,indxa,values,valuesa,ata,atd)
   elseif(writeata.eq.0)then
      deallocate(indx,indxa,values,valuesa)
   endif

   close(31)
   close(32)
   info=0
   return

! -------------------------------------------------------------
! error messages 

153   print*,"ERROR while reading pointer", icol
      stop
154   print*,"ERROR while reading data vector"
      stop
155   print*,"ERROR while reading index",  icol, jj, nrec
      stop
156   print*,"ERROR while reading value", icol, jj, nrec
      stop

end subroutine project_matrix    
        

! -------------------------------------------------------------
subroutine contribution_ata(row,index,ata,d,b,n,nata,nz)
! given a row of a and corresponding datum augments ata accordingly
! j it's unimportant for ata which row, just index of value in row important 
  real ata(nata),b(n),row(n),d
  integer*4 index(n)
  do i=1,nz
     do j=i,nz
        if(index(j).ge.index(i))then
           ind=(((index(j)-1)*index(j))/2)+index(i)
        else
           ind=(((index(i)-1)*index(i))/2)+index(j)
        endif
        ata(ind)=ata(ind)+row(i)*row(j)
     enddo
     ! add to rhs vector the contribution of j-th row:
     b(index(i))=b(index(i))+(row(i)*d)
  enddo
  return
end subroutine contribution_ata
