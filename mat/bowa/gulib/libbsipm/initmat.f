      subroutine initmat(atd,ata,natd)
c
c initializing matrix
c
      dimension atd(1), ata(1)

      k=1
      do i=1, natd
	 atd(i)=0.0
	 do j=1,i
	 	ata(k)=0.0
		k=k+1		
	 enddo
      enddo
      return
      end

