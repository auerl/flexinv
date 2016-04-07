      subroutine addmat_spline(atd1,ata1,atd,ata,natd,weight)
c
c add two ATA matrices
c
      dimension atd1(1),atd(1),ata1(1),ata(1)

      write(*, *) 'processing --- adding ATA matrix to the final ATA...'
      k=1
      do i=1, natd
	 atd(i)=atd(i)+atd1(i)*weight
	 do j=1, i
		ata(k)=ata(k)+ata1(k)*weight
		k=k+1
	 enddo
      enddo
      return
      end
