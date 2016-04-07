      subroutine join_surfacewave_ata(atd1,ata1,atd,ata,natd,weight,ngpt,nrad,imat)
c
c add two surfacewave ATA matrices
c Add two ATA (particularly made for dispersion data) and form a 
c large ATA with both one matrix at top left corner and the other 
c at the bottom right 2.c
      dimension atd1(1),atd(1),ata1(1),ata(1)

      nparm=ngpt*nrad
      if(imat.ne.1) goto 200
      write(*, *) 'processing --- adding Love wave ATA matrix to the final ATA...'
      write(*, *) 'position --- left top half of matrix'
      k=1
      do i=1, nparm
	 atd(i)=atd(i)+atd1(i)*weight
	 do j=1, i
		ata(k)=ata(k)+ata1(k)*weight
		k=k+1
	 enddo
      enddo
      goto 1000
200   write(*, *) 'processing --- adding Rayleigh wave ATA matrix to the final ATA...'
      write(*, *) 'position --- right bottom half of matrix'
      nn=1
      ind=1
      k=1
      l=1
      do i=1, nparm
	i1=i+nparm
	atd(i1)=atd(i1)+atd1(i)*weight
	iblock=(i-1)/ngpt+1
	i0=mod(i-1, ngpt) + 1 + (ki-1)*ngpt
	do j=1, i
		jblock=(j-1)/ngpt+1
		ata(ki)=ata(ki)+ata1(k)*weight
	enddo
      enddo
1000  return
      end
