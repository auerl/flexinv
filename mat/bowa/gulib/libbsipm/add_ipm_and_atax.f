	subroutine add_ipms_and_atax(ata,atd,numatd,ata_add,atd_add,numatd1,weight,coef,itype)
c ... This routine add ata to 
c

	dimension ata(1),atd(1),ata_add(1),atd_add(1),coef(1)

	print*, 'Add initial model(1=yes  0=no)?'
	read(*,*)ifaddinit
	if(ifaddinit.eq.1) then
		call add_atax(ata_add, atd_add, numatd1, coef)
	endif
	numata1=numatd1*(numatd1+1)/2
	numata=numatd*(numatd+1)/2
	print*, 'numatd =', numatd, '  numatd1=', numatd1
	if(itype.eq.1) then
c  P-wave, first part of the anisotropic inversion'
		do i=1,numatd1
           		atd(i)=atd(i)+weight*atd_add(i)
		enddo
		do i=1,numata1
			ata(i)=ata(i)+weight*ata_add(i)
		enddo
	else if(itype.eq.2) then
c  S-wave, second part of the anisotropic inversion
		do i=1,numatd1
			ii=numatd1+i
			atd(ii)=atd(ii)+weight*atd_add(i)
			atd_add(i)=0.0
		enddo
		k=1
		ind=1
		do i=1,numatd
			do j=1, i
				if(i.gt.numatd1.and.j.gt.numatd1) then
					ata(ind)=ata(ind)+weight*ata_add(k)
					if(k.lt.10) then
						print*, 'i, j', i, j, ' ind, k', ind, k
					endif
					k=k+1
				endif
				ind=ind+1
			enddo
		enddo
	else
c  combined anisotropic input, simply add to the final matrices
		do i=1, numatd
			atd(i)=atd(i)+atd_add(i)
		enddo
		do i=1, numata
			ata(i)=ata(i)+ata_add(i)
		enddo
	endif
	return
	end
