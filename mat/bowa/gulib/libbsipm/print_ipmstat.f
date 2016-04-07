	subroutine print_ipmstat(ata,atd,numatd,ngpt,n2d,ifaniso)
	dimension ata(1), atd(1)

	numatd2d=n2d*ngpt
	numatd3d=numatd-numatd2d
	numatd1=numatd3d	! first half of matrix
	if(ifaniso.ne.0) numatd1=numatd3d/2
	numata=numatd*(numatd+1)/2
	print*, ' '
	print*, ' '
	print *,".......Statistical Analysis of ATA and ATD......."
	sumata=0.0
	sumata1=0.0
	atamin=1.0e+37
	atamax=-1.0e+37
	atamax1=-1.0e+37
	atamax2=-1.0e+37
	k=1
	do i=1, numatd3d
		do j=1, i
			if(i.eq.j) then
				sumata=sumata+ata(k)
				if(atamin.gt.ata(k)) atamin=ata(k)
				if(atamax.lt.ata(k)) atamax=ata(k)
				if(i.eq.numatd1) then
					sumata1=sumata
					atamax1=atamax
				endif
				if(i.gt.numatd1) then
					if(atamax2.lt.ata(k)) atamax2=ata(k)
				endif
			endif
			k=k+1
		enddo
	enddo
	sumatd=0.0
	atdmin=1.0e+37
	atdmax=-1.0e+37
	atdmax1=-1.0e+37
	atdmax2=-1.0e+37
	do i=1,numatd3d
		sumatd=sumatd+atd(i)
		if(i.eq.numatd1) then
			atdmax1=atdmax
			sumatd1=sumatd
		endif
		if(atdmin.gt.atd(i)) atdmin=atd(i)
		if(atdmax.lt.atd(i)) atdmax=atd(i)
		if(i.gt.numatd1) then
			if(atdmax2.lt.atd(i)) atdmax2=atd(i)		
		endif
	enddo
	print *, 'number of element in entire atd = ', numatd
	print *, 'number of element in entire ata = ', numata
	print*,  'minimum  diagonal of entire ata = ', atamin
	print*,  'minimum  element  of entire atd = ', atdmin

	if(ifaniso.eq.0) then
c isotropic matrix
		print*, 'maximum diagonal of ata = ', atamax
		print*, 'average diagonal of ata = ', sumata/numatd3d
		print*, 'minimum of atd = ', atdmin
		print*, 'maximum of atd = ', atdmax
		print*, 'average of atd = ', sumatd/numatd3d
		print*, '..................................................'
	else
c anisotropic matrix
		print*, '-------first half of the ATA matrix-------'
		print*, 'maximum diagonal of ata = ', atamax1
		print*, 'average diagonal of ata = ', sumata1/numatd1
		print*, 'maximum of atd = ', atdmax1
		print*, 'average of atd = ', sumatd1/numatd1
		write(*,*)	
		print*, '-------second half of the ATA matrix-------'
		print*, 'maximum diagonal of ata = ', atamax2
		print*, 'average diagonal of ata = ', (sumata-sumata1)/numatd1
		print*, 'maximum of atd = ', atdmax2
		print*, 'average of atd = ', (sumatd-sumatd1)/numatd1
	endif
c
c  2D topography
c
	sumata=0.0
	sumatd=0.0
	itopo=1
	indtopo=numatd3d+1
	do i=indtopo, numatd
		do j=1, i
			if(i.eq.j) then
				sumata=sumata+ata(k)
			endif
			k=k+1
		enddo
		sumatd=sumatd+atd(i)
		if(mod(i,ngpt).eq.0) then
			print*, 'Topography layer', itopo
			print*, '..................................................'
		   	print*, 'average ATA diagonal of Topography =', sumata/ngpt
		   	print*, 'average ATD diagonal of Topography =', sumatd/ngpt
			write(*,*)
			sumata=0.0
			sumatd=0.0
			itopo=itopo+1
		endif
	enddo
	return
	end
