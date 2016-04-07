	subroutine print_ipmavg(ata,atd,ngpt,nst)
c
c print the diagonal average of each block for the purpose of damping.
c Input:	
c	  ata, atd  ---   as the names suggest
c	  ngpt      ---   number of horizontal splines
c	  nst       ---   number of radial parameters (3D+2D)
c
	dimension ata(1),atd(1)

	numatd=nst*ngpt
	ind=1
	sumata=0.0
	sumatd=0.0
	write(*,*)
	print*, '****diagonal averages****'
	do i=1,numatd
	   iblock=(i-1)/ngpt+1
	   sumata=sumata+ata(ind)
	   sumatd=sumatd+atd(i)
	   if(mod(i,ngpt).eq.0) then
		print*, 'layer ', iblock-1, ' ata= ',sumata/ngpt,  ' atd= ', sumatd/ngpt
		sumata=0.0
		sumatd=0.0
	   endif
	   ind=ind+i+1
	enddo
	return
	end
