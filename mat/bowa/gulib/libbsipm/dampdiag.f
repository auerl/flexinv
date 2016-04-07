	subroutine dampdiag(ata,atd,ngpt,nst0,n2d,ispl,nsplit,aal,aau,aal2,aau2,xs,xt,ifaniso)
c
c damp the diagonal of a matrix.  This routine also adds starting and target model
c to the ATD.
c Input:	
c	  ata, atd  ---   as the names suggest
c	  ngpt      ---   number of horizontal splines
c	  nst0      ---   number of radial splines for an isotropic inversion
c	  n2d       ---   number of topography layers (2 for 410+660)
c	  ispl      ---   spline number that represents the spliting radius, 8 for U6L8
c			  keep in mind the model assumes CMB-->SURFACE order
c	  aal       ---   constant damping for the lower mantle of the first structure
c	  aau       ---   constant damping for the upper mantle of the first structure
c	  aal2      ---   constant damping for the lower mantle of the second structure (if anisotropic)
c	  aau2      ---   constant damping for the upper mantle of the second structure (if anisotropic)
c	  xs        ---   starting model
c	  xt        ---   target model
c	  ifaniso   ---   1=anisotropic inversion   2=isotropic inversion
c Output:
c	  ata, atd  ---   ATA and ATD matrices after damping
c
	dimension ata(1),atd(1),xs(1),xt(1)

	if(ifaniso.ne.1) then
		nstruct=nst0+n2d
	else
		nstruct=nst0*2+n2d
	endif
	numatd=nstruct*ngpt
	indlm=ngpt*ispl
	indum=ngpt*nst0
	indlm2=ngpt*nst0+indlm
	indum2=ngpt*nst0*2
	print*, '+++++ diagonal, start and target model damping ++++'
	print*, 'indlm=', indlm
	print*, 'indum=', indum
	print*, 'indlm2=', indlm2
	print*, 'indum2=', indum2
	print*, 'number of split =', nsplit
	
	ind=1
	if(nsplit.eq.0) then
c ... continuous model parameterization
		do i=1,numatd
			if(i.le.indum) then
				add=aau
			else
				add=aau2
			endif
	   		ata(ind)=ata(ind)+add
	   		atd(i)=atd(i)+add*(xt(i)-xs(i))
	   		ind=ind+i+1
		enddo
	else if(nsplit.eq.1) then
c ... split model parameterization, must be in the order of CMB-SURFACE
		do i=1,numatd
	   		if(i.le.indlm) then
	   			add=aal
	   		else if(i.le.indum) then
	   			add=aau
	   		else if(i.le.indlm2) then
	   			add=aal2
	   		else if(i.le.indum2) then
	   			add=aau2
	   		endif
	   		ata(ind)=ata(ind)+add
	   		atd(i)=atd(i)+add*(xt(i)-xs(i))
	   		ind=ind+i+1
		enddo
	else
		write(*,"('Multiple split has not been implemented!')")
		stop
	endif
	return
	end
