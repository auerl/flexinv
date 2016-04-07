	subroutine damp_norm_grad(ata,atd,ngpt,nst0,n2d,dmscale,
     +				dmnorm,dmgrad,dmat2d,xs,xt,ifaniso)
c
c damp the ATA matrix using.  This routine also adds starting and target model
c to the ATD.  This version also allows gradient damping for 2D topography layers.
c Input:	
c	  ata, atd  ---   as the names suggest
c	  ngpt      ---   number of horizontal splines
c	  nst0      ---   number of radial splines for an isotropic inversion
c	  n2d       ---   number of topography layers (2 for 410+660)
c	  dmscale   ---   array containing weighting factors for the overall smoothing
c	  dmnorm    ---   norm damping factors
c	  dmgrad    ---   gradient damping factors
c	  dmat2d    ---   2D gradient damping matrix for splines
c	  xs        ---   starting model
c	  xt        ---   target model
c	  ifaniso   ---   1=anisotropic inversion   2=isotropic inversion
c Output:
c	  ata, atd  ---   ATA and ATD matrices after damping
c
	dimension ata(1),atd(1),xs(1),xt(1),dmat2d(1)
	dimension dmscale(1),dmnorm(1),dmgrad(1)

	if(ifaniso.ne.1) then
		nstruct=nst0+n2d
	else
		nstruct=nst0*2+n2d
	endif
	numatd=nstruct*ngpt
	nelement=ngpt*(ngpt+1)/2
	k=1
	l=1
	ind=1
	nn=1
	do i=1,numatd
	   iblock=(i-1)/ngpt+1
	   do j=1, i
		jblock=(j-1)/ngpt+1
		if(nn.gt.nelement) then
		    print*, 'damping block ',   iblock
		    nn=1	! restarting the counter
		endif
		if(iblock.eq.jblock) then
c... diagonal blocks
c		   if(iblock.ne.1.and.nn.le.3) then
c			print*, i,j, iblock,nn,ind
c			print*, dmgrad(iblock), dmat2d(nn), dmnorm(iblock), dmscale(iblock)
c		   endif
			add=dmgrad(iblock)*dmat2d(nn)
			if(i.eq.j) then
c... diagonal elements
				add=add+dmnorm(iblock)
			endif
			ata(ind)=ata(ind)+dmscale(iblock)*add
			nn=nn+1
		endif
	   	ind=ind+1
	   enddo
c
c... damping toward target model
	   if(l.gt.ngpt) then
		k=1
		l=1	
	   endif
c	   if(k.le.6) then
c		print*, 'iblock=', iblock,i,k,dmscale(iblock),dmnorm(iblock),dmgrad(iblock),dmat2d(k)
c	   endif
c... damping factors added to the diagonal terms are added to the atd
	   atd(i)=atd(i)+dmscale(iblock)*(dmnorm(iblock)+dmgrad(iblock)*dmat2d(k))*(xt(i)-xs(i))
	   if(k.gt.nelement.or.iblock.gt.nstruct) stop 'error in calculating the damping matrix!'
	   l=l+1
	   k=k+l
	enddo
	return
	end
