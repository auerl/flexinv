      subroutine stat_lower(ata,atd,numatd,amax,bmax,suma,sumb,ifaniso)
      dimension ata(1),atd(1)
c
c  this checks the diagonal of a lower triangular ATA matrix
c  Variables with subscript 0 (e.g. suma0) are the results of 
c  the first half of anisotropic matrices.
c  
      amax0 = -999999999999.
      amax = -999999999999.
      amin0 = 9999999999.0
      amin = 9999999999.0
      bmax0 = -9999999999999.
      bmax = -9999999999999.
      suma0 = 0.0
      suma = 0.
      sumb0 = 0.
      sumb = 0.
      ind = 0
      numatd0 = numatd/2
      do  i=1,numatd
         do j=1,i
	    ind=ind+1
	    if(i.eq.j) then
			add = ata(ind)
                        if(add.lt.0.) stop 'error 1 in subroutine maxrms, negative diagonal'
            		if(amax.lt.add) amax = add
            		if(amin.gt.add) amin = add
			if(ifaniso.eq.1.or.ifaniso.eq.3) then
				if(i.le.numatd0) suma0=suma0+add
            			if(amax0.lt.add) amax0 = add
            			if(amin0.gt.add) amin0 = add
			endif
                        suma=suma+add
            endif
	 enddo
         add=atd(i)**2
         sumb=sumb+add
         bmax=amax1(bmax,add)
      enddo
      if(ifaniso.eq.1.or.ifaniso.eq.3) then
	write(*, "('****diagonal elements:')")
      	write(*,"('First  half --- max=',e12.4, ' min= ', e12.4,' mean=',e12.4)")
     .  	amax0, amin0, suma0/float(numatd0)
      	write(*,"('Second half --- mean=',e12.4)") (suma-suma0)/float(numatd0) 
      endif
      suma=suma/float(numatd)
      sumb=sqrt(sumb/float(numatd))
      bmax=sqrt(bmax)
      write(*,"('Total diagonal elements:  max=',e12.4, ' min= ', e12.4,' mean=',e12.4)")
     .  amax,amin, suma
      write(*,"('Total b-vector elements:  max=',1pe12.4,'  rms=',e12.4)")
     .  bmax,sumb
      return
      end
