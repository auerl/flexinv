      subroutine get_dampfac(io,nl,dmscale,dmnorm,dmgrad)
      dimension dmscale(1),dmnorm(1),dmgrad(1)
c
c read a damping factor file to get the factors for block diagonal damping
c for each sub-layer of the model, for example, 14 basis function means
c 14 sublayers (or radial splines) for the inversion.
c Input:	io ---- file number
c		nl ---- number of expected layers
c Output:	
c		dmwate ----  scaling factor for damping matrix
c		dmnorm ----  norm damping factor
c		dmgrad ----  gradient damping factor
c
      character*256  dampfile
      write(*, "('enter the file that contains the damping factors')")
      write(*, "('contents of the file:')")
      write(*, "('-----------------------------------------------------------------------------------')")
      write(*, "('first line[mandatory, 4 entries]:  #_of_layers  multiplier   %norm   %gradient  ')")
      write(*, "('other lines [optional,4 entries]:  spec_layer#  %_of_multplier   %norm   %gradient ')")
      write(*, "('-----------------------------------------------------------------------------------')")
      read(*,"(a)") dampfile
      open(io, file=dampfile, iostat=ios)
      if(ios.ne.0) then
	 stop 'problem with this file!'
      endif
      do i=1, nl
	 dmscale(i)=0.0
	 dmnorm(i)=0.0
	 dmgrad(i)=0.0
      enddo
      read(io,*) nlayer, scale,facnorm,facgrad
      if(nlayer.ne.nl) stop 'number of expected radial parm not consistent!'
      do i=1, nl
	 dmscale(i)=scale
	 dmnorm(i)=facnorm
	 dmgrad(i)=facgrad
      enddo
10    read(io, *, end=20) ilayer,sc,facnorm,facgrad
	 dmscale(ilayer)=dmscale(ilayer)*sc
	 dmnorm(ilayer)=facnorm
	 dmgrad(ilayer)=facgrad
	 goto 10
20    continue
      close(io)
      write(*,*)
      print*, 'SCALING AND DAMPING FACTORS:'	
      print*, '-----------------------------------------------------'	
      print*, 'layer#  scaling_fac  norm_damp   grad_damp'	
      print*, '-----------------------------------------------------'	
      do idmp=1, nl
	 print*, idmp, dmscale(idmp), dmnorm(idmp), dmgrad(idmp)
      enddo
      return
      end
