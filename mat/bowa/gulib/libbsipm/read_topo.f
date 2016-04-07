      subroutine read_topo(lu,topofile,nnode,bslat,bslon,ntopo,x)
c   this subroutine reads in a B-spline topography model
c   INPUT:    lu  --- b-spline file number
c   OUTPUT:
c	 topofile --- file name
c	  nnode   --- number of horizontal splines (362 for frequency 6)
c	  ntopo   --- number of topography layers
c	  bslat   --- array for spline node lattitude
c	  bslon   --- array for spline node longitude
c	      x   --- array for model coefficients
c
      character*(*) topofile
      logical    exists
      dimension	 x(1),idumbuf(12)
      dimension  bslon(1),bslat(1)
      integer    n, nadj, ncoefleft, ncoefrow

      inquire(file=topofile(1:lnblnk(topofile)),exist=exists)
      if(exists) then
          open(lu,file=topofile(1:lnblnk(topofile)),status='old')
          write(*,*) 'model file=', topofile(1:lnblnk(topofile))
      else
	  write(*,*) 'spline node file do not exist, exit...'
	  stop
      endif
c
c  read number of coef, number of radial knots, number of inversion
c  structure and mask for structure parameters
      read(lu,*) nnode,ntopo
      print*, 'nnode=', nnode, '  ntopo=', ntopo
      if(ntopo.ne.1) stop 'more than 1 layers, error!'
      nfreq = sqrt((real(nnode)-2.0)/10.0)+0.45  ! find out the frequency
c
c  below reads in the tesselation part
c
      do i=1, nnode
           read(lu, *) n, bslat(i), bslon(i),
     +                     nadj, (idumbuf(j),j=1, nfreq)
      enddo
c
c  below reads in the coefficients
c        
      ncoefrow = nnode/5
      ncoefleft = mod(nnode, 5)
      k=1
c
      write(*,*) '2-D B-spline model'
      do i=1, ntopo
           do j=1,ncoefrow
               read(lu, *) x(k),x(k+1), x(k+2), x(k+3), x(k+4)
               k=k+5
           enddo
           read(lu,*) x(k), x(k+1)
           k=k+2
      enddo
10    close(lu)
      return
      end
