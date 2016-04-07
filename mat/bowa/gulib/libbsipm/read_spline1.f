      subroutine read_spline1(lu,bsfile,nnode,bslat,bslon,ifspltmd,nradnd,nradnd_lm,
     +		radnd,x)
c   this subroutine reads in a B-spline model
c   INPUT:    lu    --- b-spline file number
c   OUTPUT:
c	   bsfile   --- file name
c	  ifsplit   --- (1. split B-spline model  0. continuous model)
c	    nnode   --- number of horizontal splines (362 for frequency 6)
c	        x   --- array for model coefficients
c	  bslat     --- array for spline node lattitude
c	  bslon     --- array for spline node longitude
c 	   nradnd   --- number of radial knots for B-spline function
c 	nradnd_lm   --- number of radial knots lower mantle, if continuous, then set to 0
c	    radnd   --- radial knot depths
c
      character*(*) bsfile
      dimension	 x(1)
      dimension  bslon(1),bslat(1), radnd(1)
      integer    ifspltmd,nnode,nradnd,nradnd_lm
      integer    nfreq,lu
      integer    dumbuf(12)  !hardwired to frequency 12
      logical    exists	
      integer    n, nadj, ncoefleft, ncoefrow
      character  firstline*80, modeltype*11

      inquire(file=bsfile(1:lnblnk(bsfile)),exist=exists)
      if(exists) then
          open(lu,file=bsfile(1:lnblnk(bsfile)),status='old')
          write(*,*) 'model file=', bsfile(1:lnblnk(bsfile))
      else
	  write(*,*) 'spline node file do not exist, exit...'
	  stop
      endif
c
c  read number of coef, number of radial knots, number of inversion
c  structure and mask for structure parameters
c  if nradnd=14, means model has no crust and 670 portion
c  then add 2 layers of zeros to be consistent with spherical harm models
c
      read(lu, "(a80)") firstline
      read(firstline, "(a11)") modeltype
      if(modeltype.eq."UPPERLOWERM") then
	 ifspltmd =1
         read(firstline,*) modeltype, nnode , nradnd_lm, nradnd
	 write(*,"('SPLIT MODEL')")
	 write(*,"('nnode=',i4,' upperman=',i2,' lowerman=',i2)") nnode,nradnd,nradnd_lm
      else
	 ifspltmd = 0
         read(firstline,*) nnode , nradnd
	 nradnd_lm = 0
	 write(*,"('CONTINUOUS MODEL')")
	 write(*,"('nnode=',i4,' # radial spline=',i2)") nnode,nradnd
      endif
c
c this makes sure that all the components are being used for B-splines
c        
      nfreq = sqrt((real(nnode)-2.0)/10.0)  ! find out the frequency
      write(*,*) 'B-spline frequency = ', nfreq
c
c  below reads in the tesselation part
c
      do i=1, nnode
           read(lu, *) n, bslat(i), bslon(i),
     1                     nadj, (dumbuf(j),j=1, nfreq)
c		write(*,*) 'test', bslat(i), bslon(i)
      enddo
c
c  below reads in the coefficients
c        
      ncoefrow = nnode/5
      ncoefleft = mod(nnode, 5)
      k=1
c
      write(*,*) '3-D B-spline model'
      do i=1, nradnd+nradnd_lm
           read(lu, *) n, radnd(i)
           write(*,*) 'n=', n, ' radius=', radnd(i)
           do j=1,ncoefrow
               read(lu, *) x(k),x(k+1), x(k+2), x(k+3), x(k+4)
               k=k+5
           enddo
           read(lu,*) x(k), x(k+1)
           k=k+2
      enddo
      close(lu)
      return
      end
