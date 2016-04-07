      subroutine writeheader(io,mxnode,ifreq,xlat,xlon,nnb,near,
     +			nknot,rknot,ifspltmd,nparm,n2d)
c
c Write A matrix parameterization information.
c
c Input:	
c	  io  ---   unformatted A matrix file
c Output:
c      mxnode ---   maximum number of horizontal splines (required)
c       ifreq ---   frequency of tesselation (6 = 362 splines)
c	 xlat ---   latitude array of b-splines
c	 xlon ---   longitude array of b-splines
c	 nnb  ---   number of neighbors for a given spline (5 or 6)
c	 near ---   nearest neighbors of splines
c	 nknot ---   number of radial splines
c	 rknot ---   radius of radial splines
c     ifspltmd ---   if model is a split or continuous
c	 nparm ---   number of anisotropic parameters
c	   n2d ---   layers of 2D topography
c
      integer nnb(mxnode),near(mxnode,1)
      real    rknot(1),xlat(1),xlon(1)

      ngpt=ifreq*ifreq*10+2
      write(io) ngpt,nknot,ifspltmd,nparm,n2d
	print*, 'ngpt =', ngpt, ' nknot=', nknot
	print*, near(361,1), near(361,2), near(361,5), near(362,5)
	print*, nnb(1), nnb(2), nnb(361), nnb(362)
      do i=1,ngpt
		write(io) i,xlat(i),xlon(i),nnb(i),
     #		(near(i,j), j=1,ifreq)
      enddo
      write(io) (rknot(j),j=1,nknot)
      return
      end
