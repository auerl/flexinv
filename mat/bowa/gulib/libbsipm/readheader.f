      subroutine readheader(io,ngpt,xlat,xlon,nnb,near,nknot,rknot,
     +			ifsplit,nparm,n2d)
c
c Read A matrix parameterization information.
c Input:	
c	  io  ---   unformatted A matrix file
c Output:
c	 ngpt ---   number of horizontal splines
c	 xlat ---   latitude array of b-splines
c	 xlon ---   longitude array of b-splines
c	 nnb  ---   number of neighbors for a given spline (5 or 6)
c	 near ---   nearest neighbors of splines
c	 nknot ---   number of radial splines
c	 rknot ---   radius of radial splines
c	 ifsplit ---   if model is a split or continuous
c	 nparm ---   number of anisotropic parameters
c	 n2d ---   layers of 2D topography
c
      include "../aniso_spline_parallel.h"
      real xlat(1),xlon(1),rknot(1)
      integer nnb(1),near(mxleny,mxfreq)

      read(io) ngpt,nknot,ifsplit,nparm,n2d
	print*, 'test here', ngpt, nknot, ifsplit, nparm, n2d
      ifreq=sqrt(real((ngpt-2.0)/10.0))+0.45
      do i=1,ngpt
		read(io) iwhatever,xlat(i),xlon(i),nnb(i),
     #		(near(i,j), j=1,ifreq)
      enddo
      read(io) (rknot(j),j=1,nknot)
      return
      end
