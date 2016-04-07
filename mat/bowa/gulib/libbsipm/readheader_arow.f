      subroutine readheader_arow(io,ngpt,nknot,nparm,n2d)
c
c Read A matrix parameterization information.
c Input:	
c	  io  ---   unformatted A matrix file
c Output:
c	 ngpt ---   number of horizontal splines
c	 nknot ---   number of radial splines
c	 nparm ---   number of anisotropic parameters
c	 n2d ---   layers of 2D topography
c
      real xlat(1),xlon(1),rknot(1)
      integer nnb(1),near(1,1)

      read(io) ngpt,nknot,nparm,n2d
      return
      end
