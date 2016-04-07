      subroutine getseg_path(rpath,nth,rknt,nknt, kradpath)
c...Find the radial b-spline segment for a given raypath.
c
      implicit double precision (a-h,o-z)
      dimension rpath(1),rknt(1)
      dimension kradpath(1)  ! spline segments for the path

      kradpath(i)=0
      do i=1, nth
	 call getbsreg(rpath(i), rknt, nknt, ii)
	 kradpath(i) = ii
      enddo
      return
      end
