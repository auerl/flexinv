      subroutine writekin(io,outfl,lmax,kmax,iradfun,dd,dxdp,ttime,slow, iflag)
      implicit double precision (a-h,o-z)
c  write out kIn and other info
c  Input:	io --- file number
c	     outfl --- output file name
c	      lmax --- maximum spherical harmonics (number of splines)
c	      kmax --- maximum radial order
c	       sum --- kernel kIn
c	   iradfun --- choice of radial function
c	        dd --- distance
c             dxdp --- distance derivative with ray parameter
c	     ttime --- travel time
c	      slow --- slowness of the ray
c	     iflag --- if last record = 2, else = 1
c  Output:     outfl (binary) of all above info
      character*256 outfl
      common/kernel/sum(362,20,10),ifsplit,ifdiff
      
c      write(*,*) 'dd = ', dd, ' time = ', ttime, ' slowness=', slow, ' dxdp=', dxdp
      
      if(iflag.eq.2)then
         write(io) -1.0,0.0,0.0,0.0
         close(io)
         return
      endif
      write(io) sngl(dd), sngl(ttime), sngl(slow), sngl(dxdp)
c      write(*,*) sngl(dd), sngl(ttime), sngl(slow), sngl(dxdp)
      do m=1, 5
      	  write(io) ((sngl(sum(l,k+1,m)),l=1,2*lmax+1),k=0,kmax)
      enddo
      iii=1
      return
      end
