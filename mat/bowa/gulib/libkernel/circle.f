      subroutine circle(radi, inum)
      implicit real*8(a-h,o-z)
      data pi/3.14159265358979d0/
      rad=180./pi
      do i=1,361
        ang=i/rad
        write(inum,*) radi*sin(ang),radi*cos(ang)
      enddo
      return
      end
