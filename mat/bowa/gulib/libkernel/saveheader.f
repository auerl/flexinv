      subroutine saveheader(io,ngpt,nknot,npar,n2d)
c
c Output A matrix headers and model parameters information.
c
      implicit double precision (a-h,o-z)

      write(io) ngpt,nknot,npar,n2d
      return
      end
