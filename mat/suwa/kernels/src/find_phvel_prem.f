c        parameter(pi=3.1415926536)
	pi=2.*asin(1.)
cTEST
	print*,"pi=",pi
      print*,"what period (seconds)?"
      read*,period
      omega=2*pi/period
      gvelr=premgephvelo(2,omega) ! rayleigh
      gvell=premgephvelo(1,omega) ! love
      print*,"love phase velocity=",gvell
      print*,"rayleigh phase velocity=",gvelr
      end
