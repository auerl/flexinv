      subroutine getq(y,iq,qalpha,qbeta)
c*** returns qalpha and qbeta at normalized radius y which
c*** shell iq. Note that this assumes that the 5 and 4
c*** model coefficients relate to inverse Qmu and inverse Qkappa
c*** which is the case for aniprmc
c*** this code is a little cavalier about anisotropic layers!
c*** also we assume constant Q in layers (as is done in readmod)
      implicit real*8(a-h,o-z)
      common/coeff/coef(4,8,20)
      vpv=coef(1,2,iq)+y*(coef(2,2,iq)+y*(coef(3,2,iq)+y*coef(4,2,iq)))
      vsv=coef(1,3,iq)+y*(coef(2,3,iq)+y*(coef(3,3,iq)+y*coef(4,3,iq)))
      qka=coef(1,4,iq)
      qmu=coef(1,5,iq)
      qbeta=qmu
      qalpha=qka+4.d0/3.d0*(qmu-qka)*(vsv/vpv)**2
      return
      end

