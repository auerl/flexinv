      double precision function cder(y,m,iq)
c*** finds radial derivative of model parameter 'm' where
c  m=1 rho; m=2 vpv, =3 vsv; =4 Qm; =5 Qk; =6 vph, =7 vsh; =8 eta
      implicit real*8(a-h,o-z)
      common/coeff/coef(4,8,20)
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      cder=(coef(2,m,iq)+y*(2.d0*coef(3,m,iq)+3.d0*coef(4,m,iq)
     +      *y))/rnorm
      return
      end
