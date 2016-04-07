      subroutine corfac(iq,wcom0,jcom0,xac,xf,xln)
      include 'gemodl.h'
c
c---- correction of elastic parameters for attenutation.
c
      data pi/3.14159265358979/
      data fot/1.33333333333333333333/
c
      fct=2.*alog(.5*wcom0/pi)/pi
      xmu=qshear(iq)*fct
      xln=1.+xmu
      if(jcom0.eq.2) return
      xka=qkappa(iq)*fct
      rt=fot*(acon(iq)+ccon(iq)-2.*fcon(iq)+5.*ncon(iq)+6.*lcon(iq))
     1   /(8.*acon(iq)+3.*ccon(iq)+4.*fcon(iq)+8.*lcon(iq))
      xf=1.+((1.-rt)*xka-.5*rt*xmu)/(1.-1.5*rt)
      xac=1.+(1.-rt)*xka+rt*xmu
      return
      end
