      subroutine dxt(x,n,d,f)
c     accepts as input ray parameter (through common/in/),
c     turning radius x (real or fictitious), and layer
c     number n of layer containing x.
c     returns the derivative d of the turning radius w.r.t. p
c
c     the f.p. in the calling routine (integr) are:
c              x -- xnorm
c              n -- n (layer number)
c              d -- dxtdp
c              f -- fxttp
c
      implicit double precision (a-h,o-z)
      logical isotrp
      dimension v(3),dv(3),ds(5),bigv(3,3),bigd(3,3)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/parm/rho,vpv,vsv,vph,vsh,eta,s(5),r
      common/parm1/xa,xc,xn,xl,xf,tdif(5)
      common/isot/isotrp
c
      y=x/rnorm
      call getmod(y,n,rho,vpv,vsv,vph,vsh,eta)
      xc=rho*vpv*vpv
      xl=rho*vsv*vsv
      if(.not.isotrp) then
        xa=rho*vph*vph
        xn=rho*vsh*vsh
        xf=eta*(xa-2.d0*xl)
        if(imod.ne.0) then
          xkapa=(4.d0*xa+xc+4.d0*xf-4.d0*xn)/9.d0
          xmu=(xa+xc-2.d0*xf+5.d0*xn+6.d0*xl)/15.d0
          xa=xkapa+4.d0*xmu/3.d0
          xc=xa
          xf=xkapa-2.d0*xmu/3.d0
          xn=xmu
          xl=xn
          vph=dsqrt(xa/rho)
          vpv=dsqrt(xc/rho)
          vsh=dsqrt(xn/rho)
          vsv=dsqrt(xl/rho)
          eta=1.d0
        end if
      end if
      d=0.
      f=0.
      if(x.lt.1.d-3) return
      px2=(p/x)**2
c
      if(itype.eq.1) then
        if(vsv.eq.0.d0) return
c**** SH
        if(.not.isotrp) then
          dvsh=cder(y,7,n)/vsh
          dvsv=cder(y,3,n)/vsv
          top=px2*xn/(xl*p)
          bot=(px2/x-dvsv/(vsh*vsh)+px2*(dvsv-dvsh))*xn/xl
        else
          top=px2/p
          bot=px2/x-cder(y,3,n)/(vsv**3)
        end if
        d=top/bot
        f=0.5d0*d*dsqrt(2.d0*bot)
      else if(vsv.eq.0.d0) then
c**** fluid
        top=px2/p
        bot=px2/x-cder(y,2,n)/(vpv**3)
        d=top/bot
        f=0.5d0*d*dsqrt(2.d0*bot)
      else
c**** P/SV
        v(2)=vpv
        v(3)=vsv
        dv(2)=cder(y,2,n)
        dv(3)=cder(y,3,n)
        a=dv(2)/(v(2)**3)
        b=dv(3)/(v(3)**3)
        ds(1)=-a-b
        ds(2)=a-b
        if(isotrp) then
          top=2.d0*px2/p
          bot=ds(1)+2.d0*px2/x
          if(itype.eq.2) bot=bot-ds(2)
          if(itype.eq.3) bot=bot+ds(2)
        else
          v(1)=vph
          dv(1)=cder(y,6,n)
          deta=cder(y,8,n)
          do 1 i=1,3
          do 1 j=1,3
    1       bigv(i,j)=v(i)*v(i)/(v(j)*v(j))
          do 2 i=1,3
          do 2 j=1,3
    2       bigd(i,j)=2.d0*bigv(i,j)*(dv(i)/v(i)-dv(j)/v(j))
          aa=0.5d0/(vsv*vsv)
          bb=0.5d0/(vpv*vpv)
          s(1)=aa+bb
          s(2)=aa-bb
          s(3)=(xa*xc-xf*xf-2.d0*xf*xl)/(2.d0*xc*xl)
          s(4)=s(3)*s(3)-xa/xc
          s(5)=0.5d0*(1.d0+xa/xl)/(vpv*vpv)-s(1)*s(3)
          r=sqrt((s(4)*px2+2.d0*s(5))*px2+s(2)*s(2))
c*** ds(1--5) are derivs of s() wrt x (similarly with rdx and r)
          b1=-0.5d0*eta*bigd(1,3)+deta*(2.d0-0.5d0*bigv(1,3))
          a2=-1.d0+eta*(2.d0-0.5d0*bigv(1,3))
          b2=eta*bigd(1,2)+deta*bigv(1,2)
          a3=-2.d0*eta*deta*bigv(3,2)
          a4=2.d0*(1.d0-eta)*(eta*bigd(3,2)+deta*bigv(3,2))
          ds(3)=0.5d0*bigd(1,3)+eta*bigv(1,2)*b1+a2*b2+a3+a4
          ds(4)=2.d0*s(3)*ds(3)-bigd(1,2)
          ds(5)=(0.5d0*bigd(1,3)-dv(2)*(1.d0+bigv(1,3))/vpv)/(vpv*vpv)
     +          -s(1)*ds(3)-s(3)*ds(1)
          a=0.5d0*px2*px2*(ds(4)-4.d0*s(4)/x)
          b=px2*(ds(5)-2.d0*s(5)/x)+s(2)*ds(2)
          rdx=(a+b)/r
          rdp=2.d0*px2*(s(4)*px2+s(5))/(r*p)
          top=2.d0*s(3)*px2/p
          bot=ds(1)+px2*(2.d0*s(3)/x-ds(3))
          if(itype.eq.2)then
            top=top+rdp
            bot=bot-rdx
          else
            bot=bot+rdx
            top=top-rdp
          end if
        end if
        d=top/bot
        f=0.5d0*d*dsqrt(bot)
      end if
      return
      end

