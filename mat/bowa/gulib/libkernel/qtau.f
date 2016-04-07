      subroutine qtau(x,iq,q1)
c  returns the radicand of qtau for radius x
c  xa, xc, xn, xl, eta are A, C, N, L, eta respectively in PREM,
      implicit double precision (a-h,o-z)
      logical isotrp
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/parm/rho,vpv,vsv,vph,vsh,eta,s(5),r
      common/parm1/xa,xc,xn,xl,xf,tdif(5)
      common/isot/isotrp

      y=x/rnorm
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)

      xc=rho*vpv*vpv	! C (check PREM paper)
      xl=rho*vsv*vsv	! L
      xa=xc		! A = C if isotropic
      if(.not.isotrp) then
        xa=rho*vph*vph	! A
        xn=rho*vsh*vsh	! C
        xf=eta*(xa-2.d0*xl)	! F = eta*(A-2*L)
c imod=1 for isotropic, and imod=0 for anisotropic
c if isotropic, does the block code.
        if(imod.ne.0) then
          xkapa=(4.d0*xa+xc+4.d0*xf-4.d0*xn)/9.d0
          xmu=(xa+xc-2.d0*xf+5.d0*xn+6.d0*xl)/15.d0
          xa=xkapa+4.d0*xmu/3.d0
          xc=xa  
          xf=xkapa-2.d0*xmu/3.d0
          xn=xmu 
          xl=xn
          vph=dsqrt(xa/rho)
          vpv=vph
          vsh=dsqrt(xn/rho)
          vsv=vsh
          eta=1.d0
        end if
      end if
c*** avoid the center of the earth!
      if(x.lt.1.d-3) x=1.d-3
c this is p/r term in the square root
      px2=(p/x)**2
      if(itype.eq.1) then
c**** SH
c  reflect at fluid boundary
        if(vsv.eq.0.d0) then
          q1=-1.d0
          return
        end if
        if(isotrp) then
          q1=1.d0/(vsv*vsv)-px2
        else 
          q1=(1.d0-vsh*vsh*px2)/(vsv*vsv)
        end if
      else if(vsv.eq.0.d0) then
c  reflect SV at fluid boundary
        if(itype.eq.3) then
           q1=-1.d0
           return
        end if
c*** P in fluid
        q1=1.d0/(vpv*vpv)-px2
      else
c*** P/SV in solid
        if(isotrp) then
          if(itype.eq.2) q1=1.d0/(vpv*vpv)-px2
          if(itype.eq.3) q1=1.d0/(vsv*vsv)-px2
        else
          aa=0.5d0/(vsv*vsv)
          bb=0.5d0/(vpv*vpv)
          s(1)=aa+bb
          s(2)=aa-bb
          s(3)=(xa*xc-xf*xf-2.d0*xf*xl)/(2.d0*xc*xl)
          s(4)=s(3)*s(3)-xa/xc
          s(5)=0.5d0*(1.d0+xa/xl)/(vpv*vpv)-s(1)*s(3)
          r=sqrt((s(4)*px2+2.d0*s(5))*px2+s(2)*s(2))
          if(itype.eq.2) q1=s(1)-s(3)*px2-r
          if(itype.eq.3) q1=s(1)-s(3)*px2+r
        end if
      end if
c	print*, xl, vpv, vsv, xc, xl, rho, q1
      return
      end
