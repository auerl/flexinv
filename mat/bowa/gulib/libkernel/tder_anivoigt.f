      subroutine tder(x,qq)
c*** computes the derivatives of qtau wrt vph,vpv,vsh,vsv,eta
c*** at constant p, they are returned in array tdif in common parm1
c*** note that qtau must have been called first and qq=qtau is input
      implicit double precision (a-h,o-z)
      logical isotrp
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/parm/rho,vpv,vsv,vph,vsh,eta,s(5),r
      common/parm1/xa,xc,xn,xl,xf,tdif(5)
      common/isot/isotrp
      dimension sd(5,5)


c this version of tder just computes kernels for absolute
c perturbations in ACLNF and


c *** Initialize array tdif and tdifxi
      do i=1,5
        tdif(i)=0.d0
      enddo

      if(x.lt.1.d-3) return ! Escape seq.
      px2=(p/x)**2          ! Predefined var


c *** If the wave is of type SH do...
      if(itype.eq.1) then

c -> I comment out all isotrp things, as i only do anisotropic invs
c        if(isotrp) then
c          tdif(4)=-1.d0/(qq*vsv**3)
c          tdif(3)=tdif(4)
c        else

c -> parameterize in absolute L N Kernels, not directly the deriv wrt vsh
          tdif(3)=-px2/(2*xl*qq)
          tdif(4)=(xn*px2-rho)/(2*xl*xl*qq)

c          tdif(3)=-px2*vsh/(qq*vsv*vsv)
c          tdif(4)=-qq/vsv
c        end if
      else if(vsv.eq.0.d0) then
c -> not relevant for me as i dont parameterize the core, but could include
c -> derv wrt CNL
c*** P in fluid
        tdif(2)=-1.d0/(qq*vpv**3)
        tdif(1)=tdif(2)
      else
c*** P/SV in solid
c        if(isotrp) then
c          if(itype.eq.2) then
c            tdif(1)=-1.d0/(qq*vpv**3)
c            tdif(2)=tdif(1)
c          else
c            tdif(3)=-1.d0/(qq*vsv**3)
c            tdif(4)=tdif(3)
c          end if
c        else
c*** note that sd() are derivs of s() wrt A,C,N,L,F resp
c*** similarly rd is deriv of r
          sd(1,1)=0.d0
          sd(1,2)=-rho/(2.d0*xc*xc)
          sd(1,4)=-rho/(2.d0*xl*xl)
          sd(1,5)=0.d0
          sd(2,1)=0.d0
          sd(2,2)=-sd(1,2)
          sd(2,4)=sd(1,4)
          sd(2,5)=0.d0
          sd(3,1)=1.d0/(2.d0*xl)
          sd(3,2)=(xf*xf+2.d0*xl*xf)/(2.d0*xl*xc*xc)
          sd(3,4)=(xf*xf-xa*xc)/(2.d0*xc*xl*xl)
          sd(3,5)=-(1.d0+xf/xl)/xc
          sd(4,1)=s(3)/xl-1.d0/xc
          sd(4,2)=2.d0*s(3)*sd(3,2)+xa/(xc*xc)
          sd(4,4)=2.d0*s(3)*sd(3,4)
          sd(4,5)=2.d0*s(3)*sd(3,5)
          sd(5,1)=rho/(2.d0*xl*xc)-s(1)*sd(3,1)
          sd(5,2)=-rho*(1.d0+xa/xl)/(2.d0*xc*xc)
          sd(5,2)=sd(5,2)-s(1)*sd(3,2)-s(3)*sd(1,2)
          sd(5,4)=-rho*xa/(2.d0*xc*xl*xl)-s(1)*sd(3,4)-s(3)*sd(1,4)
          sd(5,5)=-s(1)*sd(3,5)
          do j=1,5
            if(j.ne.3) then
              rd=((px2*sd(4,j)+2.d0*sd(5,j))*px2+
     +                2.d0*s(2)*sd(2,j))/(2.d0*r)
              if (itype.eq.3) rd=-rd
              tdif(j)=(sd(1,j)-sd(3,j)*px2-rd)/(2.d0*qq)
            end if
          enddo
c -> dont anymore change parametization, will do this myself
c*** now we change from A,C,N,L,F to A,C,N,L,eta
c          tdif(1)=tdif(1)+eta*tdif(5)
c          tdif(4)=tdif(4)-2.d0*eta*tdif(5)
c          tdif(5)=tdif(5)*(xa-2.d0*xl)
c*** now we change to vph,vpv,vsh,vsv,eta
c          tdif(1)=2.d0*rho*vph*tdif(1)
c          tdif(2)=2.d0*rho*vpv*tdif(2)
c          tdif(4)=2.d0*rho*vsv*tdif(4)
c        end if -> isotrp endif
      end if

c now change parameterization from A,C,N,L,F to vp,xp,vs,xs,eta
c my parameterization should be in ABSOLUTE Vs and RELATIVE Xp, so
c I dont need to multiply by a model vector in the Xp case!



      return
      end
