      subroutine tder(x,qq)
c*** computes the derivatives of qtau wrt vph,vpv,vsh,vsv,eta
c*** at constant p, they are returned in array tdif in common parm1
c*** note that qtau must have been called first and qq=qtau is input
      implicit double precision (a-h,o-z)
      logical isotrp
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/parm/rho,vpv,vsv,vph,vsh,eta,s(5),r
      common/parm1/xa,xc,xn,xl,xf,tdif(5),tdifxi(5)
      common/isot/isotrp
      dimension sd(5,5)

c *** Testing parameters
      parameter(xivoigt=0)

c *** Initialize array tdif and tdifxi
      do i=1,5
        tdif(i)=0.d0
        tdifxi(i)=0.d0
      enddo

      if(x.lt.1.d-3) return ! Escape seq.
      px2=(p/x)**2          ! Predefined var


c *** If the wave is of type SH do...
      if(itype.eq.1) then
        if(isotrp) then                 ! If isotropic
          tdif(4)=-1.d0/(qq*vsv**3)
          tdif(3)=tdif(4)
        else if(xivoigt.eq.1) then      ! If xi and vsvoigt, we need der wrt to N and L, not vsv
          tdif(3)=-px2/(2*xl*qq)
          tdif(4)=(xn*px2-rho)/(2*xl*xl*qq)
          goto 77
        else                            ! If vsv and vsh, we need der wrt to vsv vsh
          tdif(3)=-px2*vsh/(qq*vsv*vsv) ! wrt vsh
          tdif(4)=-qq/vsv               ! wrt vsv
        end if
c *** Special case of a P wave in a fluid, do...
      else if(vsv.eq.0.d0) then
        tdif(2)=-1.d0/(qq*vpv**3)
        tdif(1)=tdif(2)
c *** All other cases, P or SV phases in a solid, do...
      else
        if(isotrp) then                 ! Special case of isotropy
          if(itype.eq.2) then
            tdif(1)=-1.d0/(qq*vpv**3)
            tdif(2)=tdif(1)
          else
            tdif(3)=-1.d0/(qq*vsv**3)
            tdif(4)=tdif(3)
          end if
        else            
c *** Anisotropic P-SV case. Note that sd() are derivs of s() 
c *** wrt A,C,N,L,F resp similarly rd is deriv of r. Analytical! Should be checked
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

c *** Compute sensitivities for A,C,N,L,F. Note that here basically the chain rule
c *** was applied, this is where the factors (1/2r) and (1/2qq) are coming from 
          do j=1,5
            if(j.ne.3) then ! Note: Derivatives wrt N are always 0, see s fcts in WG1981
              rd=((px2*sd(4,j)+2.d0*sd(5,j))*px2+
     +                2.d0*s(2)*sd(2,j))/(2.d0*r)   ! ERROR? px2*px2 needed? 
              if (itype.eq.3) rd=-rd ! In case of P (itype=3) otherwise take +rd
              tdif(j)=(sd(1,j)-sd(3,j)*px2-rd)/(2.d0*qq) 
            end if
          enddo

          if (xivoigt.eq.0) then! If inversion is NOT in terms of xivoigt
            ! Change from A,C,N,L,F to A,C,N,L,eta
              tdif(1)=tdif(1)+eta*tdif(5)
              tdif(4)=tdif(4)-2.d0*eta*tdif(5)
              tdif(5)=tdif(5)*(xa-2.d0*xl)
            ! Change to vph,vpv,vsh,vsv,eta
              tdif(1)=2.d0*rho*vph*tdif(1)
              tdif(2)=2.d0*rho*vpv*tdif(2)
              tdif(4)=2.d0*rho*vsv*tdif(4)
c              print*,"Kernels for vsv and vsh computed"
           else ! If inversion IS in terms of xivoigt
              goto 77
           end if
        end if
      end if

      return


c *** Compute sensitivities for vsvoigt and xi as linear
c *** combinations of A,C,L,N and F Kernels according to Panning 2006
c *** tdif(1-5) is going to be Kvs Kvp Kxi Kphi Keta
77    tdifxi(1)=2.d0*(tdif(3)+tdif(4)-((2.d0*xl)/(xa-2.d0*xl))*tdif(5))  ! vsvoigt
      tdifxi(2)=2.d0*(tdif(1)+tdif(2)+(xa/(xa-2.d0*xl))*tdif(5))         ! vpvoigt
      tdifxi(3)=(1/(2.d0*xl+xn))*(2.d0*xl*tdif(4)-xn*tdif(3)+(2.d0*xl*xn*tdif(5))/(xa-2.d0*xl))  ! xi 
      tdifxi(4)=(1/(xc+4.d0*xa))*(4.d0*xa*tdif(2)-xc*tdif(1)-(xa*xc*tdif(5))/(xa-2.d0*xl))       ! phi
      tdifxi(5)=tdif(5)         ! eta
c      print*,'Kernels for Vsvoigt and Xi computed!'
      do j=1,5
         tdif(j)=tdifxi(j)
      end do


      return
      end
