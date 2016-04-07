      subroutine fqanis(x,iq,vals)
c   builds up the discretized integrands (with appropriate
c   renormalization) for use in subroutine integr.
c              x = modified radius
c             iq = layer number
c         vals() = returned values of desired integrands
c
c  (rewritten to handle vertical incidence, p=0)
      implicit double precision (a-h,o-z)
      dimension vals(*),zk(4)
      logical isotrp
      common/isot/isotrp
      common/parm/rho,vpv,vsv,vph,vsh,eta,s(5),r
      common/parm1/xa,xc,xn,xl,xf,tdif(5)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/xnorm/xnorm
      common/dx/dxtdp,fxtpp,ifl
c*** here, y is true radius, xnorm is turning point and x=sqrt(y-xnorm) 
c*** so y=xnorm+x*x and dy=2xdx, this change of variable eliminates 
c*** turning point singularities
      y=x*x+xnorm
      call qtau(y,iq,q)
      q=dsqrt(q)
      p2=p/(y*y)
c
c  The following is p*p/(r*r), JG. 2002
c
      py2=p2*p
c*** delta
      q1=p2/q
      f=1.d0
      if(.not.isotrp) then
        if(itype.eq.1) then
          f=xn/xl
        else
          if(vsv.ne.0.d0) then
            a=(s(4)*py2+s(5))/r
            if(itype.eq.3) a=-a
            f=s(3)+a
          end if
        end if
      end if
      q1=q1*f
      vals(1)=2.d0*x*q1
      if(itim.eq.1) return
c*** dX/dp
      q2=f*(1.d0+p*q1/q)/(y*y*q)
      if(.not.isotrp) then
        if(itype.ne.1)then
          a=2.d0*py2*(s(4)*s(2)*s(2)-s(5)*s(5))/(q*(r**3)*(y**2))
          if(itype.eq.3) a=-a
          q2=q2+a
        end if
      end if
      if(ifl.ge.0) q2=q2-fxtpp*dxtdp/(2.d0*x**3)
      vals(2)=2.d0*x*q2
      if(itim.eq.2) return
c*** time and tstar     
      q3=q+p*q1
      vals(3)=2.d0*x*q3
      z=y/rnorm
      call getq(z,iq,qalpha,qbeta)
      if(itype.eq.2) then
        vals(4)=vals(3)*qalpha
      else
        vals(4)=vals(3)*qbeta
      end if
      if(ider.eq.0) return
c*** derivs
c  the tdifs are the derivs wrt vph,vpv,vsh,vsv,eta
      call tder(y,q)
c  we want to integrate these over each layer multiplied
c  by the appropriate powers of radius for a prem 
c  parameterization -- this will need to be changed for
c  other parametizations
cc      zk(1)=1.d0
cc      do k=2,4
cc        zk(k)=z*zk(k-1)
cc      enddo
cc      do 10 j=1,5
cc      do 10 k=1,4
cc        ind=4*(j-1)+k+4
cc   10   vals(ind)=2.d0*x*tdif(j)*zk(k)
c  the following is used as a test for 3D structure using spherical harmonics
c  indices 1 to 4 are reserved for time, etc..
c
      do 10 j=1,5
        ind=j+4
   10   vals(ind)=2.d0*x*tdif(j)
      return
      end
