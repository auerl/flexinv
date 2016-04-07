c
c ***
c *** 
c *** 
c
      subroutine integrate(layers,xivoigt)
      common /plevel/ iprtlv
      include 'gemodl.h'
      
      real ll,nn
c
      parameter (maxker=100)
      parameter (maxpar=8)
      common/mdpert/omegaref,gvelref,numker,numpar,parm(maxker,maxpar),ixi
      double precision temp(maxker,maxpar),acum(maxker,maxpar),add
      character*80 xivoigt
c
      dimension eif(nknts,6)
      dimension layers(30)
c
      dimension q(3),qp(3)
      dimension a(3),b(3)
      equivalence (q(1),ur),(q(2),vr),(q(3),pr)
     #          ,(qp(1),upr),(qp(2),vpr),(qp(3),ppr)
c 
c---- Gaussian integration (Abramowitz and Stegun, p. 916).
c
      dimension xi(5),wt(5)
      data xi/-0.90617 98459 38664,-0.53846 93101 05683
     #       , 0.00000 00000 00000, 0.53846 93101 05683
     #       , 0.90617 98459 38664/
      data wt/0.23692 68850 56189, 0.47862 86704 99366
     #       ,0.56888 88888 88889, 0.47862 86704 99366
     #       ,0.23692 68850 56189/
c
      data twth/0.66666 66666 66667/
      data thrd,fot/.33333 33333 33333,1.33333 33333 33333/
c
c---- copy the eigenfunctions into eif for convenience
c
      do i=1,nknts
        eif(i,1)=u(i)
        eif(i,2)=up(i)
        eif(i,3)=v(i)
        eif(i,4)=vp(i)
        eif(i,5)=p(i)
        eif(i,6)=pp(i)
      enddo

c
      nbot=noc+1
      ntop=moho
      if(jcom.eq.0) return
      fl=float(lord)
      fl3=fl*(fl+1.)
      sfl3=sqrt(fl3)
      con1=fl3-3.
      con2=fl3*(fl3-2.)+3.*(4.*fl3-8.)
      con3=fl3*(fl3-2.)
      con4=3.*(4.*fl3-8.)
      omnd=wcom/wn
      omn2=omnd**2
c
      j1=1
      j2=3
      if (jcom.eq.2) then
        j1=2
        j2=2
      endif 
c
c---- set up the accumulated perturbations
c
      numpar=8
      omegaref=wcom
      gvelref=cgp
      do ipar=1,numpar
        do iker=1,numker
          acum(iker,ipar)=0.
        enddo
      enddo
c
c---- loop 100 - main loop from bottom to top of model
c
      do 100 iq=nbot,ntop
        do ipar=1,numpar
          do iker=1,numker
            temp(iker,ipar)=0.
          enddo
        enddo
c
        if(iq.eq.n) goto 100
c
c---- call corfac to correct for dispersion
c
        iiq=iq
        call corfac(iiq,wcom,jcom,xac,xf,xln)
c
        iq1 = iq+1
        r1  = r(iq)
        r2  = r(iq1)
        hn  = r2-r1
        hnh = hn*.5
        if (hn.lt.1.e-5) goto 100
        hr  = 1./hn
        hsq = hr*hr
        hcu = hr*hsq
c
c---- cubic spline of eigen function
c
        do i=j1,j2
          i1 = 2*i - 1
          i2 = i1 + 1
          a(i) = (eif(iq,i2)+eif(iq1,i2))*hsq
     #      + 2.*(eif(iq,i1)-eif(iq1,i1))*hcu
          b(i) = -(2.*eif(iq,i2)+eif(iq1,i2))*hr
     #      - 3.*(eif(iq,i1)-eif(iq1,i1))*hsq
        enddo
c
        gd = fot*rho(iq)
        if (iq.ne.1) then
          gd = 4.*rho(iq)-2.*g(iq)/r(iq)
        endif
        gd1 = 4.*rho(iq1)-2.*g(iq1)/r(iq1)
        ag  = (gd+gd1)*hsq+2.*(g(iq)-g(iq1))*hcu
        bg  = -(2.*gd+gd1)*hr-3.*(g(iq)-g(iq1))*hsq
c
c---- loop - 5 point integration between r1 and r2
c
        do il=1,5
          t   = .5 * hn * (xi(il)+1.)
          rr  = r1 + t
          rr2 = rr * rr
c
c---- evaluate the eigenfunctions at depth rr
c
          do i=j1,j2
            i1=2*i-1
            i2=i1+1
            q(i)=eif(iq,i1)+t*(eif(iq,i2)+t*(b(i)+t*a(i)))
            qp(i)=(eif(iq,i2)+t*(2.*b(i)+t*3.*a(i)))*rr
          enddo
c
c---- evaluate the elastic parameters and density at depth rr
c
          aa=xac*(acon(iq)+t*(qacon(1,iq)+t*(qacon(2,iq)+t*qacon(3,iq))))
          cc=xac*(ccon(iq)+t*(qccon(1,iq)+t*(qccon(2,iq)+t*qccon(3,iq))))
          ff=xf*(fcon(iq)+t*(qfcon(1,iq)+t*(qfcon(2,iq)+t*qfcon(3,iq))))
          ll=xln*(lcon(iq)+t*(qlcon(1,iq)+t*(qlcon(2,iq)+t*qlcon(3,iq))))
          nn=xln*(ncon(iq)+t*(qncon(1,iq)+t*(qncon(2,iq)+t*qncon(3,iq))))
          rrho=rho(iq)+t*(qrho(1,iq)+t*(qrho(2,iq)+t*qrho(3,iq)))
c
c---- evaluate additional parameters at rr
c
          gr=g(iq)+t*(gd+t*(bg+t*ag))
          etan=ff/(aa-2.*ll)
          delv2=(ll+nn)/rrho
          delvp2=(aa+cc)/rrho
c
          delvph2=(aa)/rrho
          delvpv2=(cc)/rrho
          delvsh2=(nn)/rrho
          delvsv2=(ll)/rrho

c---- voigt and anisotropy as defined by Ferrera et al. (2009)
          vsvoigt2=(ll+nn)/(2*rrho)
          vpvoigt2=(aa+cc)/(2*rrho)
          anis=(delvsh2-delvsv2)/(2*vsvoigt2)
          anip=(delvph2-delvpv2)/(2*vpvoigt2)

c
c---- calculate the kernels, i.e. abar0, cbar0, etc.
c
          if(jcom.ne.2) then
            f=2.*ur-fl3*vr
            xx=vpr-vr+ur
            xx2=xx*xx
            vr2=vr*vr
            ur2=ur*ur
            rka0=f*f
            rkc0=upr*upr
            rkf0=2.*upr*f
            rkl0=fl3*xx2
            rkn0=-rka0+con3*vr2
            t1=(-omn2*vr2*rr+vr*(2.*pr+gr*ur))*rr
            rkr0=((8.*rrho-omn2)*ur2*rr-gr*ur*(f+2.*ur)+2.*ur*ppr)*rr+fl3*t1
          else
            xx=vpr-vr
            xx2=xx*xx
            vr2=vr*vr
            rkl0=fl3*xx2
            rkn0=con3*vr2
            t1=-omn2*vr2*rr2
            rkr0=fl3*t1
          endif
c
c---- accumulate the total perturbation: cmb to moho
c
          radius4=6371.*rr
          
          if (xivoigt.eq.'n') then ! vsh,vsv,vph,vpv mode
             if (iq.gt.noc.and.iq.lt.moho) then
                do iker=1,numker
                   call radial_basis(radius4,'F',iker,value4,istatus,layers)

                   if(istatus.ne.0) then
                      write(6,*) radius4,rr,iker
                      stop 'error in radial_basis'
                   endif
c---- rho
                   if (jcom.eq.2) then
                      add=wt(il)*0.5*rrho*value4*hnh*
     #                (rkr0+rkl0*delvsv2+rkn0*delvsh2)
                   else
                      add=wt(il)*0.5*rrho*value4*hnh*(rkr0+delvpv2*rkc0+
     #                (rka0+etan*rkf0)*delvph2+
     #                (rkl0-2.*etan*rkf0)*delvsv2+rkn0*delvsh2)
                   endif
                   temp(iker,1) = temp(iker,1) + add
c---- Vph
                   if (jcom.eq.2) then
                      add=0.
                   else
                      add=wt(il)*rrho*delvph2*value4*(rka0+etan*rkf0)*hnh
                   endif
                   temp(iker,2) = temp(iker,2) + add
c---- Vpv
                   if (jcom.eq.2) then
                      add=0.
                   else
                      add=wt(il)*rrho*delvpv2*value4*(rkc0)*hnh
                   endif
                   temp(iker,3) = temp(iker,3) + add
c---- Vsh
                   if (jcom.eq.2) then
                      add=wt(il)*rrho*delvsh2*value4*(rkn0)*hnh
                   else
                      add=wt(il)*rrho*delvsh2*value4*(rkn0)*hnh
                   endif
                   temp(iker,4) = temp(iker,4) + add
c---- Vsv
                   if (jcom.eq.2) then
                      add=wt(il)*rrho*delvsv2*value4*(rkl0)*hnh
                   else
                      add=wt(il)*rrho*delvsv2*value4*(rkl0-2.*etan*rkf0)*hnh
                   endif
                   temp(iker,5) = temp(iker,5) + add
c---- eta
                   if (jcom.eq.2) then
                      add=0.
                   else
                      add=wt(il)*0.5*rrho*(delvph2-2.*delvsv2)*value4*etan*(rkf0)*hnh
                   endif
                   temp(iker,6) = temp(iker,6) + add
c---- qmu
                   if (jcom.eq.2) then
                      add=0.
                   else
                      add=0.
                   endif
                   temp(iker,7) = temp(iker,7) + add
c---- qkappa
                   if (jcom.eq.2) then
                      add=0.
                   else
                      add=0.
                   endif
                   temp(iker,8) = temp(iker,8) + add
                enddo
             endif
          elseif(xivoigt.eq.'y') then! voigt average/xi modus
             if (iq.gt.noc.and.iq.lt.moho) then
                do iker=1,numker
                   call radial_basis(radius4,'F',iker,value4,istatus,layers)

                   if(istatus.ne.0) then
                      write(6,*) radius4,rr,iker
                      stop 'error in radial_basis'
                   endif
c---- rho
                   if (jcom.eq.2) then
                      add=wt(il)*0.5*rrho*value4*hnh*
     #                (rkr0+rkl0*delvsv2+rkn0*delvsh2)
                   else
                      add=wt(il)*0.5*rrho*value4*hnh*(rkr0+delvpv2*rkc0+
     #                (rka0+etan*rkf0)*delvph2+
     #                (rkl0-2.*etan*rkf0)*delvsv2+rkn0*delvsh2)
                   endif
                      temp(iker,1) = temp(iker,1) + add
c---- vpvoigt kernel
                   if (jcom.eq.2) then
                      add=0.
                   else
                      add=wt(il)*rrho*vpvoigt2*value4*(anip*rka0-anip*rkc0+
     #                etan*anip*rkf0+etan*rkf0+rka0+rkc0)*hnh
                   endif
                   temp(iker,2) = temp(iker,2) + add
c---- anip kernel
                   if (jcom.eq.2) then
                      add=0.
                   else
                      add=wt(il)*rrho*vpvoigt2*value4*(rka0+
     #                etan*rkf0-rkc0)*hnh
                   endif
                   temp(iker,3) = temp(iker,3) + add
c---- vsvoigt kernel
                   if (jcom.eq.2) then
                      add=wt(il)*rrho*vsvoigt2*value4*(anis*rkn0-anis*rkl0+
     #                rkn0+rkl0)*hnh
                   else
                      add=wt(il)*rrho*vsvoigt2*value4*(anis*rkn0-anis*rkl0+
     #                2.*etan*anis*rkf0-2.*etan*rkf0+rkn0+rkl0)*hnh
                   endif
                   temp(iker,4) = temp(iker,4) + add
c---- anis kernel
                   if (jcom.eq.2) then
                      add=wt(il)*rrho*vsvoigt2*value4*(rkn0-rkl0)*hnh
                   else
                      add=wt(il)*rrho*vsvoigt2*value4*(rkn0-2.*etan*rkf0-rkl0)*hnh
                   endif
                   temp(iker,5) = temp(iker,5) + add
c---- eta
                   if (jcom.eq.2) then
                      add=0.
                   else
                      add=wt(il)*0.5*rrho*(delvph2-2.*delvsv2)*value4*etan*(rkf0)*hnh
                   endif
                   temp(iker,6) = temp(iker,6) + add
c---- qmu
                   if (jcom.eq.2) then
                      add=0.
                   else
                      add=0.
                   endif
                   temp(iker,7) = temp(iker,7) + add
c---- qkappa
                   if (jcom.eq.2) then
                      add=0.
                   else
                      add=0.
                   endif
                   temp(iker,8) = temp(iker,8) + add
                enddo
             endif
          else
             print*,"ERROR: I don't know which physical parameterization to use"
             stop
          endif
       enddo
c
c
c
       do ipar=1,numpar
          do iker=1,numker
             acum(iker,ipar)=acum(iker,ipar)+temp(iker,ipar)
          enddo
       enddo
c
  100 continue
c
      do ipar=1,numpar
        do iker=1,numker
          parm(iker,ipar)=acum(iker,ipar) * wcom
        enddo
      enddo
      return
      end
