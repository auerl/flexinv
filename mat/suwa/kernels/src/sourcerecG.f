      subroutine sourcerecG(omega,lora,stazi,stdip,delta,takeoff,backazi,
     #           ospli0,dospli,nspl,splefs,splefr,maxspl,maxefs,
     #           cecoeff,cecoeffrad,cecoeffns,cecoeffew)
c
c---- this subroutine calculates the product of the source excitation
c---- term and the receiver term
c---- all input angles are in radians
c---- the following assumes that the modes in splefs and splefr are 
c---- normalized in the way Jeroen does it
c
      dimension splefs(maxspl,maxefs,2)
      dimension splefr(maxspl,maxefs,2)
c
      complex cecoeff(10)
      complex cecoeffrad(10)
      complex cecoeffns(10)
      complex cecoeffew(10)
      complex cecoefftran(10)
      complex cecoefflong(10)
c
      complex acc
      parameter (twopi=6.2831853072)
      complex clongfac
c
      stdipdeg=stdip*360./twopi
      stazideg=stazi*360./twopi
c
c---- calculate projection onto seismometer axis
c
      if(abs(stdipdeg).gt.1.) then
        prtr=0.
        prlod=0.
        ifve=1
        if(stdipdeg-90..lt.1.) then
          prlo=-1.
        else if(stdipdeg+90..lt.1.) then
          prlo=1.
        else
          write(6,"('stazi, stdip:',2g15.5)") stazi,stdip
          stop 'peculiar dip of instrument '
        endif
      else
        ifve=0
        prlo=cos(backazi-stazi)
        prtr=sin(backazi-stazi)
        prlod=-sin(backazi-stazi)
        prtrd=cos(backazi-stazi)
      endif
c
c---- we have some conversion from the angles defined in getpaparm
c
      costaz=cos(0.25*twopi-takeoff)
      sintaz=sin(0.25*twopi-takeoff)
      cos2taz=cos(2.*(0.25*twopi-takeoff))
      sin2taz=sin(2.*(0.25*twopi-takeoff))
c
c---- tanfac is still not quite clear to me
c
      tanfac=-1./(tan(delta))
      sinfac=1./(sin(delta))
c
c---- toroidal modes
c
      if(lora.eq.1) then
        xksou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,3,lora))
        wsou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,7,lora))
        wpsou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,8,lora))
        wppsou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,12,lora))
        wddrec=ecbspl(omega,ospli0,dospli,nspl,splefr(1,10,lora))
c
        acc=cmplx(wddrec,0.)
        cecoeff(1)=cmplx(0.,0.)
        cecoeff(2)=acc*cmplx(0.,-0.5*xksou*wsou*sin2taz*prtr)/cmplx(0.,-omega**3)
        cecoeff(3)=acc*cmplx(0.,0.5*xksou*wsou*sin2taz*prtr)/cmplx(0.,-omega**3)
        cecoeff(4)=acc*cmplx((wpsou-wsou)*costaz*prtr,0.)/cmplx(0.,-omega**3)
        cecoeff(5)=acc*cmplx((wpsou-wsou)*sintaz*prtr,0.)/cmplx(0.,-omega**3)
        cecoeff(6)=acc*cmplx(0.,xksou*wsou*cos2taz*prtr)/cmplx(0.,-omega**3)
c
c---- radial derivative
c
        wsourad=(-wsou+wpsou)
        wpsourad=(wppsou)
c
        cecoeffrad(1)=cmplx(0.,0.)
        cecoeffrad(2)=acc*cmplx(0.,-0.5*xksou*wsourad*sin2taz*prtr)/cmplx(0.,-omega**3)
        cecoeffrad(3)=acc*cmplx(0.,0.5*xksou*wsourad*sin2taz*prtr)/cmplx(0.,-omega**3)
        cecoeffrad(4)=acc*cmplx((wpsourad-wsourad)*costaz*prtr,0.)/cmplx(0.,-omega**3)
        cecoeffrad(5)=acc*cmplx((wpsourad-wsourad)*sintaz*prtr,0.)/cmplx(0.,-omega**3)
        cecoeffrad(6)=acc*cmplx(0.,xksou*wsourad*cos2taz*prtr)/cmplx(0.,-omega**3)
c
c---- epicentral derivatives
c
        cecoefftran(1)=cmplx(0.,0.)
        cecoefftran(2)=tanfac*acc*cmplx(0.,xksou*wsou*cos2taz*prtr)/cmplx(0.,-omega**3)
        cecoefftran(3)=tanfac*acc*cmplx(0.,-xksou*wsou*cos2taz*prtr)/cmplx(0.,-omega**3)
        cecoefftran(4)=tanfac*acc*cmplx((wpsou-wsou)*sintaz*prtr,0.)/cmplx(0.,-omega**3)
        cecoefftran(5)=tanfac*acc*cmplx(-(wpsou-wsou)*costaz*prtr,0.)/cmplx(0.,-omega**3)
        cecoefftran(6)=tanfac*acc*cmplx(0.,2.*xksou*wsou*sin2taz*prtr)/cmplx(0.,-omega**3)
        if(prtr.gt.0.00001) then
          do i=1,6
            cecoefftran(i)=cecoefftran(i)+cmplx(sinfac*prtrd/prtr,0.)*cecoeff(i)
          enddo
        endif
c
c---- spheroidal modes
c
      else if(lora.eq.2) then
        xksou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,3,lora))
        usou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,5,lora))
        upsou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,6,lora))
        uppsou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,11,lora))
        vsou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,7,lora))
        vpsou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,8,lora))
        vppsou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,12,lora))
        uddrec=ecbspl(omega,ospli0,dospli,nspl,splefr(1,9,lora))
        vddrec=ecbspl(omega,ospli0,dospli,nspl,splefr(1,10,lora))
c
        if(ifve.eq.1) then
          acc=cmplx(uddrec,0.)/cmplx(0.,-omega**3)
        else
          acc=cmplx(0.,vddrec*prlo)/cmplx(0.,-omega**3)
        endif
c
        cecoeff(1)=acc*cmplx(upsou,0.)
        cecoeff(2)=acc*cmplx((usou-0.5*vsou*xksou
     #           +0.5*xksou*vsou*cos2taz),0.)
        cecoeff(3)=acc*cmplx((usou-0.5*vsou*xksou
     #           -0.5*xksou*vsou*cos2taz),0.)
        cecoeff(4)=acc*cmplx(0.,-(vpsou-vsou+usou*xksou)*sintaz)
        cecoeff(5)=acc*cmplx(0.,(vpsou-vsou+usou*xksou)*costaz)
        cecoeff(6)=acc*cmplx(xksou*vsou*sin2taz,0.)
c
c---- radial derivative
c
        usourad=(-usou+upsou)
        upsourad=(uppsou)
        vsourad=(-vsou+vpsou)
        vpsourad=(vppsou)
c
        cecoeffrad(1)=acc*cmplx(upsourad,0.)
        cecoeffrad(2)=acc*cmplx((usourad-0.5*vsourad*xksou
     #           +0.5*xksou*vsourad*cos2taz),0.)
        cecoeffrad(3)=acc*cmplx((usourad-0.5*vsourad*xksou
     #           -0.5*xksou*vsourad*cos2taz),0.)
        cecoeffrad(4)=acc*cmplx(0.,-(vpsourad-vsourad+usourad*xksou)*sintaz)
        cecoeffrad(5)=acc*cmplx(0.,(vpsourad-vsourad+usourad*xksou)*costaz)
        cecoeffrad(6)=acc*cmplx(xksou*vsourad*sin2taz,0.)
c
c---- epicentral derivative
c
        cecoefftran(1)=cmplx(0.,0.)
        cecoefftran(2)=tanfac*acc*cmplx(xksou*vsou*sin2taz,0.)
        cecoefftran(3)=tanfac*acc*cmplx(-xksou*vsou*sin2taz,0.)
        cecoefftran(4)=tanfac*acc*cmplx(0.,(vpsou-vsou+usou*xksou)*costaz)
        cecoefftran(5)=tanfac*acc*cmplx(0.,(vpsou-vsou+usou*xksou)*sintaz)
        cecoefftran(6)=tanfac*acc*cmplx(-2.*xksou*vsou*cos2taz,0.)
        if(prlo.gt.0.00001) then
          do i=1,6
            cecoefftran(i)=cecoefftran(i)+cmplx(sinfac*prlod/prlo,0.)*cecoeff(i)
          enddo
        endif
c
      else
        stop 'wrong mode type in sourcerec'
      endif
c
      clongfac=cmplx(0.,-xksou)
c      clongfac=cmplx(-1./tan(delta),-xksou)
c
      cecoefflong(1)=cecoeff(1)*clongfac
      cecoefflong(2)=cecoeff(2)*clongfac
      cecoefflong(3)=cecoeff(3)*clongfac
      cecoefflong(4)=cecoeff(4)*clongfac
      cecoefflong(5)=cecoeff(5)*clongfac
      cecoefflong(6)=cecoeff(6)*clongfac
c
      cecoeffns(1)=cmplx(cos(takeoff),0.)*cecoefflong(1)+cmplx(-sin(takeoff),0.)*cecoefftran(1)
      cecoeffns(2)=cmplx(cos(takeoff),0.)*cecoefflong(2)+cmplx(-sin(takeoff),0.)*cecoefftran(2)
      cecoeffns(3)=cmplx(cos(takeoff),0.)*cecoefflong(3)+cmplx(-sin(takeoff),0.)*cecoefftran(3)
      cecoeffns(4)=cmplx(cos(takeoff),0.)*cecoefflong(4)+cmplx(-sin(takeoff),0.)*cecoefftran(4)
      cecoeffns(5)=cmplx(cos(takeoff),0.)*cecoefflong(5)+cmplx(-sin(takeoff),0.)*cecoefftran(5)
      cecoeffns(6)=cmplx(cos(takeoff),0.)*cecoefflong(6)+cmplx(-sin(takeoff),0.)*cecoefftran(6)
c
      cecoeffew(1)=cmplx(-sin(takeoff),0.)*cecoefflong(1)+cmplx(-cos(takeoff),0.)*cecoefftran(1)
      cecoeffew(2)=cmplx(-sin(takeoff),0.)*cecoefflong(2)+cmplx(-cos(takeoff),0.)*cecoefftran(2)
      cecoeffew(3)=cmplx(-sin(takeoff),0.)*cecoefflong(3)+cmplx(-cos(takeoff),0.)*cecoefftran(3)
      cecoeffew(4)=cmplx(-sin(takeoff),0.)*cecoefflong(4)+cmplx(-cos(takeoff),0.)*cecoefftran(4)
      cecoeffew(5)=cmplx(-sin(takeoff),0.)*cecoefflong(5)+cmplx(-cos(takeoff),0.)*cecoefftran(5)
      cecoeffew(6)=cmplx(-sin(takeoff),0.)*cecoefflong(6)+cmplx(-cos(takeoff),0.)*cecoefftran(6)
c
      return
      end
