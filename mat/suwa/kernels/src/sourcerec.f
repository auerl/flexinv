      subroutine sourcerec(omega,lora,stazi,stdip,takeoff,backazi,
     #           ospli0,dospli,nspl,splefs,splefr,maxspl,maxefs,cecoeff)
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
      complex acc
      parameter (twopi=6.2831853072)
c
      stdipdeg=stdip*360./twopi
      stazideg=stazi*360./twopi
c
c---- calculate projection onto seismometer axis
c
      if(abs(stdipdeg).gt.1.) then
        prtr=0.
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
      endif
c
c---- we have some conversion from the angles defined in getpaparm
c
      costaz=cos(0.25*twopi-takeoff)
      sintaz=sin(0.25*twopi-takeoff)
      cos2taz=cos(2.*(0.25*twopi-takeoff))
      sin2taz=sin(2.*(0.25*twopi-takeoff))
c
c---- toroidal modes
c
      if(lora.eq.1) then
        xksou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,3,lora))
        wsou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,7,lora))
        wpsou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,8,lora))
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
c---- spheroidal modes
c
      else if(lora.eq.2) then
        xksou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,3,lora))
        usou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,5,lora))
        upsou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,6,lora))
        vsou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,7,lora))
        vpsou=ecbspl(omega,ospli0,dospli,nspl,splefs(1,8,lora))
        uddrec=ecbspl(omega,ospli0,dospli,nspl,splefr(1,9,lora))
        vddrec=ecbspl(omega,ospli0,dospli,nspl,splefr(1,10,lora))
c
        if(ifve.eq.1) then
          acc=cmplx(uddrec,0.)/cmplx(0.,-omega**3)
        else
          acc=cmplx(0.,vddrec*prlo)/cmplx(0.,-omega**3)
        endif
        cecoeff(1)=acc*cmplx(upsou,0.)
        cecoeff(2)=acc*cmplx((usou-0.5*vsou*xksou
     #           +0.5*xksou*vsou*cos2taz),0.)
        cecoeff(3)=acc*cmplx((usou-0.5*vsou*xksou
     #           -0.5*xksou*vsou*cos2taz),0.)
        cecoeff(4)=acc*cmplx(0.,-(vpsou-vsou+usou*xksou)*sintaz)
        cecoeff(5)=acc*cmplx(0.,(vpsou-vsou+usou*xksou)*costaz)
        cecoeff(6)=acc*cmplx(xksou*vsou*sin2taz,0.)
      else
        stop 'wrong mode type in sourcerec'
      endif
      return
      end
