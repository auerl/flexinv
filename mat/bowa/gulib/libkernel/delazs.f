      subroutine delazs(eplat,eplong,stlat,stlong,delta,azep,azst)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      data hpi,twopi,rad,reprad/1.5707963268
     1   ,6.2831853,.017453293,57.2957795/
      arcos(x)=atan2(sqrt(1.-x*x),x)
      el=eplat*rad
      el=hpi-el
      stl=stlat*rad
      stl=hpi-stl
      elon=eplong*rad
      slon=stlong*rad
      as=cos(stl)
      bs=sin(stl)
      cs=cos(slon)
      ds=sin(slon)
      a=cos(el)
      b=sin(el)
      c=cos(elon)
      d=sin(elon)
      co=1.0
      cdel=a*as+b*bs*(c*cs+d*ds)
      if(abs(cdel).gt.1.0) cdel=sign(co,cdel)
      delt=arcos(cdel)
      delta=delt*reprad !distance
      sdel=sin(delt)
      caze=(as-a*cdel)/(sdel*b)
      if(abs(caze).gt.1.0) caze=sign(co,caze)
      aze=arcos(caze)
      if(bs.gt.0.0) cazs=(a-as*cdel)/(bs*sdel)
      if(bs.eq.0.0) cazs=sign(co,cazs)
      if(abs(cazs).gt.1.0) cazs=sign(co,cazs)
      azs=arcos(cazs)
      dif=ds*c-cs*d
      if(dif.lt.0.0) aze=twopi-aze
      azep=reprad*aze
      if(dif.gt.0.0) azs=twopi-azs
      azst=reprad*azs
      return
      end
