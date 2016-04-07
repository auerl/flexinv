      subroutine splseismogramG(eplagc,eplo,stlagc,stlo,lora,xm,
     #              stdip,stazi,iorbit,
     #              ospli0,dospli,nspl,splefs,splefr,
     #              maxspl,maxefs,
     #              spherical,
     #              distance,phavelspl,souphaspl,ampspl)
c
      parameter (twopi=6.2831853072)
      common /plevel/iprtlv
c
c---- for phase velocities
c
      parameter (maxsplpv=30)
      dimension omarr(maxsplpv)
      dimension delcgc(maxsplpv)
      dimension delc(maxsplpv)
      dimension delcgcspl(3,maxsplpv)
      dimension delcspl(3,maxsplpv)
c
c---- for the second splining
c
      parameter (maxnew=100)
      dimension amplitude(maxnew,10)
      dimension addphase(maxnew,10)
      dimension phavel(maxnew)
      dimension omnew(maxnew)
      dimension amat(maxnew*maxnew)
      dimension work(maxnew*2)
      real*8 dbata(maxnew,maxnew)
c
c---- for eigenfunctions
c
      dimension splefs(maxspl,maxefs,2)
      dimension splefr(maxspl,maxefs,2)
c
c---- for the source-receiver amplitude
c
      dimension xm(6)
      dimension ecoeff(2,10)
      dimension ecoeffrad(2,10)
      dimension ecoeffns(2,10)
      dimension ecoeffew(2,10)
c
c---- for the output
c
      dimension phavelspl(maxspl)
      dimension souphaspl(maxspl,10)
      dimension ampspl(maxspl,10)
c
      logical spherical
c
c---- get the path angles
c
      call delazgc(eplagc,eplo,stlagc,stlo,delta0deg,azatepdeg,azatstdeg)
      delta0=delta0deg*twopi/360.
      azatep=azatepdeg*twopi/360.
      azatst=azatstdeg*twopi/360.
c
      if(.not.spherical) then
        call getdelcspl(eplagc,eplo,stlagc,stlo,lora,omarr,
     #       delc,delcgc,delcspl,delcgcspl,maxsplpv,nsplpv,ierror)
      if(iprtlv.gt.1) then
        write(6,"('phase velocity spline points:',i5)") nsplpv
        write(6,"('after getdelcspl ierror:',i5)") ierror
        write(6,"('omarr',8f10.5)") (omarr(i),i=1,nsplpv)
        write(6,"('delc',8f10.5)") (delc(i),i=1,nsplpv)
      endif
      endif
c
      ifoddorbit=mod(iorbit,2)
c      
      if(ifoddorbit.eq.1) then
        distance=delta0*6371.+float(iorbit/2)*twopi*6371.
        takeoffgc=azatep
        backazigc=azatst
        weighthat=float(iorbit/2)*twopi/(delta0+float(iorbit/2)*twopi)
        weighttilde=delta0/(delta0+float(iorbit/2)*twopi)
      else
        distance=-delta0*6371.+float(iorbit/2)*twopi*6371.
        takeoffgc=azatep-twopi/2.
        if(takeoffgc.lt.0.) takeoffgc=takeoffgc+twopi
        backazigc=azatst-twopi/2.
        if(backazigc.lt.0.) backazigc=backazigc+twopi
        weighthat=float(iorbit/2)*twopi/(float(iorbit/2)*twopi-delta0)
        weighttilde=-delta0/(float(iorbit/2)*twopi-delta0)
      endif
c
c---- initialize the angles at source and receiver
c
      takeoff=takeoffgc
      backazi=backazigc
c
c---- loop on frequencies to calculate new spline point values
c
      phasestatic=twopi/8.-float(iorbit-1)*twopi/4.
c
c---- sample existing spline functions at twice the frequency
c---- calculate the ATA^-1 matrix for a given number of spline points
c
      newpt=2*nspl-1
      do ipt=1,newpt
        omega=ospli0+float(ipt-1)*dospli*0.5
        phasecum=0.
        pvel=premgephvelo(lora,omega)
        if(.not.spherical) then
            pert=rsple(1,nsplpv,omarr,delc,delcspl,omega)
            pertgc=rsple(1,nsplpv,omarr,delcgc,delcgcspl,omega)
      if(iprtlv.gt.1) then
            write(6,"('pert,pertgc,weighttilde,weighthat',4g14.4)") 
     #         pert,pertgc,weighthat,weighttilde
      endif
            pvel=pvel*(1.+0.01*weighthat*pertgc+0.01*weighttilde*pert)
        endif
c
c---- attenuation and spreading depend on two additional
c---- path parameters - the group velocity and the attenuation
c---- here both of these are taken as the PREM values
c
        gvelprem=premgegrvelo(lora,omega)
c
c---- modified 9/27/1999 to use QL6 Q values
c
        qprem=dureksmallq(lora,omega)
ccc        qprem=premgesmallq(lora,omega)
c
c---- the wavenumber at the receiver enters into the spreading
c
        xkrec=ecbspl(omega,ospli0,dospli,nspl,splefr(1,3,lora))
c
        phasecum=phasecum+omega*distance/pvel
c
        attenuation=exp(-omega*distance*qprem/(2.*gvelprem))
c
        spreading=1./sqrt(4.*twopi*xkrec*sin(delta0))
c
c
c---- calculate the source and receiver phase and amplitude
c
        call sourcerecG(omega,lora,stazi,stdip,delta0,takeoff,backazi,
     #       ospli0,dospli,nspl,splefs,splefr,maxspl,maxefs,
     #       ecoeff,ecoeffrad,ecoeffns,ecoeffew)
c
        if(iprtlv.gt.0) then
          write(6,"(6e12.3)") omega,stazi,stdip,delta0,takeoff,backazi
          write(6,"(3i5)") lora,nspl,maxspl
          write(6,"(10e12.3)") ecoeff
          write(6,"(10e12.3)") ecoeffns
        endif
c
        do ifact=1,6
          sophase=-atan2(ecoeff(2,ifact),ecoeff(1,ifact))
          soamp=sqrt(ecoeff(2,ifact)**2+ecoeff(1,ifact)**2)
          addphase(ipt,ifact)=sophase+phasestatic
          amplitude(ipt,ifact)=soamp*attenuation*spreading*0.95179*0.01*6371.**2
        enddo
        realsource=0.
        aimasource=0.
        realsourcerad=0.
        aimasourcerad=0.
        realsourcens=0.
        aimasourcens=0.
        realsourceew=0.
        aimasourceew=0.
        do j=1,6
          realsource=realsource+xm(j)*1.e-30*ecoeff(1,j)
          aimasource=aimasource+xm(j)*1.e-30*ecoeff(2,j)
          realsourcerad=realsourcerad+xm(j)*1.e-30*ecoeffrad(1,j)
          aimasourcerad=aimasourcerad+xm(j)*1.e-30*ecoeffrad(2,j)
          realsourcens=realsourcens+xm(j)*1.e-30*ecoeffns(1,j)
          aimasourcens=aimasourcens+xm(j)*1.e-30*ecoeffns(2,j)
          realsourceew=realsourceew+xm(j)*1.e-30*ecoeffew(1,j)
          aimasourceew=aimasourceew+xm(j)*1.e-30*ecoeffew(2,j)
        enddo
        sophase=-atan2(aimasource,realsource)
        soamp=sqrt(aimasource**2+realsource**2)
        if(iprtlv.gt.2) then
          write(6,"(i5,f10.4,e12.3)") ipt,sophase,soamp
        endif
c
        sophaserad=-atan2(aimasourcerad,realsourcerad)
        soamprad=sqrt(aimasourcerad**2+realsourcerad**2)
c
        sophasens=-atan2(aimasourcens,realsourcens)
        soampns=sqrt(aimasourcens**2+realsourcens**2)
c
        sophaseew=-atan2(aimasourceew,realsourceew)
        soampew=sqrt(aimasourceew**2+realsourceew**2)
c
        addphase(ipt,7)=sophase+phasestatic
        amplitude(ipt,7)=soamp*attenuation*spreading*0.95179*0.01*6371.**2
        addphase(ipt,8)=sophaserad+phasestatic
        amplitude(ipt,8)=soamprad*attenuation*spreading*0.95179*0.01*6371.**2
        addphase(ipt,9)=sophasens+phasestatic
        amplitude(ipt,9)=soampns*attenuation*spreading*0.95179*0.01*6371.**2
        addphase(ipt,10)=sophaseew+phasestatic
        amplitude(ipt,10)=soampew*attenuation*spreading*0.95179*0.01*6371.**2
c
        phavel(ipt)=(omega*distance/phasecum)-premgephvelo(lora,omega)
        omnew(ipt)=omega
      enddo
      
      call findcbsplmat(omnew,newpt,ospli0,dospli,nspl,amat,dbata,work)
      call findcbsplcoef(phavel,newpt,amat,dbata,nspl,phavelspl,work)
c
      do j=1,10
        do ipt=2,newpt
          diff=addphase(ipt,j)-addphase(ipt-1,j)
          if(2.*diff.gt.twopi) then
            addphase(ipt,j)=addphase(ipt,j)-twopi
          else if(2.*diff.lt.-twopi) then
            addphase(ipt,j)=addphase(ipt,j)+twopi
          endif
        enddo
        call findcbsplcoef(addphase(1,j),newpt,amat,dbata,nspl,souphaspl(1,j),work)
        call findcbsplcoef(amplitude(1,j),newpt,amat,dbata,nspl,ampspl(1,j),work)
        if(iprtlv.gt.0.and.j.eq.9) then
          write(6,"('amp',10e12.3)") (ampspl(i,j),souphaspl(i,j),i=1,nspl)
        endif
      enddo
c
      if(iprtlv.gt.1) then
        write(6,"('splined the phase and amp:',2i5)") nspl,newpt
      endif
      return
      end
