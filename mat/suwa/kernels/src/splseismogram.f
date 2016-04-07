      subroutine splseismogram(eplagc,eplo,stlagc,stlo,lora,xm,
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
      dimension amplitude(maxnew)
      dimension addphase(maxnew)
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
c
c---- for the output
c
      dimension phavelspl(maxspl)
      dimension souphaspl(maxspl)
      dimension ampspl(maxspl)
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
        qprem=premgesmallq(lora,omega)
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
        call sourcerec(omega,lora,stazi,stdip,takeoff,backazi,
     #       ospli0,dospli,nspl,splefs,splefr,maxspl,maxefs,ecoeff)
c
        realsource=0.
        aimasource=0.
        do j=1,6
          realsource=realsource+xm(j)*1.e-30*ecoeff(1,j)
          aimasource=aimasource+xm(j)*1.e-30*ecoeff(2,j)
        enddo
        sophase=-atan2(aimasource,realsource)
        soamp=sqrt(aimasource**2+realsource**2)
        if(iprtlv.gt.2) then
          write(6,"(i5,f10.4,e12.3)") ipt,sophase,soamp
        endif
c
c---- sum source, propagation, and receiver phase and amplitude
c
        addphase(ipt)=sophase+phasestatic
        phavel(ipt)=(omega*distance/phasecum)-premgephvelo(lora,omega)
        amplitude(ipt)=soamp*attenuation*spreading*0.95179*0.01*6371.
c
c---- another factor of r0 -- not clear where it comes from
c
        amplitude(ipt)=amplitude(ipt)*6371.
        omnew(ipt)=omega
      enddo
      
      call findcbsplmat(omnew,newpt,ospli0,dospli,nspl,amat,dbata,work)
c
      call findcbsplcoef(phavel,newpt,amat,dbata,nspl,phavelspl,work)
      call findcbsplcoef(addphase,newpt,amat,dbata,nspl,souphaspl,work)
      call findcbsplcoef(amplitude,newpt,amat,dbata,nspl,ampspl,work)
c
      if(iprtlv.gt.1) then
        write(6,"('splined the phase and amp:',2i5)") nspl,newpt
      endif
      return
      end
