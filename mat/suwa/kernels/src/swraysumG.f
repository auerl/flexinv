      subroutine swraysumG(sourcemodes,receivermodes,dispersion,
     #         spherical,norbits,
     #         epla,eplo,dep,xm,ityphdur,hdur,stla,stlo,stazideg,stdipdeg,
     #         nps,dt,t0,ommin,ommax,buffer,itypkern)
c
      character*80 sourcemodes
      character*80 receivermodes
      character*80 dispersion

      parameter (twopi=6.2831853072)
      logical fundonly /.true./
      logical acconly /.false./
      logical overonly /.false./
      logical first /.true./
      logical spherical
c
      common /plevel/ iprtlv
c
      parameter (maxspl=50)
      parameter (maxefs=12)
      dimension splefs(maxspl,maxefs,2)
      dimension splefr(maxspl,maxefs,2)
c
c---- more splines
c
      dimension phavelspl(maxspl)
      dimension souphaspl(maxspl,10)
      dimension ampspl(maxspl,10)
c
c
      parameter (maxpts=100000)
      dimension array(maxpts,11)
      dimension buffer(1)
      complex carray(maxpts/2,11)
      equivalence (array(1),carray(1))
c
      logical spherical1 /.false./
      data depth /1000./
      data omax /0./
      save spherical1, omax, depth
c
c---- spline parameters
c
      ospli0=twopi/500.
      dospli=twopi/750.
c
c---- check if it is necessary to read modes again
c
      if(depth.eq.dep.and.omax.ge.ommax.and.
     #   spherical.eq.spherical1) then
      else
c
        depth=dep
        spherical1=spherical
        omax=ommax
c
c---- get the eigenfunctions at the source location
c
        call openfl(1,sourcemodes,1,0,0,ierr,-1)
        if(iprtlv.gt.0) then
        write(6,"('getting modes at source depth',f7.2)") dep
        endif
        call getgemodes(1,1,fundonly,overonly,acconly,dep)
        call closfl(1,ierr)
        do lora=1,2
          call spleigen(1,lora,ommin,ommax,ospli0,dospli,
     #              nspl,splefs,maxspl,maxefs)
        enddo
c
c---- get the receiver structure
c---- this could be modified to change for different receivers
c
        call openfl(1,receivermodes,1,0,0,ierr,-1)
        call getgemodes(1,2,fundonly,overonly,acconly,dep)
        call closfl(1,ierr)
        do lora=1,2
          call spleigen(2,lora,ommin,ommax,ospli0,dospli,
     #         nspl,splefr,maxspl,maxefs)
        enddo
c
c---- initialize the dispersion 
c
        if(.not.spherical.and.first) then
          lu=13
          call getemodl(lu,dispersion,ierror)
          first=.false.
          if(ierror.ne.0) then
            stop 'error reading the dispersion files'
          endif
        endif
      endif
c
      eplagc=atand(0.993277*tand(epla))
      stlagc=atand(0.993277*tand(stla))
c
      if(abs(stdipdeg+90.0).lt.1.) then
        isvert=1
      else if(abs(stdipdeg-90.0).lt.1.) then
        isvert=1
      else
        isvert=0
      endif
      stazi=stazideg*twopi/360.
      stdip=stdipdeg*twopi/360.
c
c---- determine the length of the output time series and some factors
c
      numbsamples=nps
      n2logneeded=1
      do while(2**n2logneeded.lt.numbsamples)
        n2logneeded=n2logneeded+1
      enddo
      numbsamples=2**n2logneeded
      if(numbsamples.gt.maxpts) then
            stop 'too many points'
      endif
      df=1./(dt*float(numbsamples))
      domega=twopi*df
      fftfact=2.*df
c
c---- t0 is time of first sample w.r.t. CMT origin time
c---- the program returns ground acceleration 
c
      do i=1,nps*11
        buffer(i)=0.
      enddo
      nkerns=10
c
c---- loop on orbits and loop on Love and Rayleigh waves
c
c      do iorbit=1,norbits
      do iorbit=1,norbits
        do lora=1,2
          if(isvert.eq.1.and.lora.eq.1) then
          else
            call splseismogramG(eplagc,eplo,stlagc,stlo,lora,xm,
     #              stdip,stazi,iorbit,
     #              ospli0,dospli,nspl,splefs,splefr,
     #              maxspl,maxefs,
     #              spherical,
     #              distance,phavelspl,souphaspl,ampspl)
c
            do i=2,1+numbsamples/2
              omega=float(i-1)*domega
c
              if(hdur.gt.0.)then
                if(ityphdur.eq.2) then
                  factor=2.*(1.-cos(omega*hdur))/((hdur*omega)**2)
                else if(ityphdur.eq.1) then
                  factor=sin(omega*hdur)/(omega*hdur)
                else
                  stop 'ityphdur not given'
                endif
              endif
c
              if(i.gt.1.and.omega.lt.ommax) then
                    pvel=premgephvelo(lora,omega)
                    pvel=pvel+ecbspl(omega,ospli0,dospli,nspl,phavelspl)
                    do ikern=1,nkerns
                      iindex=ikern+1
                      ampl=ecbspl(omega,ospli0,dospli,nspl,ampspl(1,ikern))
                      phase=ecbspl(omega,ospli0,dospli,nspl,souphaspl(1,ikern))
                      if(iprtlv.gt.0.and.ikern.eq.9) then
c                        write(6,"('amp,phase',2e12.3)") ampl,phase
                      endif 
                      phase=phase+omega*distance/pvel
                      carray(i,iindex)=cmplx(ampl*cos(-phase),ampl*sin(-phase))
                      carray(i,iindex)=carray(i,iindex)*cmplx(cos(t0*omega),sin(t0*omega))
                      carray(i,iindex)=carray(i,iindex)*cmplx(factor,0.)
                      if(ikern.eq.7) then
                        carray(i,1)=carray(i,iindex)
                        carray(i,iindex)=carray(i,iindex)*cmplx(0.,-omega)
                      endif
                      if(ikern.eq.8) then
                        carray(i,iindex)=-carray(i,iindex)
                      endif
                    enddo
              else
                    do iindex=1,11
                      carray(i,iindex)=cmplx(0.,0.)
                    enddo
              endif
            enddo
            do iindex=1,11
                carray(1,iindex)=cmplx(0.,0.)
            enddo
c
            do iindex=1,11
              call rfour(array(1,iindex),n2logneeded-1,-1)
              do i=1,nps
                  buffer(i+(iindex-1)*nps)=buffer(i+(iindex-1)*nps)+fftfact*array(i,iindex)
                      if(iprtlv.gt.0.and.iindex.eq.10.and.i.lt.30) then
                        write(6,"(10e13.3)") array(i,iindex),buffer(i+(iindex-1)*nps)
                      endif
              enddo
            enddo
          endif
        enddo
      enddo
      return
      end
