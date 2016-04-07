      subroutine swraysum(sourcemodes,receivermodes,dispersion,
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
      logical spherical
c
c
      parameter (maxspl=40)
      parameter (maxefs=12)
      dimension splefs(maxspl,maxefs,2)
      dimension splefr(maxspl,maxefs,2)
c
c---- more splines
c
      dimension phavelspl(maxspl)
      dimension souphaspl(maxspl)
      dimension ampspl(maxspl)
c
c
      parameter (maxpts=100000)
      dimension array(maxpts)
      dimension buffer(maxpts)
      complex carray(maxpts/2)
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
        if(.not.spherical) then
          lu=13
          call getemodl(lu,dispersion,ierror)
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
      do i=1,nps
        buffer(i)=0.
      enddo
c
c---- loop on orbits and loop on Love and Rayleigh waves
c
      do iorbit=1,norbits
        do lora=1,2
          if(isvert.eq.1.and.lora.eq.1) then
          else
            call splseismogram(eplagc,eplo,stlagc,stlo,lora,xm,
     #              stdip,stazi,iorbit,
     #              ospli0,dospli,nspl,splefs,splefr,
     #              maxspl,maxefs,
     #              spherical,
     #              distance,phavelspl,souphaspl,ampspl)
c
            do i=2,numbsamples/2
              omega=float(i-1)*domega
              if(i.gt.1.and.omega.lt.ommax) then
                    pvel=premgephvelo(lora,omega)
                    pvel=pvel+ecbspl(omega,ospli0,dospli,nspl,phavelspl)
                    ampl=ecbspl(omega,ospli0,dospli,nspl,ampspl)
                    phase=ecbspl(omega,ospli0,dospli,nspl,souphaspl)
                    phase=phase+omega*distance/pvel
                    carray(i)=cmplx(ampl*cos(-phase),ampl*sin(-phase))
                    carray(i)=carray(i)*cmplx(cos(t0*omega),sin(t0*omega))
              else
                    carray(i)=cmplx(0.,0.)
              endif
            enddo
c
            call rfour(array,n2logneeded-1,-1)
            do i=1,nps
                  buffer(i)=buffer(i)+fftfact*array(i)
            enddo
          endif
        enddo
      enddo
      return
      end
