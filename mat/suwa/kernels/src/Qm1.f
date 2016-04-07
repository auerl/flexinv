      character*80 getunx
      character*80 string
      character*80 sysstring
c
      character*80 cmtsolution
      character*80 recordlist
      character*80 outdir
      character*80 sourcemodes
      character*80 receivermodes
      character*80 dispersion
      character*80 dbsdir
c
      logical spherical
c
      character*80 newfile
      character*80 recordfile
      logical exists
      integer system
c
      character*8 stn
      character*8 chn
      character*8 cmnt
      character*8 stn0 /'        '/
      real*4 stla0 /-99.9/
      real*4 stlo0 /-99.9/
c
      parameter (twopi=6.2831853072)
      dimension xm(6)
      real*8 dtlast
      real*8 dt0
      parameter (maxpts=100000)
      dimension array(maxpts)
      complex carray(maxpts/2)
      dimension buffer(maxpts)
      equivalence (array(1),carray(1))
      dimension resparr(1000)
      dimension zerarr(4)
      data zerarr /0.,0.,0.,0./
c
      logical fundonly /.true./
      logical acconly /.false./
      logical overonly /.false./
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
c---- spline parameters
c
      ospli0=twopi/500.
      dospli=twopi/750.

c
      call chekcl('|   -C:o:1:CMT source file [CMTSOLUTION]'//
     #            '|   -L:o:1:List with records [DATA/RECORDHEADERS]'//
     #            '|   -D:o:1:Directory for output files [SYNT]'//
     #            '|   -S:o:1:Mode catalog (source) [PREM.FUND]'//
     #            '|   -R:o:1:Mode catalog (receiver) [PREM.FUND]'//
     #            '|   -M:o:1:Phase velocity maps [DISPERSION]'//
     #            '|   -O:o:1:Minimum period to consider [25.0]'//
     #            '|   -P:o:1:Maximum period to consider [333.0]'//
     #            '|   -o:o:1:Maximum orbit [1]'//
     #            '|   -s:o:0:Spherical'//
     #            '|   -v:o:0,1:Verbosity level|')
c
c---- get option flags
c
      spherical=.false.
      string=getunx('-s',0,nbyts)
      if(nbyts.eq.0) spherical=.true.
c
      string=getunx('-o',0,nbyts)
      if(nbyts.gt.0) then
        read(string,*) norbits
        write(6,"('will add',i2,' orbits')") norbits
      else
        norbits=1
      endif
c
      call getenv('HRVDBS',dbsdir)
      ldbs=lnblnk(dbsdir)
c
c---- get the filename for moment tensor
c
      cmtsolution=getunx('-C',1,nbyts)
      inquire(file=cmtsolution,exist=exists)
      if(.not.exists) then
        write(6,"('MT file does not exist:',a)") cmtsolution(1:nbyts)
        call exit(1)
      endif
c
c---- get the filename for list of record/path information
c
      recordlist=getunx('-L',1,nbyts)
      inquire(file=recordlist,exist=exists)
      if(.not.exists) then
        write(6,"('path list does not exist:',a)") recordlist(1:nbyts)
        call exit(1)
      endif
c
c---- get the name of the output directory
c
      outdir=getunx('-D',1,loutdir)
      if(loutdir.gt.0) then
        inquire(file=outdir,exist=exists)
        if(.not.exists) then
          sysstring='mkdir -p '//outdir(1:loutdir)
          ierror=system(sysstring)
          if(ierror.ne.0) then
            write(6,"('unable to create:',a)") outdir(1:loutdir)
            call exit(1)
          endif
        endif
      endif
c
c---- get the name of the mode catalog (source)
c
      sourcemodes=getunx('-S',1,nbyts)
      if(sourcemodes(1:9).eq.'PREM.FUND') then
        sourcemodes=dbsdir(1:ldbs)//'/PREM.FUND'
      endif
      lf=lnblnk(sourcemodes)
      inquire(file=sourcemodes,exist=exists)
      if(.not.exists) then
        write(6,"('mode file does not exist:',a)") sourcemodes(1:lf)
        call exit(1)
      endif
      write(6,"('source model:      ',a)") sourcemodes(1:lf)
c
c---- get the name of the mode catalog (receiver)
c
      receivermodes=getunx('-R',1,nbyts)
      if(receivermodes(1:9).eq.'PREM.FUND') then
        receivermodes=dbsdir(1:ldbs)//'/PREM.FUND'
      endif
      lf=lnblnk(receivermodes)
      inquire(file=receivermodes,exist=exists)
      if(.not.exists) then
        write(6,"('mode file does not exist:',a)") receivermodes(1:lf)
        call exit(1)
      endif
      write(6,"('receiver model:    ',a)") receivermodes(1:lf)
c
c---- get the name of the dispersion maps
c
      if(.not.spherical) then
        dispersion=getunx('-M',1,nbyts)
        if(dispersion(1:nbyts).eq.'DISPERSION') then
          dispersion=dbsdir(1:ldbs)//'/DISPERSION'
        endif
        lf=lnblnk(dispersion)
        inquire(file=dispersion,exist=exists)
        if(.not.exists) then
          write(6,"('map file does not exist:',a)") dispersion(1:lf)
          call exit(1)
        endif
        write(6,"('phase velocities:  ',a)") dispersion(1:lf)
      endif
c
c---- minimum and maximum period for synthetics
c
      read(getunx('-O',1,nbyts),*) permin
      ommax=twopi/permin
      read(getunx('-P',1,nbyts),*) permax
      ommin=twopi/permax
      write(6,"('minimum and maximum period:',2f8.2)") permin,permax
c
c---- read the CMT information
c
      open(1,file=cmtsolution)
      call rcmt(1,iyr,imo,ida,iho,imi,fsec,
     #          tcmt,hdur,epla,eplo,dep,xm,ierr)
      jda=julday(iyr,imo,ida)
      write(6,"('time shift:',f12.4)") tcmt
      write(6,"('half duration:',f9.4)") hdur
      write(6,"('latitude:  ',f12.4)") epla
      write(6,"('longitude: ',f12.4)") eplo
      write(6,"('depth:     ',f12.4)") dep
      write(6,"('Mrr:    ',e15.6)") xm(1)
      write(6,"('Mtt:    ',e15.6)") xm(2)
      write(6,"('Mpp:    ',e15.6)") xm(3)
      write(6,"('Mrt:    ',e15.6)") xm(4)
      write(6,"('Mrp:    ',e15.6)") xm(5)
      write(6,"('Mtp:    ',e15.6)") xm(6)
      close(1)
      eplagc=atand(0.993277*tand(epla))
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
c---- initialize the dispersion 
c
      if(.not.spherical) then
        lu=13
        call getemodl(lu,dispersion,ierror)
        if(ierror.ne.0) then
          stop 'error reading the dispersion files'
        endif
      endif
c
c---- start loop on records to be synthesized
c
      open(1,file=recordlist)
c
      ios=0
      do while(ios.eq.0)
        read(1,"(a8,1x,a8,1x,f8.4,1x,f9.4,1x,f6.1,1x,f6.1,f6.1,1x,f6.1,
     #    1x,f12.4,1x,i7,1x,i4,1x,i3,1x,i2,1x,i2,1x,f6.3)",iostat=ios)
     #    stn,chn,stla,stlo,stele,stbur,stazideg,stdipdeg,spsec,nps,
     #    iyr1,ijd1,iho1,imi1,fsc1
          stazi=stazideg*twopi/360.
          stdip=stdipdeg*twopi/360.
        write(6,"(4i5,f6.1,4i5,f6.1)") iyr1,ijd1,iho1,imi1,fsc1,
     #    iyr,jda,iho,imi,fsec+tcmt
        if(ios.eq.0) then
          if(iyr1.eq.0) then
            call tadder(iyr,jda,iho,imi,fsec,
     #           iyr1,ijd1,iho1,imi1,fsc1,dble(tcmt))
            t0=0.
          else
            call tdiffer(iyr1,ijd1,iho1,imi1,fsc1,
     #           iyr,jda,iho,imi,fsec+tcmt,dt0)
            t0=sngl(dt0)
          endif
          stlagc=atand(0.993277*tand(stla))
c
          if(abs(stdipdeg+90.0).lt.1.) then
            isvert=1
          else
            isvert=0
          endif
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
          df=spsec/float(numbsamples)
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
c---- get the receiver structure
c
          if(stn.ne.stn0.or.stla.ne.stla0.or.stlo.ne.stlo0) then
            call openfl(1,receivermodes,1,0,0,ierr,-1)
            call getgemodes(1,2,fundonly,overonly,acconly,dep)
            call closfl(1,ierr)
            do lora=1,2
              call spleigen(2,lora,ommin,ommax,ospli0,dospli,
     #             nspl,splefr,maxspl,maxefs)
            enddo
            stn0=stn
            stla0=stla
            stlo0=stlo
          endif
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
c
c---- the seismogram is done - create a new output file 
c
          call mnam(stn,chn,'record',outdir,loutdir,newfile)
          if(loutdir.gt.0) then
            recordfile=outdir(1:loutdir)//'/'//newfile
          else
            recordfile=newfile
          endif
          write(6,"('new file:',a)") recordfile(1:lnblnk(recordfile))
          open(13,file=recordfile)
          close(13)
          call openfl(3,recordfile,4,0,0,istat,-1)
c
c---- add a signature for the seismogram header and write it out
c
          cmnt(1:4)='RAY1'
          if(.not.spherical) then
            cmnt(5:6)='-H'
          else
            cmnt(5:6)='-S'
          endif
          dtlast=dfloat(nps-1)/dble(spsec)
          call wrec(3,stn,chn,dble(stla),dble(stlo),stele,stbur,stazideg,
     #        stdipdeg,spsec,nps,iyr1,ijd1,iho1,imi1,fsc1,dtlast,cmnt)
c
c---- write out the response (acceleration response)
c
          call packresponse(resparr,resparr,resparr,lres,1.,0,zerarr,0,
     #        zerarr)
          call wres(3,resparr,lres)
c
c---- write out a comment
c
          call wcom(3,'Qm1',3)
c
c---- write out the seismogram
c
          call wsei(3,buffer,nps)
          call closfl(3,ierror)
        endif
      enddo
      end     
      
