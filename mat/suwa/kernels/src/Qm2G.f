      character*80 getunx
      character*80 string
      character*80 sysstring
      common /plevel/ iprtlv
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
c
      parameter (twopi=6.2831853072)
      dimension xm(6)
      real*8 dtlast
      real*8 dt0
c
      parameter (maxpts=1000000)
      dimension buffer(maxpts)
      dimension resparr(1000)
      dimension zerarr(4)
      data zerarr /0.,0.,0.,0./

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
      string=getunx('-v',0,nbyts)
      if(nbyts.gt.0) then
        read(string,*) iprtlv
      else
        iprtlv=0
      endif
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
c
c---- call the surface wave ray summation subroutine
c
          dt=1./spsec
          ityphdur=1
          itypkern=0
          call swraysumG(sourcemodes,receivermodes,dispersion,
     #         spherical,norbits,
     #         epla,eplo,dep,xm,ityphdur,hdur,stla,stlo,stazideg,stdipdeg,
     #         nps,dt,t0,ommin,ommax,buffer,itypkern)
          
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
          cmnt(1:4)='RAY2'
          if(.not.spherical) then
            cmnt(5:6)='-H'
          else
            cmnt(5:6)='-S'
          endif
          dtlast=dfloat(nps-1)/dble(spsec)
          call wrec(3,stn,chn,dble(stla),dble(stlo),stele,stbur,stazideg,
     #        stdipdeg,spsec,nps,iyr1,ijd1,iho1,imi1,fsc1,dtlast,cmnt)
c
c---- write out the response (displacement response)
c
          call packresponse(resparr,resparr,resparr,lres,1.,0,zerarr,0,
     #        zerarr)
          call wres(3,resparr,lres)
c
c---- write out a comment
c
          call wcom(3,'Qm2',3)
c
c---- write out the seismogram
c
          call wsei(3,buffer,nps)
          call closfl(3,ierror)
        endif
      enddo
      end     
      
