      character*80 getunx
      character*80 modecat
      parameter (twopi=6.2831853072)
      dimension array (2000,4)
c
      call chekcl('|   -M:o:1:Mode catalog [SPRM1.BIN ]'//
     #            '|   -v:o:0,1:Verbosity level|')
c
c---- get the name of the mode catalog
c
      write(6,"('type 1 for phase velocity')")
      write(6,"('type 2 for group velocity')")
      write(6,"('type 3 for small q')")
      read(5,"(i1)") isel
      
      modecat=getunx('-M',1,lmodecat)
      open(1,file=modecat(1:lmodecat)//'.love')
      ios=0
      ipl=0
      do while(ios.eq.0)
        read(1,*,iostat=ios) f1,f2,i1,i2,i3,f3,f4,f5,f6
        if(ios.eq.0) then
          ipl=ipl+1
          array(ipl,1)=f3
          if(isel.eq.1) then
            array(ipl,2)=f4
          else if(isel.eq.2) then
            array(ipl,2)=f5
          else if(isel.eq.3) then
            array(ipl,2)=f6
          endif
        endif
      enddo
      close(1)
      open(1,file=modecat(1:lmodecat)//'.rayl')
      ios=0
      ipr=0
      do while(ios.eq.0)
        read(1,*,iostat=ios) f1,f2,i1,i2,i3,f3,f4,f5,f6
        if(ios.eq.0) then
          ipr=ipr+1
          array(ipr,3)=f3
          if(isel.eq.1) then
            array(ipr,4)=f4
          else if(isel.eq.2) then
            array(ipr,4)=f5
          else if(isel.eq.3) then
            array(ipr,4)=f6
          endif
        endif
      enddo
      ntspl=ipl
      nsspl=ipr
c
c---- now calculate the splines
c
      write(6,"('there are',i4,' toroidal modes to spline')") ntspl
      write(6,"('there are',i4,' spheroidal modes to spline')") nsspl
      write(6,"('      dimension tomega(',i4,')')") ntspl
      write(6,"('      dimension tpvel(',i4,')')") ntspl
      write(6,"('      dimension somega(',i4,')')") nsspl
      write(6,"('      dimension spvel(',i4,')')") nsspl
      write(6,"('      data ntspl,nsspl /',i4,',',i4,'/')") ntspl,nsspl
      write(6,"('c')") 
      do j=1,4
        if(j.eq.1.or.j.eq.2) then
          nmodes=ntspl
        else
          nmodes=nsspl
        endif
        if(j.eq.1) write(6,"('      data tomega /')")
        if(j.eq.2) write(6,"('      data tpvel /')")
        if(j.eq.3) write(6,"('      data somega /')")
        if(j.eq.4) write(6,"('      data spvel /')")
c
        do i=1,nmodes
          if(i.eq.nmodes) then
            if(mod(i,6).eq.1) then
              write(6,"('     # ',f10.8,'/')") array(i,j)
            else
              write(6,"(f10.8,'/')") array(i,j)
            endif
          else 
            if(mod(i,6).eq.1) then
              write(6,"('     # ',f10.8,',',$)") array(i,j)
            else if(mod(i,6).eq.0) then
              write(6,"(f10.8,',')") array(i,j)
            else
              write(6,"(f10.8,',',$)") array(i,j)
            endif
          endif
        enddo
        write(6,"('c')") 
      enddo
      end
