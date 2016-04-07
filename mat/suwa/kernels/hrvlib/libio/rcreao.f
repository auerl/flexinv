      subroutine rcreao(path,ichan)
      character*1 null
      parameter (null=char(0))
      character*(*) path
      ichan=-1
      itry=0
      ifail=0
      ilen=istlen(path)
c read-write
      isunop=2
   30 ip2=ilen
      if(path(ip2:ip2).eq.'/'.and.ifail.ne.0) return
      call copen(path(1:ip2)//null,ichan,isunop,iopner,2,z'1a4')
      if(iopner.eq.0) then
c       call cclose(ichan,ires,ier)
        return
      else if(ifail.ne.0) then
        call check('copen in recreat failed')
      else if(iopner.eq.2) then
        itry=itry+1
        ifail=1
      else





        call check('copen in rcreat')
      endif

   10 continue
      do while(ip2.gt.0.and.path(ip2:ip2).ne.'/') 
        ip2=ip2-1
      enddo

   20 call cmkdir(path(1:ip2-1)//null,z'1ed',ires,ierrno)
      if(ierrno.eq.2) then
        itry=1+itry
        ip2=ip2-1
        goto 10
      else if(ierrno.eq.0) then
        itry=itry-1
        if(itry.gt.0) then
          ip2=1+ip2
          do while(path(ip2:ip2).ne.'/')
            ip2=1+ip2
          enddo

          goto 20
        else 




          goto 30
        endif

      else





        call check('rcreat')
      endif

      return
      end
