      subroutine oassend(command,response,lr,iwt)
      character*(*) command,response
      parameter (maxoas=200)
      character*1 cr,nl
      parameter (cr=char(13))
      parameter (nl=char(10))
      character*(maxoas) msg
      character*(maxoas) temp
      logical condition
      include "oasstatus.h"
c
      ilc=0
      ip=0
      lc=lnblnk(command)
c
      do while (ilc.lt.lc)
        ilc=ilc+1
        if(command(ilc:ilc).eq.'^') then
          ip=ip+1
          ilc=ilc+1
          temp(ip:ip)=char(and(z'0f',ichar(command(ilc:ilc))))
        else
          ip=ip+1
          temp(ip:ip)=command(ilc:ilc)
        endif
      enddo

      ip=1+ip
      if (ip.gt.79) then
         write(*,*) '**** WARNING: OASSEND: ip = ',ip
         write(*,'(a)') temp(1:ip)
      endif
      temp(ip:ip)=cr
      if(iwt.ge.0) then
        call cnoblock(ischan,1,ires,ierrno)
        if(ierrno.ne.0) call check('cnoblock in oassend')
      endif
      if(ioasverbose.gt.0) then
        write(6,"(a,a)") 'sending to OAS: ',temp(1:ip-1)
      endif
c
ccccccc trying to fix reading 3/22/93
      call cnoblock(ischan,1,ires,ierrno)
      iwrote=0
c
   22 continue
      call cwrite(ischan,temp,ip,ires,ierrno)
      iwrote=iwrote+1
      if(ierrno.ne.0) call check('cwrite in oassend')
cccccc removed the following line to not intefere with reading? 3/22/93
c      call cflush(ischan,1,ires,ierrno)
cccccc-----------------------------------------------------------------
      if(iwt.gt.0) call wait(iwt)
      lr=0

      if(ioasverbose.gt.0) then
        write(6,"(a)") 'begin reading from OAS: '
      endif
      ierrno=0
      condition=.true.
      do while(ierrno.eq.0.and.condition)
        lr=lr+1
        millis=0
   11   continue
	if(command(7:11).eq.'blank'.or.command(8:12).eq.'blank') then
	  if(millis.eq.0) write(6,"('searching for a blank area on disk')")
	  millis=10
        else
	  if(millis.gt.0.and.mod(millis,2000).eq.0) then
	    write(6,"('waiting',i3,' seconds for OAS response')") millis/1000
          endif
	endif
        if(millis.lt.30000) then
          call cread(ischan,msg(lr:lr),1,igot,ierrno)
          if(ierrno.ne.0) then
            call wait(10)
            millis=millis+10
            go to 11
          endif
        else
          write(6,"('giving up on reading from OAS:',i4)") lr
          write(6,"(a,a)") 'previously sent to OAS: ',temp(1:ip-1)
          if(iwrote.eq.1) go to 22
          condition=.false.
        endif
        condition=msg(lr:lr).ne.'>'
        if(.not.condition.and.lr.gt.3) then
          condition=(msg(lr:lr).eq.'>'.and.msg(lr-3:lr-3).eq.'<')
        endif
      enddo
      if(ioasverbose.gt.0) then
        write(6,"(a)") 'after reading from OAS: '
      endif

      if(ierrno.ne.0) lr=lr-1
      do i=1,lr
        if(msg(i:i).eq.cr) msg(i:i)=nl
        response(i:i)=msg(i:i)
      enddo

      if(msg(lr:lr).eq.nl) lr=lr-1
      if(response(1:1).eq.nl) then
        lr=lr-1
        do i=1,lr
          i1=i+1
          response(i:i)=response(i1:i1)
        enddo
      endif
cccccc trying to fix reading 3/22/93
      call cnoblock(ischan,0,ires,ierrno)
c

      if(iwt.ge.0) call cnoblock(ischan,0,ires,ierrno)
      if(lr.lt.1.or.response(lr:lr).ne.'>') then
        lr=-1
      else
        call oasstrip(response,lr)
      endif
      if (lr.gt.0.and.ioasverbose.eq.2)
     &     write(*,'(a,a)') 'recvd from OAS: ',response(1:lr)
      return
      end
