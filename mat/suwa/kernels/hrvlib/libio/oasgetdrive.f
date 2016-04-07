      subroutine oasgetdrive(idrive,ierror)
      include "oasstatus.h"
      character*1 null
      character*200 mess
      character*80 string
      parameter (null=char(0))

      ierror = 0
      itry   = 0
c
c---- attempt to acquire the appropriate tape drive
c
  10  continue
      if(idrive.eq.1.or.idrive.eq.2) then
	write(6,"('trying to open',a)") nameblock(idrive)
        if(ioasverbose.ge.2) write(6,"('trying to open: ',a)")nameblock(idrive)
        call copen(nameblock(idrive)//null,itch,0,ierrno,0,0)
        if(ierrno.eq.0) then
          if(ioasverbose.ge.2)write(6,"('closing: ',a)")nameblock(idrive)
          call cclose(itch,ires,ierr)
          ierror=0
        else if( (ierrno.eq.6.or.ierrno.eq.16).and.ioaswait.eq.1) then
          if (itry.eq.0) then
            write(6,'(a,i1,a)') 'oasgetdrive: tape device # ',idrive,
     &                          ' in use: waiting ...'
            itry = 1
          endif
          call wait(4000)
          goto 10
        else if( (ierrno.eq.6.or.ierrno.eq.16).and.ioaswait.eq.0) then
          ierror = ierrno
        else if(ierrno.eq.5) then
          if(ioasverbose.ge.2) write(6,"('already off-line: ',a)")nameblock(idrive)
          ierror=0
        else
          write(6,"('oasgetdrive: ierrno',i5,1x,a)") ierrno,nameblock(idrive)
          ierror=ierrno
          call exit(2)
        endif
c
c---- if the tape was opened and closed (ierror=0) or off-line (ierror=5) we can
c---- now manipulate the optical drive
c
        if(ierror.eq.0) then
          iowned(idrive)=1
c
c---- set the default drive
c
          write(string,"('user/default=',i1)") idrive-1
          lstring=lnblnk(string)
          call oassend(string(1:lstring),mess,lmess,-1)
          if(lmess.ne.0) then
            write(0,*) 'From OAS after user: ',mess(1:lmess)
            pause
          endif
c
c---- go off-line and dismount the tape (if mounted)
c
          call oassend('offline',mess,lmess,-1)
          if(lmess.ne.0) then
            write(0,*) 'From OAS after offline: ',mess(1:lmess)
            pause
          endif
          call oassend('dismount',mess,lmess,-1)
          if(lmess.ne.0) then
            if(mess(1:lmess).eq.'E-OMSG-DISK NOT MOUNTED') then
            else if(mess(1:lmess).eq.'E-OMSG-TAPE NOT MOUNTED') then
            else
              write(0,*) 'From OAS after dismount: ',mess(1:lmess)
              pause
            endif
          endif
        else
          iowned(idrive)=0
        endif
        imounted(idrive)=0
      else
        write(6,"('illegal argument for oasgetdrive',i10)") idrive
        call exit(2)
      endif
      return
      end
