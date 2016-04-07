      subroutine oasopendrive(idrive,drive,ierror)
      character*80 string
      character*(*) drive
      character*200 mess
      character*1 quote,null
      parameter (quote='"')
      parameter (null=ichar(0))
      include 'oasstatus.h'
c
      ierror=0
c
c---- first make sure that the drive is owned
c
      if(idrive.eq.1.or.idrive.eq.2) then
        if(iowned(idrive).eq.1) then
          if(imounted(idrive).eq.1) then
c
c---- set the default drive
c
            write(string,"('user/default=',i1)") idrive-1
            call oassend(string,mess,lmess,-1)
            if(lmess.ne.0) then
                write(0,*) 'From OAS after user: ',mess(1:lmess)
            endif
c
c---- go on-line
c
            call oassend('online',mess,lmess,-1)
            if(lmess.ne.0) then
                write(0,*) 'From OAS after online: ',mess(1:lmess)
            endif
            drive=namedrive(idrive)
          endif
        endif
      endif
      return
      end
