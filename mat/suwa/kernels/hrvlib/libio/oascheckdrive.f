      subroutine oascheckdrive(idrive,ierror)
      include 'oasstatus.h'
      character*80 string
      character*200 mess
c
      ierror=0
c
c---- first make sure that the drive is owned
c
      if(idrive.eq.1.or.idrive.eq.2) then
        if(iowned(idrive).eq.1) then
          write(string,"('user/default=',i1)") idrive-1
          call oassend(string,mess,lmess,-1)
          if(lmess.ne.0) then
            write(0,*) 'From OAS after user: ',mess(1:lmess)
            pause
          endif
c
          call oassend('status/system/type=9',mess,lmess,-1)
          if(lmess.ne.0) then
c              write(6,"('OAS tape status - ',a)") mess(1:lmess)
          endif
c
          string='status/volume/type=0'
          call oassend(string,mess,lmess,-1)
          if(lmess.ne.0) then
            if(mess(1:lmess).eq.volumemounted(idrive)(1:lmess).or.
     #           lnblnk(volumemounted(idrive)).eq.0) then
            else
              write(6,*) 'From OAS after volume inquiry: ',mess(1:lmess)
              ierror=1
              pause
            endif
          endif
          string='status/tape/type=0'
          call oassend(string,mess,lmess,-1)
          if(lmess.ne.0) then
            if(mess(1:lmess).eq.tapemounted(idrive)(1:lmess)) then
            else
              write(0,*) 'From OAS after tape inquiry  : ',mess(1:lmess)
              ierror=1
              pause
            endif
          endif
        endif
      endif
      return
      end
       
