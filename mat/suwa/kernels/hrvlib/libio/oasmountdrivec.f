      subroutine oasmountdrivec(idrive,volume,tape,irdwr,ierror)
      character*(*) volume
      character*(*) tape
      character*80 string
      character*200 mess
      character*1 quote
      parameter (quote='"')
      include 'oasstatus.h'
c
c---- check read/write flag
c
      ierror=0
      if(irdwr.ne.0.and.irdwr.ne.1) then
        ierror=1
        return
      endif
c
c---- check that tape name if provided
c
      lvol=lnblnk(volume)
      ltap=lnblnk(tape)
      if(ltap.le.0) then
        ierror=2
        return
      endif
c
c---- if this is a blank tape image, the volume must be given
c
      if(lvol.le.0.and.irdwr.eq.1) then
        ierror=5
        return
      endif
c
c---- first make sure that the drive is owned
c
      if(idrive.eq.1.or.idrive.eq.2) then
        if(iowned(idrive).eq.1) then
          
c
c---- set the default drive
c
          write(string,"('user/default=',i1)") idrive-1
          call oassend(string,mess,lmess,-1)
          if(lmess.ne.0) then
            write(6,*) 'From OAS after user: ',mess(1:lmess)
            pause
          endif
c
c---- go off-line and dismount the tape (if mounted)
c
          call oassend('offline',mess,lmess,-1)
          if(lmess.ne.0) then
            write(6,*) 'From OAS after offline: ',mess(1:lmess)
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
c
c---- mount the appropriate tape for reading
c
          if(irdwr.eq.0) then
            if(lvol.gt.0) then
              string='mount '//quote//volume(1:lvol)//quote//';'//quote//tape(1:ltap)//quote
            else
              string='mount '//quote//tape(1:ltap)//quote
            endif
            call oassend(string,mess,lmess,-1)
            if(lmess.ne.0) then
              write(6,*) 'From OAS after tape mount (read): ',mess(1:lmess)
              if(mess(1:lmess).eq.'E-OMSG-DISK IN USE') then
                ierror=3
                return
              else if(mess(1:lmess).eq.'E-DMGR-TAPE NOT FOUND') then
                ierror=4
                return
              else if(mess(1:lmess).eq.'E-JUKE-VOLUME NOT IN') then
                ierror=4
                return
              else if(mess(1:lmess).eq.'E-OMSG-NAME ALREADY EXISTS') then
                ierror=10
                return
              else
                pause
              endif
            endif
c
c---- mount the appropriate tape for writing
c
          else
            if(lvol.gt.0) then
              string='mount '//quote//volume(1:lvol)//quote//';'
              write(6,"(a)") string(1:lnblnk(string))
              call oassend(string,mess,lmess,-1)
              if(lmess.ne.0) then
                write(0,*) 'From OAS after volume mount (read): ',mess(1:lmess)
                if(mess(1:lmess).eq.'E-OMSG-DISK IN USE') then
                  ierror=3
                  return
                else if(mess(1:lmess).eq.'E-DMGR-TAPE NOT FOUND') then
                  ierror=4
                  return
                else if(mess(1:lmess).eq.'E-JUKE-VOLUME NOT IN') then
                  ierror=4
                  return
                else
                  pause
                endif
              endif
c
c---- test to see if compression screws things up
c
              string='mount/cblank '//quote//tape(1:ltap)//quote
c               string='mount/blank '//quote//tape(1:ltap)//quote
c
c              write(6,"(a)") string(1:lnblnk(string))
              call oassend(string,mess,lmess,-1)
              if(lmess.ne.0) then
                write(0,*) 'From OAS after tape mount (write): ',mess(1:lmess)
                if(mess(1:lmess).eq.'E-OMSG-DISK IN USE') then
                  ierror=3
                  return
                else if(mess(1:lmess).eq.'E-DMGR-TAPE NOT FOUND') then
                  ierror=4
                  return
                else if(mess(1:lmess).eq.'E-OMSG-NAME ALREADY EXISTS') then
                  ierror=10
                  return
                else
                  pause
                endif
              endif
            else
              ierror=11
              return
            endif
          endif
c
          imounted(idrive)=1
          volumemounted(idrive)=volume(1:lvol)
          tapemounted(idrive)=tape(1:ltap)
        else
          write(6,"('the optical drive is not owned by you:',i2)") idrive
          ierror=1
        endif
      else
        write(6,"('mountoasvolume: illegal drive',i2)") idrive
        call exit(2)
      endif
      return
      end
