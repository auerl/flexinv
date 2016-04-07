      subroutine oasopenserial(ierror)
      include "oasstatus.h"
      character*1 null
      parameter (null=char(0))
      character*200 mess
      character*30 dev
      character*80 safdir
      character*80 string
c
c---- get the name of the SAF directory and the serial port
c
      call getsafdir(safdir,lsafdir)
c
      open(99,file=safdir(1:lsafdir)//'/jukeconfig',status='read')
      read(99,"(i1,1x,a30)",iostat=ios) n,dev
      if(ios.ne.0) call exit(2)
      close(99)
      ldev=istlen(dev)
c
c---- try to open the serial line
c
      itry = 0
 10   call copen(dev(1:ldev)//null,ischan,2,ierrno,0,0)

      if(ierrno.ne.0) then
        if (ierrno.eq.16.and.ioaswait.eq.1) then
          if(itry.eq.0) then
             itry = 1
             write(*,"('serial port busy: waiting ...')")
           endif
           call wait(4000)
           goto 10
        endif

        ierror=ierrno
      else
        ierror=0
c
c---- make the serial line exclusive to this process
c
        call cgtflag(ischan,iflag,ires,ierrno)
        iflag=o'341'
        call cstflag(ischan,iflag,13,ires,ierrno)
        call cgtlmw(ischan,iword,ires,ierrno)
        call cstlmw(ischan,z'4000',ires,ierrno)
        call cstdtr(ischan,ires,ierrno)
        call cstexc(ischan,ires,ierrno)
c
c---- flush the serial port
c
        call oasflushserial()
c
c---- check for remaining errors
c
        do iuser=0,1
          write(string,"('user/default=',i1)") iuser
          lstring=lnblnk(string)
          call oassend(string(1:lstring),mess,lmess,-1)
          if(lmess.ne.0) then
              write(6,"('OAS - ',a)") mess(1:lmess)
          endif
          call oassend('status/system/type=9',mess,lmess,-1)
          if(lmess.ne.0) then
c              write(6,"('OAS tape status - ',a)") mess(1:lmess)
          endif
        enddo
  
c
c---- show the name legality check
c
        call oassend('show/name',mess,lmess,-1)
        if(lmess.ne.0) then
          if(ioasverbose.ge.2)  write(6,"('OAS - duplicate tape names allowed?',
     #               ' yes(1)/no(0): ',a)")mess(1:lmess)
        endif
c          if(mess(1:1).eq.'1') then
c            call oassend('set/name=0',mess,lmess,-1)
c            if(lmess.ne.0) then
c              write(6,"('OAS - ',a)")mess(1:lmess)
c            endif
c            call oassend('show/name',mess,lmess,-1)
c            if(lmess.ne.0) then
c              write(6,"('OAS - duplicate tape names allowed? yes(1)/no(0): ',a)")mess(1:lmess)
c            endif
c          endif
c        endif
      endif
      return
      end
