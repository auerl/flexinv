      subroutine oasflushserial(ierror)
      integer EWOULDBLOCK
      parameter (EWOULDBLOCK=35)
      include "oasstatus.h"
      character*1 cc
c
c---- read from serial line until no more characters
c---- first set descriptor for non-blocking i/o
c
      call cnoblock(ischan,1,ires,ierrno)
      ierrno=0
      do while(ierrno.ne.EWOULDBLOCK)
        call cread(ischan,cc,1,ires,ierrno)
      enddo
c
c---- clear for non-blocking i/o
c
      call cnoblock(ischan,0,ires,ierrno)
      end
