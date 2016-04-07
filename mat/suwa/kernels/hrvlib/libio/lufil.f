      function lufil(name)
      include "openfile.h"
      character*(*) name
      lufil=-1
      do i=1,nlu
        if(name.eq.opnnam(i).and.jchn(i).gt.0) then
          lufil=i
          goto 10
        endif
      enddo
   10 return
      end
