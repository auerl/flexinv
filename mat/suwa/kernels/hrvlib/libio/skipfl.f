c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine skipfl(lufl)
      include "openfile.h"
      if(jrecl(lufl).eq.0) then
        call cmtio(jchn(lufl),1,1,ires,ierrno)
        if(ierrno.ne.0) call check('cmtio:forwardfile in skipfl')
        jfile(lufl)=jfile(lufl)+1
      endif
      jrec(lufl)=0
      return
      end
