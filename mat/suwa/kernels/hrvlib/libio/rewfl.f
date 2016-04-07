c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine rewfl(lufl)
      include "openfile.h"
      if(jrecl(lufl).eq.0) then
        if(jfile(lufl).le.0) then
           call rewtp(lufl)
           jfile(lufl)=0
        else
          call cmtio(jchn(lufl),2,1,ires,ierrno)
          if(ierrno.ne.0) call check('cmtio:backfile in rewfl')
          call cmtio(jchn(lufl),1,1,ires,ierrno)
          if(ierrno.ne.0) call check('cmtio:forwardfile in rewfl')
        endif
      endif
      jrec(lufl)=0
      return
      end
