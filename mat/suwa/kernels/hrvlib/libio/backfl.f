c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine backfl(lufl)
      include "openfile.h"
      if(jrecl(lufl).eq.0) then
        if(jrec(lufl).eq.0) then
          if(jfile(lufl).le.1) then
            call rewtp(lufl)
            jfile(lufl)=0
            jrec(lufl)=0
            return
          else
            call cmtio(jchn(lufl),2,2,ires,ierrno)
            if(ierrno.ne.0) call check('cmtio:backfile in backfl')
            call cmtio(jchn(lufl),1,1,ires,ierrno)
            if(ierrno.ne.0) call check('cmtio:forwardfile in backfl')
            jrec(lufl)=0
            jfile(lufl)=jfile(lufl)-1
          endif
        endif
      else 
        call rewfl(lufl)
      endif
      return
      end
