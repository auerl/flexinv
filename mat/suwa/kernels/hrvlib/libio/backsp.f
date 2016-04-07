c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine backsp(lufl)
      include "openfile.h"
      if(jrecl(lufl).eq.0) then
        if(jrec(lufl).le.0) then
          call rewfl(lufl)
          jrec(lufl)=0
        else
          call cmtio(jchn(lufl),4,1,ires,ierrno)
          if(ierrno.ne.0) call check('cmtio:backrec in rewfl')
          jrec(lufl)=jrec(lufl)-1
        endif
      else
        jrec(lufl)=max0(0,jrec(lufl)-1)
      endif
      return
      end
