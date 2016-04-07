c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine rewtp(lufl)
      include "openfile.h"
      if(jrecl(lufl).eq.0) then
        call cmtio(jchn(lufl),5,1,ires,ierrno)
        if(ierrno.ne.0) call check('cmtio:rewind in rewtp')
        jrec(lufl)=0
        jfile(lufl)=0
      else
        call rewfl(lufl)
      endif
      return
      end
