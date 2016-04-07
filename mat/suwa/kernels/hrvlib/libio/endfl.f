c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine endfl(lufl)
      include "openfile.h"
      if(jrecl(lufl).ne.0) then
        leng=jrec(lufl)*jrecl(lufl)
        call clseek(jchn(lufl),leng,0,ires,ierrno)
        if(ierrno.ne.0.or.ires.ne.leng) call check('clseek in endfl')
        call ctrun(jchn(lufl),leng,ierrno)
        if(ierrno.ne.0) call check('ctrun in endfl')
      else
        call cmtio(jchn(lufl),0,1,ires,ierrno)
        if(ierrno.ne.0) call check('cmtio in endfl')
        jfile(lufl)=1+jfile(lufl)
        jrec(lufl)=0
      endif
      return
      end
