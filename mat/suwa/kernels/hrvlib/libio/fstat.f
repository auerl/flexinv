c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine fstat(lufl,lent,iftyp,isize)
      include "openfile.h"
        call cfstat(jchn(lufl),isize,iftyp,ierrno)
c        if(ierrno.ne.0) call check('cfstat on fstat')
        lent=jrecl(lufl)
        return
      end
