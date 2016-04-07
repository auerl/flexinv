
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine fstatn(name,lent,iftyp,isize)
      character*(*) name
      include "openfile.h"
        lu=-1

        call openfl(lu,name,1,0,0,istat,-1)
        call cfstat(jchn(lu),isize,iftyp,ierrno)
        lent=jrecl(lu)
        call closfl(lu,istat)
        return
      end
