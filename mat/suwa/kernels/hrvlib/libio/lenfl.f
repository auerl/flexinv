c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      function lenfl(lufl)
      include "openfile.h"
      if(jrecl(lufl).eq.0) then
        lenfl=0
      else
        call fstat(lufl,idummy,idummy1,isize)
        lenfl=isize
      endif
      return
      end
