      subroutine getflpos(lufl,iposition,istatus)
      include "openfile.h"
      iposition=jrec(lufl)
      istatus=0
      return
      end
