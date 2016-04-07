      subroutine fillu(lu,name)
      include "openfile.h"
      character*(*) name
      name=' '
      if(jchn(lu).gt.0) name=opnnam(lu)
      return
      end
