      subroutine getsafdir(string,ll)
      character*(*) string
      include "safdir.h"
      string=safdir
      ll=lnblnk(string)
      return
      end
