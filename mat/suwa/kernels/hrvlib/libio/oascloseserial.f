      subroutine oascloseserial(ierror)
      include "oasstatus.h"
      call cclose(ischan,ires,ierrno)
      ierror=ierrno
      return
      end
