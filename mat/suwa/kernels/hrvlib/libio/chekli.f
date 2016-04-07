c-------------------------------------------------------------------
      subroutine chekli(string,ierr)
      character*(*) string

      include "getgnl.h"

      include "getli.h"

      external iargli,getali

      call chekgl(string,iargli,getali,ierr
     1  ,swnam,cline,deflt
     1  ,nswt,lswnam,ldeflt,ifreq,nmin,nmax,icnt,iopn,iptr,itable)
      end
