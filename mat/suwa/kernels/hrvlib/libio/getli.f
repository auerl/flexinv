c-------------------------------------------------------------------
      character*80 function getli(ident,nni,nbyts)
      character*(*) ident
      character*80 getgnl

      include "getgnl.h"

      include "getli.h"

      getli=getgnl(ident,nni,nbyts
     1  ,swnam,cline,deflt
     1  ,nswt,lswnam,ldeflt,ifreq,nmin,nmax
     1  ,icnt,iopn,iptr,itable)

      return
      end
