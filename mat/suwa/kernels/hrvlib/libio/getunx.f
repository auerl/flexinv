c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      character*80 function getunx(ident,nni,nbyts)
      character*(*) ident
      character*80 getgnl

      include "getgnl.h"
cc      character*80 cline
cc      character*30 swnam(30)
cc      dimension lswnam(30),ifreq(30),nmin(30),nmax(30)
cc     1         ,icnt(30),iopn(30),iptr(30),itable(100,3)

      include "getunx.h"
cc      common/ccmlin/swnam,cline
cc      common/icmlin/nswt,lswnam,ifreq,nmin,nmax
cc     1         ,icnt,iopn,iptr,itable

      getunx=getgnl(ident,nni,nbyts
     1  ,swnam,cline,deflt
     1  ,nswt,lswnam,ldeflt,ifreq,nmin,nmax
     1  ,icnt,iopn,iptr,itable)

      return
      end
