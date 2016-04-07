c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine chekcl(string)
      character*(*) string

      include "getgnl.h"
cc      character*80 cline
cc      character*30 swnam(30)
cc      dimension lswnam(30),ifreq(30),nmin(30),nmax(30)
cc     1         ,icnt(30),iopn(30),iptr(30),itable(100,3)

      include "getunx.h"
cc      common/ccmlin/swnam,cline
cc      common/icmlin/nswt,lswnam,ifreq,nmin,nmax
cc     1         ,icnt,iopn,iptr,itable

      external iargc,getarg

      call chekgl(string,iargc,getarg,ierr
     1  ,swnam,cline,deflt
     1  ,nswt,lswnam,ldeflt,ifreq,nmin,nmax,icnt,iopn,iptr,itable)
      if(ierr.ne.0) call exit(2)
      return
      end
