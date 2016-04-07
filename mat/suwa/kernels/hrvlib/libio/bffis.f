c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffis(lufl,ifbin,ibuf,nbytes,istat,nread,irec)
      include "openfile.h"
      dimension ibuf(1)
      krec=irec-1
      if(irec.eq.0) krec=jrec(lufl)
      if(jrecl(lufl).eq.0) then
        if(krec.lt.jrec(lufl)) then
c	  write(6,"('before positioning tape -- backing up',i6)") 
c     #     jrec(lufl)-krec
          call cmtio(jchn(lufl),4,jrec(lufl)-krec,ires,ierrno)
c          write(6,"('positioning tape: krec=',i6,'  jrec=',i6)")
c     1       krec,jrec(lufl)
        else if(krec.gt.jrec(lufl)) then
c	  write(6,"('before positioning tape -- skipping',i6)")
c     #     krec-jrec(lufl)
          call cmtio(jchn(lufl),3,krec-jrec(lufl),ires,ierrno)
c          write(6,"('positioning tape: krec=',i6,'  jrec=',i6)")
c     1       krec,jrec(lufl)
        endif
        call cread(jchn(lufl),ibuf,nbytes,nread,ierrno)
        if(nread.eq.255) nread=0
        istat=2
        if(ierrno.eq.0) then
          if(nread.eq.0) istat=3
        else
          istat=5
          write(6,"(6i10)") lufl,jchn(lufl),nbytes,nread,krec,ierrno
          call check('cread in bffis 1')
        endif
        if(istat.eq.3) goto 30
        jrec(lufl)=1+krec
        return
   30   jfile(lufl)=1+jfile(lufl)
        jrec(lufl)=0
      else
      nb=nbytes
      if(jfile(lufl).ne.200) nb=min0(nb,jrecl(lufl))
        call clseek(jchn(lufl),jrecl(lufl)*krec,0,ires,ierrno)
        if(ierrno.ne.0) call check('clseek in bffis')
        call cread(jchn(lufl),ibuf,nb,nread,ierrno)
        if(ierrno.ne.0) then
          write(6,"(7i10)") lufl,jfile(lufl),jrecl(lufl),jchn(lufl),nb,ierrno,nread
          call check('cread in bffis 2')
        endif
        istat=2
        if(nread.eq.0) istat=3
        if(istat.ne.3) then
          jrec(lufl)=krec+(nread+jrecl(lufl)-1)/jrecl(lufl)
        else
          jrec(lufl)=krec
        endif
      endif
      return
      end
