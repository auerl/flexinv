c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffos(lufl,ifbin,ibuf,nbytes,istat,irec)
      include "openfile.h"
      dimension ibuf(1)
      krec=irec-1
      if(irec.eq.0) krec=jrec(lufl)
      kchn=jchn(lufl)
      if(jrecl(lufl).eq.0) then
c        write(6,"('krec,jrec',2i10)") krec,jrec(lufl)
        if(krec.lt.jrec(lufl)) then
          call cmtio(jchn(lufl),4,jrec(lufl)-krec,ires,ierrno)
        else if(krec.gt.jrec(lufl)) then
          call cmtio(jchn(lufl),3,krec-jrec(lufl),ires,ierrno)
        endif
        call cwrite(jchn(lufl),ibuf,nbytes,ires,ierrno)
        if(ierrno.ne.0.or.ires.ne.nbytes) then
          write(6,"('jchn,ierrno',2i10)") jchn(lufl),ierrno
          print *,hello
          write(6,"('nbytes,ires',2i10)") nbytes,ires
          call check('cwrite in bffos 1')
        endif
        istat=2
        jrec(lufl)=1+krec
        return
      else
        nb=nbytes
        if(jfile(lufl).ne.200) nb=min0(nbytes,jrecl(lufl))
        call clseek(jchn(lufl),jrecl(lufl)*krec,0,ires,ierrno)
        if(ierrno.ne.0) call check('clseek in bffos')
        call cwrite(jchn(lufl),ibuf,nb,ires,ierrno)
        if(ierrno.ne.0.or.ires.ne.nbytes) call check('cwrite in bffos 2')
        istat=2
        jrec(lufl)=krec+(nb+jrecl(lufl)-1)/jrecl(lufl)
        lenglu(lufl)=max0(lenglu(lufl),jrec(lufl))
      endif
      return
      end
