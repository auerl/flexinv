c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffi(lufl,ifbin,ibuf,nbytes,istat,nread,irec)
      include "openfile.h"
      dimension ibuf(1)
      if(jrecl(lufl).eq.1.and.jfile(lufl).eq.200) then
        nread=0
        nget=nbytes
        irc=irec
        if(irc.eq.0) irc=jrec(lufl)+1
        iad=1
   10   nb=min0(nget,65532)
        call bffis(lufl,ifbin,ibuf(iad),nb,istat,nr,irc)
        if(istat.eq.2.and.nr.eq.nb) then
          if(nr.eq.nget) then
            nread=nread+nr
          else
            nr4=nr/4
            nr=4*nr4
            irc=irc+nr
            nread=nread+nr
            iad=iad+nr4
            nget=nget-nr
            goto 10
          endif
        else
          nread=nread+nr
        endif
      else
        call bffis(lufl,ifbin,ibuf,nbytes,istat,nread,irec)
      endif
      return
      end
