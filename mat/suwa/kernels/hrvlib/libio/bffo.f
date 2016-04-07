c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffo(lufl,ifbin,ibuf,nbytes,istat,irec)
      include "openfile.h"
      dimension ibuf(1)
      if(jrecl(lufl).eq.1.and.jfile(lufl).eq.200) then
        nread=0
        nget=nbytes
        irc=irec
        if(irc.eq.0) irc=jrec(lufl)+1
        iad=1
   10   nb=min0(nget,65532)
        call bffos(lufl,ifbin,ibuf(iad),nb,istat,irc)
        nget=nget-nb
        if(nget.ne.0) then
          nb4=nb/4
          irc=irc+nb
          iad=iad+nb4
          goto 10
        endif
      else
        call bffos(lufl,ifbin,ibuf,nbytes,istat,irec)
      endif
      return
      end
