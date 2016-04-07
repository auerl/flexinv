c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffin(lufl,ifbin,ibuf,nwords,istat,nread)
      dimension ibuf(1)
      nbytes=4*nwords
      call bffi(lufl,ifbin,ibuf,nbytes,istat,nread,0)
      nread=nread/4
      return
      end
