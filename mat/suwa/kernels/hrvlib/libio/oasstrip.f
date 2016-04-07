      subroutine oasstrip(mess,lmess)
      character*1 null,nl
      parameter (null=ichar(0))
      parameter (nl=ichar(10))
      character*(*) mess
      do while (lmess.gt.0)
        if(mess(lmess:lmess).ne.'>'.and.mess(lmess:lmess).ne.nl
     1  .and.mess(lmess:lmess).ne.null.and.mess(lmess:lmess).ne.' ') goto 10
        lmess=lmess-1
      enddo
   10 continue
      return
      end
