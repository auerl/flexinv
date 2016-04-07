c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      function istlen(string)
      character*1 null
      parameter(null=char(0))
      character*(*) string
      k=len(string)
      do 10 i=1,k
      j=k+1-i
      if(string(j:j).eq.' '.or.string(j:j).eq.null) goto 10
      istlen=j
      goto 99
   10 continue
      istlen=0
   99 return
      end
