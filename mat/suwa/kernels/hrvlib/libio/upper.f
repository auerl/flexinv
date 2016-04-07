c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine upper(string)
      character*(*) string
      ll=len(string)
      do 10 i=1,ll
      ic=ichar(string(i:i))
      if(ic.gt.95) ic=ic-32
   10 string(i:i)=char(ic)
      return
      end
