c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine lower(string)
      character*(*) string
      ll=len(string)
      do 10 i=1,ll
      ic=ichar(string(i:i))
      if(ic.ge.65.and.ic.le.90) ic=ic+32
   10 string(i:i)=char(ic)
      return
      end
