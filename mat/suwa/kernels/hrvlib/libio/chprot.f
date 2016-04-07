      subroutine chprot(string)
      character*(*) string
      do i=1,len(string)
        string(i:i)=char(xor(ichar(string(i:i)),z'80'))
      enddo

      return
      end
