c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine check(text)
      character*(*) text
      call cperror(text)
      call exit(2)
      end
