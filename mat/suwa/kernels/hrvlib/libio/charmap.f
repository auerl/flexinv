      character*80 function charmap(string,nchar)
      integer*4 ia(3)
      character*(*) string
      if(string.eq.'TIME') then
        call itime(ia)
        write(charmap,"(i2,2(':',i2))") (ia(i),i=1,3)
        nchar=8
        do i=1,nchar
          if(charmap(i:i).eq.' ') charmap(i:i)='0'
        enddo
      else if(string.eq.'DATE') then
        call idate(ia)
        write(charmap,"(i2,2('/',i2))") ia(2),ia(1),mod(ia(3),100)
        nchar=8
        do i=1,nchar
          if(charmap(i:i).eq.' ') charmap(i:i)='0'
        enddo
      else
        charmap='UNKNOWN'
        nchar=7
      endif
      return
      end
