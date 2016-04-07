      subroutine ilbyte(k,buf,ibyt)
      character*1 buf(*)
      k=ichar(buf(ibyt+1))
      return
      end
