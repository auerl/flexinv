      subroutine getjukehost(string,ll)
      character*(*) string
      character*80 string2
      call getsafdir(string2,ll)
      open(99,file=string2(1:lnblnk(string2))//'/jukehost')
      read(99,"(a)") string
      close(99)
      ll=lnblnk(string)
      return
      end
