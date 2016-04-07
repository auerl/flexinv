      subroutine oasverbose(iverb)
      include 'oasstatus.h'

      if (iverb.lt.0.or.iverb.gt.2) then
        write(*,*) 'invalid verbosity option = ',iverb
        write(*,*) ' 0 = no oas echoing'
        write(*,*) ' 1 = echo strings sent to oas'
        write(*,*) ' 2 = echo strings received from oas'
        return
      endif
      ioasverbose = iverb
      return
      end

     
