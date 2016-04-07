      subroutine oaswait(iwait)
      include 'oasstatus.h'

      if (iwait.lt.0.or.iwait.gt.1) then
        write(*,*) 'invalid wait option = ',iwait
        write(*,*) ' 0 = no wait'
        write(*,*) ' 1 = wait'
        return
      endif
      ioaswait = iwait
      return
      end

     
