      logical function error(ierrno)
      call cgterr(ierrno)
      if(ierrno.eq.0) then
        error=.FALSE.
      else
        error=.TRUE.
      endif
      return
      end
