      subroutine rcreat(path)
      character*(*) path
      call rcreao(path,ichan)
      if(ichan.gt.0) call cclose(ichan,ires,ier)
      return
      end
