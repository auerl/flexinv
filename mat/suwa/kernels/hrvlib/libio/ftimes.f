      subroutine ftimes(file,isgmta,isgmtm)
      dimension itime(4)
      character*(*) file
      character*80 file1
      call adnul(file,file1)
      itime(1)=isgmta
      itime(2)=0
      itime(3)=isgmtm
      itime(4)=0
      call cutimes(file1,itime,ires,ierrno)
      if(ierrno.ne.0) call check('ftimes')
      return
      end
