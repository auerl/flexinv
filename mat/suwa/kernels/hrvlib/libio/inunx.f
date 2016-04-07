c---------------------------------------------------------------
      integer*4 function inunx(id,iseq,nbyts)
      character*(*) id
      character*80 string,getunx
      string=getunx(id,iseq,nbyts)
      if(nbyts.gt.0) then
        read(string,*) itemp
      else




        itemp=0

      endif

      inunx=itemp
      return
      end
