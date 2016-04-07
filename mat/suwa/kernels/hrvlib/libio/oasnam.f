      character*80 function oasnam(string,lname)
      character*(*) string
      character*1 cc
      logical quote
      j=0



      quote=.FALSE.


      do i=1,len(string)
        cc=string(i:i)
        if(cc.eq.'"') then
          write(0,"('oasnam: double quote in oas name (locally) illegal')")
          call exit(2)
        endif

        if(   cc.eq.' '.or.cc.eq.':'.or.cc.eq.'/'
     1                 .or.cc.eq.'['.or.cc.eq.'['
     1                 .or.cc.eq.'>'.or.cc.eq.'<') then
          if(.not.quote) then
            j=j+1


            oasnam(j:j)='"'
            quote=.TRUE.
          endif
        else




          if(quote) then
            j=j+1

            oasnam(j:j)='"'
            quote=.FALSE.
          endif
        endif
          
          

        j=j+1



        oasnam(j:j)=cc

      enddo
      if(quote) then
        j=j+1
        oasnam(j:j)='"'
      endif
      if(j.ge.2.and.oasnam(j-1:j).eq.';"') oasnam(j-1:j)='";'

      lname=j
      return
      end
