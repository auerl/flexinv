      function ichdec(chr,ktype)
c Decodes characters in ray name, returning an integer.
      character*1 chr
      character*10 okchr                                              
      okchr(1:10)=' PSIJKpsic'
      do 1 i=1,10
      ichdec=i
      if(chr.eq.okchr(i:i)) goto 2
    1 continue
      ichdec=11                 
      return
    2 ktype=1
      if(ichdec.eq.3.or.ichdec.eq.5.or.ichdec.eq.8) ktype=2
      return
      end
