c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine adnul(string1,string2)
      character*1 null
      parameter (null=char(0))
      character*(*) string1,string2
      string2=string1
      ip0=len(string2)
      ip=ip0
      do while (ip.gt.0.and.string2(ip:ip).eq.' '
     1     .or.string2(ip:ip).eq.null)
        ip=ip-1
      enddo

      ip=1+ip
      if(ip.le.ip0) string2(ip:ip)=null
      return
      end
