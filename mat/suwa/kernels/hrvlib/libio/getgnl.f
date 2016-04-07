
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      character*80 function getgnl(ident,nni,nbyts
     1  ,swnam,cline,deflt
     1  ,nswt,lswnam,ldeflt,ifreq,nmin,nmax
     1  ,icnt,iopn,iptr,itable)
      character*(*) ident
      
      include "getgnl.h"

      character*1 null
      parameter (null=char(0))

      nn=max0(1,nni)
      ll=len(ident)
      do while (ll.gt.0.and.(ident(ll:ll).eq.' '
     1          .or.ident(ll:ll).eq.null))
        ll=ll-1
      enddo
      imatch=0
      iswt=0
      ierr=0
      do i=1,nswt
        if((ll.eq.0.and.lswnam(i).eq.0)
     1   .or.(ll.gt.0.and.ident(1:ll).eq.swnam(i)(1:ll))) then
          imatch=1+imatch
          if(iswt.eq.0) iswt=i
        endif
      enddo
      if(imatch.eq.0) then
        ierr=1
        goto 99
      else if(imatch.ne.1) then
        ierr=2
        goto 99
      else




        if(icnt(iswt).lt.0) then
          nbyts=-1
          getgnl=' '
          if(ldeflt(iswt).ge.0) then
            ii=1
 



            ik=0


            do while(ii.le.ldeflt(iswt).and.ik.lt.nn)
            do while(ii.le.ldeflt(iswt).and.deflt(iswt)(ii:ii).eq.' ')
              ii=ii+1
            enddo

            ik=ik+1

            ij=0


            do while(ii.le.ldeflt(iswt).and.deflt(iswt)(ii:ii).ne.' ')
              ij=ij+1
              if(ik.eq.nn) then
                nbyts=ij
                getgnl(ij:ij)=deflt(iswt)(ii:ii)
              endif
              ii=ii+1
            enddo



            enddo
          
          endif


          return
        else if(nn.gt.icnt(iswt)) then
          nbyts=0
          getgnl=' '
          return
        else




          ipt=iptr(iswt)
          ipt0=ipt
          do i=1,nn-1
            ipt=itable(ipt,3)
            if(ipt.eq.ipt0) then
              ierr=3
              goto 99
            endif
          enddo
          ip1=itable(ipt,1)
          ip2=itable(ipt,2)
          nbyts=ip2-ip1+1
          getgnl=cline(ip1:ip2)
          return
        endif
      endif
   99 continue
      goto (1,2,3),ierr
    1 write(0,"('getgnl: no match')")
      goto 97
    2 write(0,"('getgnl: non-unique abbreviation for identifier')")
      goto 97
    3 write(0,"('getgnl: required argument not found in table')")
      goto 97
   97 call exit(2)
      end
