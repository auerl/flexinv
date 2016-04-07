c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------

      subroutine chekgl(string,iargf,gtargf,ierr
     1  ,swnam,cline,deflt
     1  ,nswt,lswnam,ldeflt,ifreq,nmin,nmax,icnt,iopn,iptr,itable)
      character*(*) string

      include "getgnl.h"
cc     character*1000 cline
cc     character*30 swnam(30),deflt(30)
cc     dimension lswnam(30),ldeflt(30),ifreq(30),nmin(30),nmax(30)
cc    1         ,icnt(30),iopn(30),iptr(30),itable(100,3)

      character*1 null
      parameter (null=char(0))
      dimension mlen(30)
      integer*4 mess(30),mss
      character*30 swt,swtest
      character*1 c
      character*8 str1
      character*21 str2
      character*10 str3
      data nbig/99999999/


      ierr=0
      nswt=-1
      swt=null
      ll=len(string)
      do ip=1,ll
        c=string(ip:ip)
        if(c.eq.'|') then
          do k=1,nswt
            if(swt.eq.swnam(k)) then
              ierr=5
              goto 52
            endif
          enddo
          nswt=1+nswt
          if(nswt.gt.0) then
            swnam(nswt)=swt
            lswnam(nswt)=lswt
            ifreq(nswt)=ifrq
            nmin(nswt)=nmn
            nmax(nswt)=nmx
            mess(nswt)=mss
            mlen(nswt)=mln

            ldef=-1
            if(ifrq.eq.0) then
              mee=mss+mln
              do while(mss.lt.mee.and.string(mss:mss).ne.'[')
                mss=mss+1
              enddo

              if(mss.le.mee) ldef=0
              mss=mss+1
              do while(mss.lt.mee.and.string(mss:mss).ne.']')
                ldef=ldef+1
                deflt(nswt)(ldef:ldef)=string(mss:mss)
                mss=mss+1
              enddo
            endif

            ldeflt(nswt)=ldef
              
        
          endif
          swt=null
          lswt=0
          ifrq=-1
          nmn=0
          nmx=0
          mode=0
          mode1=0
          num=0

          mss=-1
          mln=0
        else if(c.eq.':') then
          mode=mode+1
        else if(c.eq.','.and.mode.ne.3) then
          if(mode.ne.2) then
            write(0,*) 'chekcl: format error 1'
            call exit(2)
          endif
          num=0
          mode1=1
        else if(c.eq.' '.and.mode.ne.3) then
        else
          if(mode.eq.0) then
            lswt=1+lswt
            swt(lswt:lswt)=c
          else if(mode.eq.1) then
            if(c.eq.'r') then
              ifrq=1
            else if(c.eq.'o') then
              ifrq=0
            else
              write(0,*) 'chekcl: format error 2'
              call exit(2)
            endif
          else if(mode.eq.2) then
            ic=ichar(c)
            if(ic.lt.z'30'.or.ic.gt.z'39') then
              if(c.eq.'*') then
                num=nbig
              else
                ierr=1
                goto 52
              endif
            else
              num=num*10+(ic-z'30')
            endif
            if(mode1.eq.0) then
              nmn=num
              nmx=num
              if(num.eq.nbig) then
                if(ifrq.eq.1) then
                  nmn=1
                else
                  nmn=0
                endif
              endif
            else
              nmx=num
            endif
          else if (mode.eq.3) then
            mln=1+mln
            if(mss.lt.0) mss=ip
          else
            ierr=1
            goto 52
          endif
        endif
      enddo

      ioplv=0
      kent=0
      do i=1,nswt
        icnt(i)=-1
      enddo

      if(lswnam(1).eq.0) then
        if(nmax(1).gt.0) then
          ioplv=1+ioplv
          iopn(ioplv)=1
          icnt(1)=0
        endif
      endif
      ip2=0
      
      cline=' '
      klen=0
      llen=len(cline)
      narg=iargf()
      do iarg=1,narg
        klen=klen+1
        ip1=klen
        call gtargf(iarg,cline(klen:llen))
        k=llen
        do while(cline(k:k).eq.' '.and.k.ge.klen)
          k=k-1
        enddo
        if(iarg.eq.1.and.cline(klen:k).eq.'!') then
          ierr=6
          goto 52
        endif
        ip2=k
        klen=k+1
        if(klen.gt.llen) pause 'chekcl: command line too long'
        ln=ip2-ip1+1
        swtest=cline(ip1:ip2)
        ibreak=0
   67   continue
        if(ibreak.ne.0) then
          swtest='-'//cline(ip1+ibreak:ip1+ibreak)
          ln=2


          ibreak=1+ibreak
        endif
        imatch=0
        iswt=0
        do i=1,nswt
          if(swtest(1:ln).eq.swnam(i)(1:ln)) then
            if(iswt.eq.0) iswt=i
            imatch=1+imatch
          endif
        enddo
        if(imatch.gt.1) then
          ierr=2
          goto 52
        else if(imatch.eq.1) then
          if(nmax(iswt).gt.0) then
            ioplv=1+ioplv
            iopn(ioplv)=iswt
          endif

          icnt(iswt)=0
        else




          if(cline(ip1:ip1).eq.'-'.and.ip2.gt.1+ip1.and.ibreak.eq.0) then
            ibreak=1
            goto 67
          endif


          if(ioplv.le.0) then
            ierr=1
            goto 52
          endif

          kent=kent+1
          iswt=iopn(ioplv)
          if(icnt(iswt).eq.0) then
            iptr(iswt)=kent
          endif

          ipt=iptr(iswt)
          ipt0=ipt
          iptl=kent
          do i=1,icnt(iswt)
            iptl=ipt
            ipt=itable(ipt,3)
          enddo

          if(ipt.ne.iptr(iswt)) pause 'circle does not close'
          itable(kent,1)=ip1
          itable(kent,2)=ip2
          itable(kent,3)=iptr(iswt)
          itable(iptl,3)=kent
          icnt(iswt)=1+icnt(iswt)
          if(icnt(iswt).eq.nmax(iswt)) then
            ioplv=ioplv-1
          endif
        endif
        if(ibreak.ne.0.and.ibreak.le.ip2-ip1) goto 67
      enddo

      klen=max0(0,klen-1)


   52 continue
      if(ierr.eq.0) then
        do i=1,nswt
          if(ifreq(i).ne.0.and.icnt(i).lt.0) ierr=3
          if(icnt(i).ge.0.and.(icnt(i).lt.nmin(i)
     1                        .or.icnt(i).gt.nmax(i))) ierr=4
        enddo
      endif
      if(ierr.ne.0) then
        write(0,"(72a1)") (cline(i:i),i=1,klen)
      else
        goto 99
      endif
      goto(1,2,3,4,5,6),ierr
      goto 99
    1 write(0,"('command line error')")
      goto 97
    2 write(0,"('non-unique abbreviation for identifier')")
      goto 97
    3 write(0,"('required identifier not present')")
      goto 97
    4 write(0,"('invalid number of arguments')")
      goto 97
    5 write(0,"('duplicate identifier in command description')")
    6 goto 97
   97 write(0,"('usage:')")
      do i=1,nswt
        if(ifreq(i).ne.0) then
          str1='required'
        else
          str1='optional'
        endif
        if(nmax(i).eq.nbig) then
          write(str2,"(i3,' or more arguments')") nmin(i)
        else if(nmax(i).eq.nmin(i).and.nmax(i).eq.1) then
          write(str2,"(i3,' argument         ')") nmax(i)
        else if(nmax(i).eq.nmin(i)) then
          write(str2,"(i3,' arguments        ')") nmax(i)
        else if(nmax(i).eq.1) then
          write(str2,"(i3,' -',i3,' argument    ')") nmin(i),nmax(i)
        else
          write(str2,"(i3,' -',i3,' arguments   ')") nmin(i),nmax(i)
        endif
        if(lswnam(i).eq.0) then
          str3=' '
        else
          str3=swnam(i)(1:lswnam(i))
        endif

        m10=mess(i)
        m1=m10

        mm=mlen(i)+m1-1
        do while(m1.lt.mm)
          m20=min0(m1+35,mm)
          m2=m20
          do while(m2.lt.mm.and.m2.gt.m1+20.and.string(m2:m2).ne.' ')
            m2=m2-1
          enddo

          if(string(m2:m2).ne.' ') m2=m20
          if(m1.eq.m10) then
            write(0,"(a10,1x,a21,1x,a8,': ',80a1)") str3,str2,str1
     1       ,(string(j:j),j=m1,m2)
          else




            write(0,"(43x,80a1)") (string(j:j),j=m1,m2)
          endif
          m1=m2+1
        enddo

      enddo

   99 continue
      return
      end
