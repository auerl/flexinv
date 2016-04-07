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
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine backfl(lufl)
      include "openfile.h"
      if(jrecl(lufl).eq.0) then
        if(jrec(lufl).eq.0) then
          if(jfile(lufl).le.1) then
            call rewtp(lufl)
            jfile(lufl)=0
            jrec(lufl)=0
            return
          else
            call cmtio(jchn(lufl),2,2,ires,ierrno)
            if(ierrno.ne.0) call check('cmtio:backfile in backfl')
            call cmtio(jchn(lufl),1,1,ires,ierrno)
            if(ierrno.ne.0) call check('cmtio:forwardfile in backfl')
            jrec(lufl)=0
            jfile(lufl)=jfile(lufl)-1
          endif
        endif
      else 
        call rewfl(lufl)
      endif
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine backsp(lufl)
      include "openfile.h"
      if(jrecl(lufl).eq.0) then
        if(jrec(lufl).le.0) then
          call rewfl(lufl)
          jrec(lufl)=0
        else
          call cmtio(jchn(lufl),4,1,ires,ierrno)
          if(ierrno.ne.0) call check('cmtio:backrec in rewfl')
          jrec(lufl)=jrec(lufl)-1
        endif
      else
        jrec(lufl)=max0(0,jrec(lufl)-1)
      endif
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffi(lufl,ifbin,ibuf,nbytes,istat,nread,irec)
      include "openfile.h"
      dimension ibuf(1)
cTEST
c        print*,"in bffi, nbytes=",nbytes
      if(jrecl(lufl).eq.1.and.jfile(lufl).eq.200) then
        nread=0
        nget=nbytes
        irc=irec
        if(irc.eq.0) irc=jrec(lufl)+1
        iad=1
   10   nb=min0(nget,65532)
        call bffis(lufl,ifbin,ibuf(iad),nb,istat,nr,irc)
        if(istat.eq.2.and.nr.eq.nb) then
          if(nr.eq.nget) then
            nread=nread+nr
          else
            nr4=nr/4
            nr=4*nr4
            irc=irc+nr
            nread=nread+nr
            iad=iad+nr4
            nget=nget-nr
            goto 10
          endif
        else
          nread=nread+nr
        endif
      else
        call bffis(lufl,ifbin,ibuf,nbytes,istat,nread,irec)
      endif
cTEST
c        print*,"in bffi about to return with ibuf=",ibuf
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffin(lufl,ifbin,ibuf,nwords,istat,nread)
      dimension ibuf(1)
      nbytes=4*nwords
      call bffi(lufl,ifbin,ibuf,nbytes,istat,nread,0)
      nread=nread/4
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffis(lufl,ifbin,ibuf,nbytes,istat,nread,irec)
      include "openfile.h"
      dimension ibuf(1)
      krec=irec-1
      if(irec.eq.0) krec=jrec(lufl)
      if(jrecl(lufl).eq.0) then
        if(krec.lt.jrec(lufl)) then
	  write(6,"('before positioning tape -- backing up',i6)") 
     #     jrec(lufl)-krec
          call cmtio(jchn(lufl),4,jrec(lufl)-krec,ires,ierrno)
          write(6,"('positioning tape: krec=',i6,'  jrec=',i6)")
     1       krec,jrec(lufl)
        else if(krec.gt.jrec(lufl)) then
	  write(6,"('before positioning tape -- skipping',i6)")
     #     krec-jrec(lufl)
	  do iii=1,krec-jrec(lufl)
            call cread(jchn(lufl),ibuf,nbytes,nread,ierrno)
          enddo
c         call cmtio(jchn(lufl),3,krec-jrec(lufl),ires,ierrno)
          write(6,"('positioning tape: krec=',i6,'  jrec=',i6)")
     1       krec,jrec(lufl)
        endif
        call cread(jchn(lufl),ibuf,nbytes,nread,ierrno)
        if(nread.eq.255) nread=0
        istat=2
        if(ierrno.eq.0) then
          if(nread.eq.0) istat=3
        else
          istat=5
          write(6,"(7i10)") lufl,jchn(lufl),nbytes,nread,krec,nread,ierrno
          call check('cread in bffis 1')
        endif
        if(istat.eq.3) goto 30
        jrec(lufl)=1+krec
        return
   30   jfile(lufl)=1+jfile(lufl)
        jrec(lufl)=0
      else
      nb=nbytes
      if(jfile(lufl).ne.200) nb=min0(nb,jrecl(lufl))
cTEST
c        print*,"input of clseek:"
c        print*,jchn(lufl),jrecl(lufl)*krec,ires,ierrno
        call clseek(jchn(lufl),jrecl(lufl)*krec,0,ires,ierrno)
        if(ierrno.ne.0) call check('clseek in bffis') !here
        call cread(jchn(lufl),ibuf,nb,nread,ierrno)
        if(ierrno.ne.0) then
          write(6,"(7i10)") lufl,jfile(lufl),jrecl(lufl),jchn(lufl),nb,nread,ierrno
          call check('cread in bffis 2')
        endif
        istat=2
        if(nread.eq.0) istat=3
        if(istat.ne.3) then
          jrec(lufl)=krec+(nread+jrecl(lufl)-1)/jrecl(lufl)
        else
          jrec(lufl)=krec
        endif
      endif
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffo(lufl,ifbin,ibuf,nbytes,istat,irec)
      include "openfile.h"
      dimension ibuf(1)
cTEST
c        print*,"in bffo",nbytes,ibuf
      if(jrecl(lufl).eq.1.and.jfile(lufl).eq.200) then
        nread=0
        nget=nbytes
        irc=irec
        if(irc.eq.0) irc=jrec(lufl)+1
        iad=1
   10   nb=min0(nget,65532)
        call bffos(lufl,ifbin,ibuf(iad),nb,istat,irc)
        nget=nget-nb
        if(nget.ne.0) then
          nb4=nb/4
          irc=irc+nb
          iad=iad+nb4
          goto 10
        endif
      else
        call bffos(lufl,ifbin,ibuf,nbytes,istat,irec)
      endif
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffos(lufl,ifbin,ibuf,nbytes,istat,irec)
      include "openfile.h"
      dimension ibuf(1)
cTEST
c        print*,"in bffos, nbytes=",nbytes
      krec=irec-1
      if(irec.eq.0) krec=jrec(lufl)
      kchn=jchn(lufl)
c---------------------------------------------------------------
      if(jrecl(lufl).eq.0) then
c        write(6,"('krec,jrec',2i10)") krec,jrec(lufl)
        if(krec.lt.jrec(lufl)) then
          call cmtio(jchn(lufl),4,jrec(lufl)-krec,ires,ierrno)
        else if(krec.gt.jrec(lufl)) then
          call cmtio(jchn(lufl),3,krec-jrec(lufl),ires,ierrno)
        endif
        call cwrite(jchn(lufl),ibuf,nbytes,ires,ierrno)
        if(ierrno.ne.0.or.ires.ne.nbytes) then
          write(6,"('jchn,ierrno',2i10)") jchn(lufl),ierrno
          print *,'hello2'
          write(6,"('nbytes,ires',2i10)") nbytes,ires
          call check('cwrite in bffos 1')
        endif
        istat=2
        jrec(lufl)=1+krec
        return
c---------------------------------------------------------------
      else
        nb=nbytes
        if(jfile(lufl).ne.200) nb=min0(nbytes,jrecl(lufl))
cTEST
c        print*,"calling clseek",ires
        call clseek(jchn(lufl),jrecl(lufl)*krec,0,ires,ierrno)
        if(ierrno.ne.0) call check('clseek in bffos')
cTEST
c        print*,"calling cwrite",ires,nbytes,nb

        call cwrite(jchn(lufl),ibuf,nb,ires,ierrno)
cTEST
c        print*,"back from cwrite",ires,nbytes
        if(ierrno.ne.0.or.ires.ne.nbytes)then
cTEST
c            print*,"error status",ierrno
c            print *,'hello3'
            write(6,"('nbytes,ires',2i10)") nbytes,ires
            call check('cwrite in bffos 2')
        endif
        istat=2
        jrec(lufl)=krec+(nb+jrecl(lufl)-1)/jrecl(lufl)
        lenglu(lufl)=max0(lenglu(lufl),jrec(lufl))
      endif
c---------------------------------------------------------------
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine bffout(lufl,ifbin,ibuf,nwords,istat)
      dimension ibuf(1)
      nbytes=4*nwords
      call bffo(lufl,ifbin,ibuf,nbytes,istat,0)
      return
      end

c      character*80 function charmap(string,nchar)
c      integer*4 ia(3)
c      character*(*) string
c      if(string.eq.'TIME') then
c        call itime(ia)
c        write(charmap,"(i2,2(':',i2))") (ia(i),i=1,3)
c        nchar=8
c        do i=1,nchar
c          if(charmap(i:i).eq.' ') charmap(i:i)='0'
c        enddo
c      else if(string.eq.'DATE') then
c        call idate(ia)
c        write(charmap,"(i2,2('/',i2))") ia(2),ia(1),mod(ia(3),100)
c        nchar=8
c        do i=1,nchar
c          if(charmap(i:i).eq.' ') charmap(i:i)='0'
c        enddo
c      else
c        charmap='UNKNOWN'
c        nchar=7
c      endif
c      return
c      end
c----function charmap commented out by lapo for linux version march 2004
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine check(text)
      character*(*) text
cTEST
        print*,"calling cperror:"
      call cperror(text)
cTEST
        print*,"terminating"
      call exit(2)
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine chekcl(string)
      character*(*) string

      include "getgnl.h"
cc      character*80 cline
cc      character*30 swnam(30)
cc      dimension lswnam(30),ifreq(30),nmin(30),nmax(30)
cc     1         ,icnt(30),iopn(30),iptr(30),itable(100,3)

      include "getunx.h"
cc      common/ccmlin/swnam,cline
cc      common/icmlin/nswt,lswnam,ifreq,nmin,nmax
cc     1         ,icnt,iopn,iptr,itable

      external iargc,getarg

      call chekgl(string,iargc,getarg,ierr
     1  ,swnam,cline,deflt
     1  ,nswt,lswnam,ldeflt,ifreq,nmin,nmax,icnt,iopn,iptr,itable)
      if(ierr.ne.0) call exit(2)
      return
      end
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
            ic=ichar_(c)
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
c-------------------------------------------------------------------
c      subroutine chekli(string,ierr)
c      character*(*) string
c      include "getgnl.h"
c      include "getli.h"
c      external iargli,getali
c      call chekgl(string,iargli,getali,ierr
c     1  ,swnam,cline,deflt
c     1  ,nswt,lswnam,ldeflt,ifreq,nmin,nmax,icnt,iopn,iptr,itable)
c      end
c--routine above commented out by Lapo for linux version march 2004
      subroutine chprot(string)
      character*(*) string
      do i=1,len(string)
        string(i:i)=char(xor(ichar_(string(i:i)),z'80'))
      enddo

      return
      end
      subroutine closfl(lufl,istat)
      include "openfile.h"
      common/optdr/ifopt(nlu)

      if(jchn(lufl).le.0) return

      if(ifopt(lufl).ne.0) then
        call cmtio(jchn(lufl),5,1,ires,ierrno)
        if(ierrno.ne.0) call check('cmtio:rewind in closfl')
        ifopt(lufl)=0
      endif

      call cclose(jchn(lufl),ires,ierrno)
      if(ierrno.ne.0) then
         call cperror('cclose in closfl')
         print*,'logical unit ',lufl,' chanel ',jchn(lufl)
      endif

      jchn(lufl)=0
      return
      end
      DOUBLE PRECISION FUNCTION DFLDG(X)
      DOUBLE PRECISION X,Y,SMALL,D1,D2,D3,D4
      DATA SMALL/1.52587890625D-05/
      INTEGER*2 I(4),IMANT,IEXP,MINUS8,JJ(2),FIFTEEN
      EQUIVALENCE (Y,I(1)),(JJ(1),J)
      DATA IMANT/Z'00FF'/,IEXP/Z'7F00'/,MINUS8/-8/,JJ/0,0/
     1 ,FIFTEEN/15/
      Y=X                                        



      D1=AND(I(1),IMANT)
      JJ(2)=I(2)

      D2=J                              


      JJ(2)=I(3)
      D3=J                            


      JJ(2)=I(4)
      D4=J                           


      IPOW=ISHFT(AND(IEXP,I(1)),MINUS8)-66
      DFLDG=(16.D0**IPOW)*(D1+SMALL*(D2+SMALL*(D3+SMALL*D4)))
      IF(BTEST(I(1),FIFTEEN)) DFLDG=-DFLDG
      RETURN
      END
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine endfl(lufl)
      include "openfile.h"
      if(jrecl(lufl).ne.0) then
        leng=jrec(lufl)*jrecl(lufl)
        call clseek(jchn(lufl),leng,0,ires,ierrno)
        if(ierrno.ne.0.or.ires.ne.leng) call check('clseek in endfl')
        call ctrun(jchn(lufl),leng,ierrno)
        if(ierrno.ne.0) call check('ctrun in endfl')
      else
        call cmtio(jchn(lufl),0,1,ires,ierrno)
        if(ierrno.ne.0) call check('cmtio in endfl')
        jfile(lufl)=1+jfile(lufl)
        jrec(lufl)=0
      endif
      return
      end
      logical function error(ierrno)
      call cgterr(ierrno)
      if(ierrno.eq.0) then
        error=.FALSE.
      else
        error=.TRUE.
      endif
      return
      end
      subroutine fillu(lu,name)
      include "openfile.h"
      character*(*) name
      name=' '
      if(jchn(lu).gt.0) name=opnnam(lu)
      return
      end
      FUNCTION FLDG(X)
      INTEGER*4 IMANT,IEXP
      DATA IMANT/Z'00FFFFFF'/,IEXP/Z'7F000000'/
      LOGICAL BTEST
      EQUIVALENCE (Y,I)
      Y=X



      FLDG=FLOAT(AND(I,IMANT))
      IF(FLDG.NE.0.) THEN
        IPOW=ISHFT(AND(IEXP,I),-24)-70
        FLDG=(16.**IPOW)*FLDG
        IF(BTEST(I,31)) FLDG=-FLDG
      ENDIF
      RETURN
      END
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
c      subroutine fstat(lufl,lent,iftyp,isize)
      subroutine fstatf(lufl,lent,iftyp,isize)
      include "openfile.h"

cTEST
      print*,"call cfstat",jchn(lufl),isize,iftyp,ierrno

        call cfstat(jchn(lufl),isize,iftyp,ierrno)
cTEST
        print*,"okay back",ierrno

c        if(ierrno.ne.0) call check('cfstat on fstat')
        lent=jrecl(lufl)
        return
      end

c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine fstatn(name,lent,iftyp,isize)
      character*(*) name
      include "openfile.h"
        lu=-1

        call openfl(lu,name,1,0,0,istat,-1)
        call cfstat(jchn(lu),isize,iftyp,ierrno)
        lent=jrecl(lu)
        call closfl(lu,istat)
        return
      end
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
      subroutine getflpos(lufl,iposition,istatus)
      include "openfile.h"
      iposition=jrec(lufl)
      istatus=0
      return
      end

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
      subroutine getjukehost(string,ll)
      character*(*) string
      character*80 string2
      call getsafdir(string2,ll)
      open(99,file=string2(1:lnblnk(string2))//'/jukehost')
      read(99,"(a)") string
      close(99)
      ll=lnblnk(string)
      return
      end
c-------------------------------------------------------------------
      character*80 function getli(ident,nni,nbyts)
      character*(*) ident
      character*80 getgnl

      include "getgnl.h"

      include "getli.h"

      getli=getgnl(ident,nni,nbyts
     1  ,swnam,cline,deflt
     1  ,nswt,lswnam,ldeflt,ifreq,nmin,nmax
     1  ,icnt,iopn,iptr,itable)

      return
      end
      subroutine getolinkdir(string,ll)
      character*(*) string
      character*80 string2
      call getsafdir(string2,ll)
      string=string2(1:ll)//'/olinks'
      ll=ll+7
      return
      end
      subroutine getotapedir(string,ll)
      character*(*) string
      character*80 string2
      call getsafdir(string2,ll)
      string=string2(1:ll)//'/otapes'
      ll=ll+7
      return
      end
      subroutine getsafdir(string,ll)
      character*(*) string
      include "safdir.h"
      string=safdir
      ll=lnblnk(string)
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      character*80 function getunx(ident,nni,nbyts)
      character*(*) ident
      character*80 getgnl

      include "getgnl.h"
cc      character*80 cline
cc      character*30 swnam(30)
cc      dimension lswnam(30),ifreq(30),nmin(30),nmax(30)
cc     1         ,icnt(30),iopn(30),iptr(30),itable(100,3)

      include "getunx.h"
cc      common/ccmlin/swnam,cline
cc      common/icmlin/nswt,lswnam,ifreq,nmin,nmax
cc     1         ,icnt,iopn,iptr,itable

      getunx=getgnl(ident,nni,nbyts
     1  ,swnam,cline,deflt
     1  ,nswt,lswnam,ldeflt,ifreq,nmin,nmax
     1  ,icnt,iopn,iptr,itable)

      return
      end
      subroutine glinkd(linki,value,lvalue,istat)
      character*80 link,olinkdir
      character*(*) linki,value
c
      call getsafdir(link,ll)
      olinkdir=link(1:ll)//'/olinks'
c
      link=linki
      llink=istlen(link)
      ldir=istlen(olinkdir)
      value=link(1:llink)
      lvalue=llink
      istat=3
   10 continue
      call creadlink(olinkdir(1:ldir)//'/'//link(1:llink)//char(0)
     1   ,value,len(value),ires,ierrno)
      if(ierrno.eq.0) then
        link=value
        llink=ires
        lvalue=ires
        goto 10

      else if(ierrno.eq.2) then
        istat=1
      else if(ierrno.eq.22) then
        istat=0
      else
        call check('glink')
      endif

      return
      end
      subroutine glink(link,value,istat)
      character*(*) link,value
      llink=istlen(link)
      value=link
      istat=3
   10 call creadlink(link(1:llink)//char(0),value,len(value),ires,ierrno)
      write(6,"(2i5,':',80a1)") ires,ierrno,(value(i:i),i=1,ires)
      if(ierrno.eq.0) then
        link=value
        llink=ires
        goto 10
      else if(ierrno.eq.2) then
        istat=1
      else if(ierrno.eq.22) then
        istat=0
      else




        call check('glink')
      endif

      value(llink+1:len(value))=char(0)
      return
      end
      subroutine ilbyte(k,buf,ibyt)
      character*1 buf(*)
      k=ichar_(buf(ibyt+1))
      return
      end
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
	subroutine isbyte(k,ibuff,j)
c		needed in subroutine trnslt.
	dimension ibuff(*),mask(4)
 	data mask/z'00ffffff',z'ff00ffff',z'ffff00ff',z'ffffff00'/
	iw=1+j/4	
	ib=1+j-4*(iw-1)
   	ii=and(ibuff(iw),mask(ib))
 	ibuff(iw)=or(ii,ishft(k,8*(4-ib)))
 	return
	end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      function istlen(string)
      character*1 null
      parameter(null=char(0))
      character*(*) string
      k=len(string)
      do 10 i=1,k
      j=k+1-i
      if(string(j:j).eq.' '.or.string(j:j).eq.null) goto 10
      istlen=j
      goto 99
   10 continue
      istlen=0
   99 return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      function lenfl(lufl)
      include "openfile.h"
      if(jrecl(lufl).eq.0) then
        lenfl=0
      else
c        call fstat(lufl,idummy,idummy1,isize)
        call fstatf(lufl,idummy,idummy1,isize)
        lenfl=isize
      endif
      return
      end
      function lenlu(lufl)
      include "openfile.h"
      lenlu=lenglu(lufl)
      return
      end
      subroutine lockfl(lufl,iopt,isize,ierrno)
      include 'openfile.h'
cc     PARAMETER (NLU=40)
cc     PARAMETER (NLUUSE=20)
cc     COMMON/OPENFILE/JCHN(NLU),JREC(NLU),JFILE(NLU),JRECL(NLU)
cc    1  ,LENGLU(NLU)
cc     CHARACTER*32 OPNNAM
cc     COMMON/OPENNAME/OPNNAM(NLU)

      if(jchn(lufl).eq.0) then
        pause 'lock on file which is not open'
      else
        call clockf(jchn(lufl),iopt,isize,ires,ierrno)
      endif
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine lower(string)
      character*(*) string
      ll=len(string)
      do 10 i=1,ll
      ic=ichar_(string(i:i))
      if(ic.ge.65.and.ic.le.90) ic=ic+32
   10 string(i:i)=char(ic)
      return
      end
      function lufil(name)
      include "openfile.h"
      character*(*) name
      lufil=-1
      do i=1,nlu
        if(name.eq.opnnam(i).and.jchn(i).gt.0) then
          lufil=i
          goto 10
        endif
      enddo
   10 return
      end
      subroutine oascheckdrive(idrive,ierror)
      include 'oasstatus.h'
      character*80 string
      character*200 mess
c
      ierror=0
c
c---- first make sure that the drive is owned
c
      if(idrive.eq.1.or.idrive.eq.2) then
        if(iowned(idrive).eq.1) then
          write(string,"('user/default=',i1)") idrive-1
          call oassend(string,mess,lmess,-1)
          if(lmess.ne.0) then
            write(0,*) 'From OAS after user: ',mess(1:lmess)
            pause
          endif
c
          call oassend('status/system/type=9',mess,lmess,-1)
          if(lmess.ne.0) then
c              write(6,"('OAS tape status - ',a)") mess(1:lmess)
          endif
c
          string='status/volume/type=0'
          call oassend(string,mess,lmess,-1)
          if(lmess.ne.0) then
            if(mess(1:lmess).eq.volumemounted(idrive)(1:lmess).or.
     #           lnblnk(volumemounted(idrive)).eq.0) then
            else
              write(6,*) 'From OAS after volume inquiry: ',mess(1:lmess)
              ierror=1
              pause
            endif
          endif
          string='status/tape/type=0'
          call oassend(string,mess,lmess,-1)
          if(lmess.ne.0) then
            if(mess(1:lmess).eq.tapemounted(idrive)(1:lmess)) then
            else
              write(0,*) 'From OAS after tape inquiry  : ',mess(1:lmess)
              ierror=1
              pause
            endif
          endif
        endif
      endif
      return
      end
       
      subroutine oascloseserial(ierror)
      include "oasstatus.h"
      call cclose(ischan,ires,ierrno)
      ierror=ierrno
      return
      end
      subroutine oasflushserial(ierror)
      integer EWOULDBLOCK
      parameter (EWOULDBLOCK=35)
      include "oasstatus.h"
      character*1 cc
c
c---- read from serial line until no more characters
c---- first set descriptor for non-blocking i/o
c
      call cnoblock(ischan,1,ires,ierrno)
      ierrno=0
      do while(ierrno.ne.EWOULDBLOCK)
        call cread(ischan,cc,1,ires,ierrno)
      enddo
c
c---- clear for non-blocking i/o
c
      call cnoblock(ischan,0,ires,ierrno)
      end
      subroutine oasgetdrive(idrive,ierror)
      include "oasstatus.h"
      character*1 null
      character*200 mess
      character*80 string
      parameter (null=char(0))

      ierror = 0
      itry   = 0
c
c---- attempt to acquire the appropriate tape drive
c
  10  continue
      if(idrive.eq.1.or.idrive.eq.2) then
	write(6,"('trying to open',a)") nameblock(idrive)
        if(ioasverbose.ge.2) write(6,"('trying to open: ',a)")nameblock(idrive)
        call copen(nameblock(idrive)//null,itch,0,ierrno,0,0)
        if(ierrno.eq.0) then
          if(ioasverbose.ge.2)write(6,"('closing: ',a)")nameblock(idrive)
          call cclose(itch,ires,ierr)
          ierror=0
        else if( (ierrno.eq.6.or.ierrno.eq.16).and.ioaswait.eq.1) then
          if (itry.eq.0) then
            write(6,'(a,i1,a)') 'oasgetdrive: tape device # ',idrive,
     &                          ' in use: waiting ...'
            itry = 1
          endif
          call wait(4000)
          goto 10
        else if( (ierrno.eq.6.or.ierrno.eq.16).and.ioaswait.eq.0) then
          ierror = ierrno
        else if(ierrno.eq.5) then
          if(ioasverbose.ge.2) write(6,"('already off-line: ',a)")nameblock(idrive)
          ierror=0
        else
          write(6,"('oasgetdrive: ierrno',i5,1x,a)") ierrno,nameblock(idrive)
          ierror=ierrno
          call exit(2)
        endif
c
c---- if the tape was opened and closed (ierror=0) or off-line (ierror=5) we can
c---- now manipulate the optical drive
c
        if(ierror.eq.0) then
          iowned(idrive)=1
c
c---- set the default drive
c
          write(string,"('user/default=',i1)") idrive-1
          lstring=lnblnk(string)
          call oassend(string(1:lstring),mess,lmess,-1)
          if(lmess.ne.0) then
            write(0,*) 'From OAS after user: ',mess(1:lmess)
            pause
          endif
c
c---- go off-line and dismount the tape (if mounted)
c
          call oassend('offline',mess,lmess,-1)
          if(lmess.ne.0) then
            write(0,*) 'From OAS after offline: ',mess(1:lmess)
            pause
          endif
          call oassend('dismount',mess,lmess,-1)
          if(lmess.ne.0) then
            if(mess(1:lmess).eq.'E-OMSG-DISK NOT MOUNTED') then
            else if(mess(1:lmess).eq.'E-OMSG-TAPE NOT MOUNTED') then
            else
              write(0,*) 'From OAS after dismount: ',mess(1:lmess)
              pause
            endif
          endif
        else
          iowned(idrive)=0
        endif
        imounted(idrive)=0
      else
        write(6,"('illegal argument for oasgetdrive',i10)") idrive
        call exit(2)
      endif
      return
      end
      subroutine oasmountdrivec(idrive,volume,tape,irdwr,ierror)
      character*(*) volume
      character*(*) tape
      character*80 string
      character*200 mess
      character*1 quote
      parameter (quote='"')
      include 'oasstatus.h'
c
c---- check read/write flag
c
      ierror=0
      if(irdwr.ne.0.and.irdwr.ne.1) then
        ierror=1
        return
      endif
c
c---- check that tape name if provided
c
      lvol=lnblnk(volume)
      ltap=lnblnk(tape)
      if(ltap.le.0) then
        ierror=2
        return
      endif
c
c---- if this is a blank tape image, the volume must be given
c
      if(lvol.le.0.and.irdwr.eq.1) then
        ierror=5
        return
      endif
c
c---- first make sure that the drive is owned
c
      if(idrive.eq.1.or.idrive.eq.2) then
        if(iowned(idrive).eq.1) then
          
c
c---- set the default drive
c
          write(string,"('user/default=',i1)") idrive-1
          call oassend(string,mess,lmess,-1)
          if(lmess.ne.0) then
            write(6,*) 'From OAS after user: ',mess(1:lmess)
            pause
          endif
c
c---- go off-line and dismount the tape (if mounted)
c
          call oassend('offline',mess,lmess,-1)
          if(lmess.ne.0) then
            write(6,*) 'From OAS after offline: ',mess(1:lmess)
            pause
          endif
          call oassend('dismount',mess,lmess,-1)
          if(lmess.ne.0) then
            if(mess(1:lmess).eq.'E-OMSG-DISK NOT MOUNTED') then
            else if(mess(1:lmess).eq.'E-OMSG-TAPE NOT MOUNTED') then
            else
              write(0,*) 'From OAS after dismount: ',mess(1:lmess)
              pause
            endif
          endif
c
c---- mount the appropriate tape for reading
c
          if(irdwr.eq.0) then
            if(lvol.gt.0) then
              string='mount '//quote//volume(1:lvol)//quote//';'//quote//tape(1:ltap)//quote
            else
              string='mount '//quote//tape(1:ltap)//quote
            endif
            call oassend(string,mess,lmess,-1)
            if(lmess.ne.0) then
              write(6,*) 'From OAS after tape mount (read): ',mess(1:lmess)
              if(mess(1:lmess).eq.'E-OMSG-DISK IN USE') then
                ierror=3
                return
              else if(mess(1:lmess).eq.'E-DMGR-TAPE NOT FOUND') then
                ierror=4
                return
              else if(mess(1:lmess).eq.'E-JUKE-VOLUME NOT IN') then
                ierror=4
                return
              else if(mess(1:lmess).eq.'E-OMSG-NAME ALREADY EXISTS') then
                ierror=10
                return
              else
                pause
              endif
            endif
c
c---- mount the appropriate tape for writing
c
          else
            if(lvol.gt.0) then
              string='mount '//quote//volume(1:lvol)//quote//';'
              write(6,"(a)") string(1:lnblnk(string))
              call oassend(string,mess,lmess,-1)
              if(lmess.ne.0) then
                write(0,*) 'From OAS after volume mount (read): ',mess(1:lmess)
                if(mess(1:lmess).eq.'E-OMSG-DISK IN USE') then
                  ierror=3
                  return
                else if(mess(1:lmess).eq.'E-DMGR-TAPE NOT FOUND') then
                  ierror=4
                  return
                else if(mess(1:lmess).eq.'E-JUKE-VOLUME NOT IN') then
                  ierror=4
                  return
                else
                  pause
                endif
              endif
c
c---- test to see if compression screws things up
c
              string='mount/cblank '//quote//tape(1:ltap)//quote
c               string='mount/blank '//quote//tape(1:ltap)//quote
c
c              write(6,"(a)") string(1:lnblnk(string))
              call oassend(string,mess,lmess,-1)
              if(lmess.ne.0) then
                write(0,*) 'From OAS after tape mount (write): ',mess(1:lmess)
                if(mess(1:lmess).eq.'E-OMSG-DISK IN USE') then
                  ierror=3
                  return
                else if(mess(1:lmess).eq.'E-DMGR-TAPE NOT FOUND') then
                  ierror=4
                  return
                else if(mess(1:lmess).eq.'E-OMSG-NAME ALREADY EXISTS') then
                  ierror=10
                  return
                else
                  pause
                endif
              endif
            else
              ierror=11
              return
            endif
          endif
c
          imounted(idrive)=1
          volumemounted(idrive)=volume(1:lvol)
          tapemounted(idrive)=tape(1:ltap)
        else
          write(6,"('the optical drive is not owned by you:',i2)") idrive
          ierror=1
        endif
      else
        write(6,"('mountoasvolume: illegal drive',i2)") idrive
        call exit(2)
      endif
      return
      end
      subroutine oasmountdrive(idrive,volume,tape,irdwr,ierror)
      character*(*) volume
      character*(*) tape
      character*80 string
      character*200 mess
      character*1 quote
      parameter (quote='"')
      include 'oasstatus.h'
c
c---- check read/write flag
c
      ierror=0
      if(irdwr.ne.0.and.irdwr.ne.1) then
        ierror=1
        return
      endif
c
c---- check that tape name if provided
c
      lvol=lnblnk(volume)
      ltap=lnblnk(tape)
      if(ltap.le.0) then
        ierror=2
        return
      endif
c
c---- if this is a blank tape image, the volume must be given
c
      if(lvol.le.0.and.irdwr.eq.1) then
        ierror=5
        return
      endif
c
c---- first make sure that the drive is owned
c
      if(idrive.eq.1.or.idrive.eq.2) then
        if(iowned(idrive).eq.1) then
          
c
c---- set the default drive
c
          write(string,"('user/default=',i1)") idrive-1
          call oassend(string,mess,lmess,-1)
          if(lmess.ne.0) then
            write(6,*) 'From OAS after user: ',mess(1:lmess)
            pause
          endif
c
c---- go off-line and dismount the tape (if mounted)
c
          call oassend('offline',mess,lmess,-1)
          if(lmess.ne.0) then
            write(6,*) 'From OAS after offline: ',mess(1:lmess)
            pause
          endif
          call oassend('dismount',mess,lmess,-1)
          if(lmess.ne.0) then
            if(mess(1:lmess).eq.'E-OMSG-DISK NOT MOUNTED') then
            else if(mess(1:lmess).eq.'E-OMSG-TAPE NOT MOUNTED') then
            else
              write(0,*) 'From OAS after dismount: ',mess(1:lmess)
              pause
            endif
          endif
c
c---- mount the appropriate tape for reading
c
          if(irdwr.eq.0) then
            if(lvol.gt.0) then
              string='mount '//quote//volume(1:lvol)//quote//';'//quote//tape(1:ltap)//quote
            else
              string='mount '//quote//tape(1:ltap)//quote
            endif
            call oassend(string,mess,lmess,-1)
            if(lmess.ne.0) then
              write(6,*) 'From OAS after tape mount (read): ',mess(1:lmess)
              if(mess(1:lmess).eq.'E-OMSG-DISK IN USE') then
                ierror=3
                return
              else if(mess(1:lmess).eq.'E-DMGR-TAPE NOT FOUND') then
                ierror=4
                return
              else if(mess(1:lmess).eq.'E-JUKE-VOLUME NOT IN') then
                ierror=4
                return
              else if(mess(1:lmess).eq.'E-OMSG-NAME ALREADY EXISTS') then
                ierror=10
                return
              else
                pause
              endif
            endif
c
c---- mount the appropriate tape for writing
c
          else
            if(lvol.gt.0) then
              string='mount '//quote//volume(1:lvol)//quote//';'
              write(6,"(a)") string(1:lnblnk(string))
              call oassend(string,mess,lmess,-1)
              if(lmess.ne.0) then
                write(0,*) 'From OAS after volume mount (read): ',mess(1:lmess)
                if(mess(1:lmess).eq.'E-OMSG-DISK IN USE') then
                  ierror=3
                  return
                else if(mess(1:lmess).eq.'E-DMGR-TAPE NOT FOUND') then
                  ierror=4
                  return
                else if(mess(1:lmess).eq.'E-JUKE-VOLUME NOT IN') then
                  ierror=4
                  return
                else
                  pause
                endif
              endif
c
c---- test to see if compression screws things up
c
c              string='mount/cblank '//quote//tape(1:ltap)//quote
               string='mount/blank '//quote//tape(1:ltap)//quote
c
c              write(6,"(a)") string(1:lnblnk(string))
              call oassend(string,mess,lmess,-1)
              if(lmess.ne.0) then
                write(0,*) 'From OAS after tape mount (write): ',mess(1:lmess)
                if(mess(1:lmess).eq.'E-OMSG-DISK IN USE') then
                  ierror=3
                  return
                else if(mess(1:lmess).eq.'E-DMGR-TAPE NOT FOUND') then
                  ierror=4
                  return
                else if(mess(1:lmess).eq.'E-OMSG-NAME ALREADY EXISTS') then
                  ierror=10
                  return
                else
                  pause
                endif
              endif
            else
              ierror=11
              return
            endif
          endif
c
          imounted(idrive)=1
          volumemounted(idrive)=volume(1:lvol)
          tapemounted(idrive)=tape(1:ltap)
        else
          write(6,"('the optical drive is not owned by you:',i2)") idrive
          ierror=1
        endif
      else
        write(6,"('mountoasvolume: illegal drive',i2)") idrive
        call exit(2)
      endif
      return
      end
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
      subroutine oasopendrive(idrive,drive,ierror)
      character*80 string
      character*(*) drive
      character*200 mess
      character*1 quote,null
      parameter (quote='"')
c      parameter (null=ichar_(0))
c---changed by Lapo 11 Feb 04
      parameter(null="0")
      include 'oasstatus.h'
c
      ierror=0
c
c---- first make sure that the drive is owned
c
      if(idrive.eq.1.or.idrive.eq.2) then
        if(iowned(idrive).eq.1) then
          if(imounted(idrive).eq.1) then
c
c---- set the default drive
c
            write(string,"('user/default=',i1)") idrive-1
            call oassend(string,mess,lmess,-1)
            if(lmess.ne.0) then
                write(0,*) 'From OAS after user: ',mess(1:lmess)
            endif
c
c---- go on-line
c
            call oassend('online',mess,lmess,-1)
            if(lmess.ne.0) then
                write(0,*) 'From OAS after online: ',mess(1:lmess)
            endif
            drive=namedrive(idrive)
          endif
        endif
      endif
      return
      end
      subroutine oasopenserial(ierror)
      include "oasstatus.h"
      character*1 null
      parameter (null=char(0))
      character*200 mess
      character*30 dev
      character*80 safdir
      character*80 string
c
c---- get the name of the SAF directory and the serial port
c
      call getsafdir(safdir,lsafdir)
c
c      open(99,file=safdir(1:lsafdir)//'/jukeconfig',status='read')
c--changed Lapo 11 Feb 04
      open(99,file=safdir(1:lsafdir)//'/jukeconfig')
      read(99,"(i1,1x,a30)",iostat=ios) n,dev
      if(ios.ne.0) call exit(2)
      close(99)
      ldev=istlen(dev)
c
c---- try to open the serial line
c
      itry = 0
 10   call copen(dev(1:ldev)//null,ischan,2,ierrno,0,0)

      if(ierrno.ne.0) then
        if (ierrno.eq.16.and.ioaswait.eq.1) then
          if(itry.eq.0) then
             itry = 1
             write(*,"('serial port busy: waiting ...')")
           endif
           call wait(4000)
           goto 10
        endif

        ierror=ierrno
      else
        ierror=0
c
c---- make the serial line exclusive to this process
c
        call cgtflag(ischan,iflag,ires,ierrno)
        iflag=o'341'
        call cstflag(ischan,iflag,13,ires,ierrno)
        call cgtlmw(ischan,iword,ires,ierrno)
        call cstlmw(ischan,z'4000',ires,ierrno)
        call cstdtr(ischan,ires,ierrno)
        call cstexc(ischan,ires,ierrno)
c
c---- flush the serial port
c
        call oasflushserial()
c
c---- check for remaining errors
c
        do iuser=0,1
          write(string,"('user/default=',i1)") iuser
          lstring=lnblnk(string)
          call oassend(string(1:lstring),mess,lmess,-1)
          if(lmess.ne.0) then
              write(6,"('OAS - ',a)") mess(1:lmess)
          endif
          call oassend('status/system/type=9',mess,lmess,-1)
          if(lmess.ne.0) then
c              write(6,"('OAS tape status - ',a)") mess(1:lmess)
          endif
        enddo
  
c
c---- show the name legality check
c
        call oassend('show/name',mess,lmess,-1)
        if(lmess.ne.0) then
          if(ioasverbose.ge.2)  write(6,"('OAS - duplicate tape names allowed?',
     #               ' yes(1)/no(0): ',a)")mess(1:lmess)
        endif
c          if(mess(1:1).eq.'1') then
c            call oassend('set/name=0',mess,lmess,-1)
c            if(lmess.ne.0) then
c              write(6,"('OAS - ',a)")mess(1:lmess)
c            endif
c            call oassend('show/name',mess,lmess,-1)
c            if(lmess.ne.0) then
c              write(6,"('OAS - duplicate tape names allowed? yes(1)/no(0): ',a)")mess(1:lmess)
c            endif
c          endif
c        endif
      endif
      return
      end
      subroutine oassend(command,response,lr,iwt)
      character*(*) command,response
      parameter (maxoas=200)
      character*1 cr,nl
      parameter (cr=char(13))
      parameter (nl=char(10))
      character*(maxoas) msg
      character*(maxoas) temp
      logical condition
      include "oasstatus.h"
c
      ilc=0
      ip=0
      lc=lnblnk(command)
c
      do while (ilc.lt.lc)
        ilc=ilc+1
        if(command(ilc:ilc).eq.'^') then
          ip=ip+1
          ilc=ilc+1
          temp(ip:ip)=char(and(z'0f',ichar_(command(ilc:ilc))))
        else
          ip=ip+1
          temp(ip:ip)=command(ilc:ilc)
        endif
      enddo

      ip=1+ip
      if (ip.gt.79) then
         write(*,*) '**** WARNING: OASSEND: ip = ',ip
         write(*,'(a)') temp(1:ip)
      endif
      temp(ip:ip)=cr
      if(iwt.ge.0) then
        call cnoblock(ischan,1,ires,ierrno)
        if(ierrno.ne.0) call check('cnoblock in oassend')
      endif
      if(ioasverbose.gt.0) then
        write(6,"(a,a)") 'sending to OAS: ',temp(1:ip-1)
      endif
c
ccccccc trying to fix reading 3/22/93
      call cnoblock(ischan,1,ires,ierrno)
      iwrote=0
c
   22 continue
      call cwrite(ischan,temp,ip,ires,ierrno)
      iwrote=iwrote+1
      if(ierrno.ne.0) call check('cwrite in oassend')
cccccc removed the following line to not intefere with reading? 3/22/93
c      call cflush(ischan,1,ires,ierrno)
cccccc-----------------------------------------------------------------
      if(iwt.gt.0) call wait(iwt)
      lr=0

      if(ioasverbose.gt.0) then
        write(6,"(a)") 'begin reading from OAS: '
      endif
      ierrno=0
      condition=.true.
      do while(ierrno.eq.0.and.condition)
        lr=lr+1
        millis=0
   11   continue
	if(command(7:11).eq.'blank'.or.command(8:12).eq.'blank') then
	  if(millis.eq.0) write(6,"('searching for a blank area on disk')")
	  millis=10
        else
	  if(millis.gt.0.and.mod(millis,2000).eq.0) then
	    write(6,"('waiting',i3,' seconds for OAS response')") millis/1000
          endif
	endif
        if(millis.lt.30000) then
          call cread(ischan,msg(lr:lr),1,igot,ierrno)
          if(ierrno.ne.0) then
            call wait(10)
            millis=millis+10
            go to 11
          endif
        else
          write(6,"('giving up on reading from OAS:',i4)") lr
          write(6,"(a,a)") 'previously sent to OAS: ',temp(1:ip-1)
          if(iwrote.eq.1) go to 22
          condition=.false.
        endif
        condition=msg(lr:lr).ne.'>'
        if(.not.condition.and.lr.gt.3) then
          condition=(msg(lr:lr).eq.'>'.and.msg(lr-3:lr-3).eq.'<')
        endif
      enddo
      if(ioasverbose.gt.0) then
        write(6,"(a)") 'after reading from OAS: '
      endif

      if(ierrno.ne.0) lr=lr-1
      do i=1,lr
        if(msg(i:i).eq.cr) msg(i:i)=nl
        response(i:i)=msg(i:i)
      enddo

      if(msg(lr:lr).eq.nl) lr=lr-1
      if(response(1:1).eq.nl) then
        lr=lr-1
        do i=1,lr
          i1=i+1
          response(i:i)=response(i1:i1)
        enddo
      endif
cccccc trying to fix reading 3/22/93
      call cnoblock(ischan,0,ires,ierrno)
c

      if(iwt.ge.0) call cnoblock(ischan,0,ires,ierrno)
      if(lr.lt.1.or.response(lr:lr).ne.'>') then
        lr=-1
      else
        call oasstrip(response,lr)
      endif
      if (lr.gt.0.and.ioasverbose.eq.2)
     &     write(*,'(a,a)') 'recvd from OAS: ',response(1:lr)
      return
      end
      subroutine oasstrip(mess,lmess)
      character*1 null,nl
c      parameter (null=ichar_(0))
c      parameter (nl=ichar_(10))
c---changed by Lapo 11 Feb 04
      parameter(null="0")   ! I wonder how this is gonna work.
      parameter (nl="10")
      character*(*) mess
      do while (lmess.gt.0)
        if(mess(lmess:lmess).ne.'>'.and.mess(lmess:lmess).ne.nl
     1  .and.mess(lmess:lmess).ne.null.and.mess(lmess:lmess).ne.' ') goto 10
        lmess=lmess-1
      enddo
   10 continue
      return
      end
      subroutine oasverbose(iverb)
      include 'oasstatus.h'

      if (iverb.lt.0.or.iverb.gt.2) then
        write(*,*) 'invalid verbosity option = ',iverb
        write(*,*) ' 0 = no oas echoing'
        write(*,*) ' 1 = echo strings sent to oas'
        write(*,*) ' 2 = echo strings received from oas'
        return
      endif
      ioasverbose = iverb
      return
      end

     
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

     


c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine openfc(lufl,name,iap,ifile,irec,istat,lrec,inewi)
      character*(*) name
      
      call opnflc(lufl,name,iap,ifile,irec,istat,lrec,inewi)
      if(istat.ne.0) then
        ilen=istlen(name)
        write(0,*) 'openfc: error for file:',name(1:ilen)
        call exit(2)
      endif

      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine openfl(lufl,namein,iap,ifile,irec,istat,lrec)
      character*(*) namein
cTEST
c      print*,lufl,namein,iap,ifile,irec,istat,lrec
C	print*,"okay, in openfl, call openfil",istat
	print*,"okay, in openfl, call openfil",
     &lufl,namein,iap,ifile,irec,jstat,lrec

      istat=0
      call opnfil(lufl,namein,iap,ifile,irec,jstat,lrec)
      if(jstat.ne.0) then
	 write(6,'("file=",a60)') namein
         write(6,"('status',i10)") jstat
	 write(6,"('from openfl: file does not exist')")
	 call exit(1)
      endif
cTEST
	print*,"okay return from openfl",jstat

      return
      end
      subroutine openti(lufl,namea,iap,ifile,irec,istat)
      integer*4 namea(5)
      character*26 name
      write(name,"('/opt//',5a4)") (namea(i),i=1,5)
      call openfl(lufl,name,iap,ifile,irec,istat,0)
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine openw(lu,ifile,iac,i1,i2,istat,lrec)
      integer*4 ifile(20)
      character*80 file,file1
      write(file,1) (ifile(i),i=1,20)
    1 format(20a4)
      ip=1
      do while (file(ip:ip).ne.' '.and.ip.lt.80)
        ip=ip+1
      enddo
      call adnul(file(1:ip),file1)
      call openfl(lu,file1,iac,i1,i2,istat,lrec)
      return
      end

c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine opnfil(lufl,namein,iap,ifile,irec,istat,lrec)
      character*(*) namein
cTEST
c      print*,"opnfil:",lufl,namein,iap,ifile,irec,istat,lrec
	print*,"okay in opnfil, call opnflc",istat

      call opnflc(lufl,namein,iap,ifile,irec,istat,lrec,0)
cTEST
c	print*,"okay back from opnflc, now return from opnfil"


      return
      end
      subroutine opnflc(lufl,namein,iap,ifile,irec,istat,lrec,inewi)
c---- iap=1  => read only
c---- iap=2  => write only
c---- iap=4  => read/write
c
c---- inew=0 => old
c---- inew=1 => new
c---- inew=2 => fresh
c---- inew=3 => unknown
c
c---- istat=3  <= new file exists
c---- istat=4  <= old file does not exist
c
c only old (0) and new (1) are legal for tape images on the junkebox
c inew is ignored for mag tape
c
      character*1 quote,null
      parameter (quote='"')
      parameter (null=char(0))
      character*(*) namein
      character*80 path,name
      character*20 oasnam,value,oasvol
c
      logical nolink
c
      include "openfile.h"
      common/optdr/ifopt(nlu)
      include "oasstatus.h"
      save ifirst
      data ifirst/0/
c
cTEST
	print*,"okay in opnflc",iap,inew,istat

      if(ifirst.eq.0) then
        ifirst=1
        do i=1,nlu
          jchn(i)=0
        enddo
      endif
c
      inew=inewi
c
      if(lufl.gt.nluuse) then
	write(6,"('opnfil: lu out of range')")
	call exit(1)
      endif
c
c---- set iopt if dynamic block size
c
      iopt=0
      if(lrec.lt.0) iopt=1
      istat=0
c
c---- find an lu if calling program uses lu<=0
c
      if(lufl.le.0) then
        do i=nluuse+1,nlu
          if(jchn(i).eq.0) then
            lufl=i
            goto 71
          endif
        enddo
        write(6,"('opnfil: no available lu')")
	call exit(1)
   71   continue
      endif
c
c---- parse the file name
c
      opnnam(lufl)=namein
      path=namein
      name=namein

c
c---- check if this is an optical tape image---------- START OAS STUFF ---------
c
      if(path(1:5).eq.'/opt/'.or.path(1:5).eq.'/op1/'.or.path(1:5).eq.'/op2/') then
        if(path(1:5).eq.'/op2/') then
          idrive=2
        else
          idrive=1
        endif
c
        nolink=.false.
        name=path
        nbynam=len(name)
        do while (nbynam.gt.0.and.name(nbynam:nbynam).eq.' ')
          nbynam=nbynam-1
        enddo
c
c---- check if last character is '!' -- if so, don't link optical names
c
        if(name(nbynam:nbynam).eq.'!') then
          nolink=.true.
          nbynam=nbynam-1
        endif
c
c---- find the volume name, if any
c
        nbys=6
        do while (nbys.lt.nbynam.and.name(nbys:nbys).ne.'/')
          nbys=1+nbys
        enddo
        nbys=1+nbys
        oasnam=name(nbys:nbynam)
        loasnam=nbynam-nbys+1
        oasvol=name(6:nbys-2)
        loasvol=lnblnk(oasvol)
        write(6,"(a)") oasvol(1:loasvol)
        write(6,"(a)") oasnam(1:loasnam)
c
c---- initialize OAS
c
        itry = 0
  10    call oaswait(1)
        call oasopenserial(ierror)
	write(6,"('after oasopenserial',i6)") ierror
        if(ierror.ne.0) then
          istat=16
          return
        endif
c
c---- try to block the tape device
c
	write(6,"('trying to open ',a)") nameblock(idrive)
        call oaswait(0)
        call oasgetdrive(idrive,ierror)
        if (ierror.eq.6.or.ierror.eq.16) then
           if (itry.eq.0) then
             write(6,'(a,i1,a)') 'openflc: tape device # ',idrive,
     &                          ' in use: waiting ...'
             itry = 1
           endif
           call oascloseserial(ierror)
           call wait(4000)
           goto 10
        elseif (ierror.ne.0) then
           istat=ierror
           call oascloseserial(ierror)
           return
        endif
c
c---- open an existing tape image
c
        if(inew.eq.0) then
          if(nolink) then
            write(6,"('Using the exact OAS name (no links): ',a20)")oasnam
          else
            call glinkd(oasnam,value,lvalue,ilstat)
            if(ilstat.ne.1) then
              write(6,"( 'opnflc: unexpected status')")
              call exit(1)
            endif
            if(oasnam(1:loasnam).ne.value(1:lvalue)) then
              write(6,"('linked name: ',a,' --> ',a)") oasnam(1:loasnam),value(1:lvalue)
            endif
            oasnam=value
            loasnam=lvalue
          endif
c
          call oasmountdrive(idrive,oasvol(1:loasvol),oasnam(1:loasnam),0,ierror)
          if (ierror.ne.0) then
            call oascloseserial(ierror)
            istat=4
            return
          endif
c
c---- open a new tape image
c
        else if(inew.eq.1) then
          call oasmountdrive(idrive,oasvol(1:loasvol),oasnam(1:loasnam),1,ierror)
          if(ierror.ne.0) then
            istat=3
            call oascloseserial(ierror)
            return
          endif
        else 
          write(0,*) 'opnflc: illegal inew for an optical tape image'
          call exit(2)
        endif
        call oasopendrive(idrive,path,ierror)
        if(ierror.ne.0) then
          istat=7
          call oascloseserial(ierror)
          return
        endif
        ifopt(lufl)=1
c
c---- not an optical tape---------------------------------------------------------
c
      else 
       ifopt(lufl)=0
      endif
c
      if(path(1:5).eq.'/dev/') then
        inew=0
      endif
c
c---- opened for read only
c
      if(iap.eq.1.or.iap.eq.5) then
        isunop=0
        if(inew.ne.0) then
          write(0,*) 'opnflc: read only file not old'
          call exit(2)
        endif
c
c---- opened for write only
c
      else if(iap.eq.2) then
        isunop=1
c
c---- opened for read/write
c
      else
        isunop=2
      endif

      call adnul(path,path)
      if(inew.eq.0.or.inew.eq.1.or.inew.eq.3) then
        call copen(path,jchn(lufl),isunop,ierrno,0,0)
        if(ierrno.ne.0)  then
          write(6,"('after copen in opnflc - ierrno ',i5)") ierrno
          write(6,"('after copen in opnflc - ierrno ',a)") path(1:lnblnk(path))
          if(inew.eq.3.or.inew.eq.1) then
              inew=2
          else
              jchn(lufl)=0
              istat=4
              return
          endif
        else 
          if(inew.eq.1) then
            call cclose(jchn(lufl),ires,ierrno)
            jchn(lufl)=0
            istat=3
            return
          endif
          if(ifopt(lufl).eq.1) then
            call oascheckdrive(idrive,ierror)
            if(ierror.ne.0) then
              call oascloseserial(ierror)
              call closfl(lufl,ierror)
              istat=7
              return
            else
              call oascloseserial(ierror)
            endif
          endif
        endif
      else if(inew.eq.2) then
        call rcreao(path,jchn(lufl))
      endif
c

cTEST
      print*,"call fstatf", lufl,lent,isuntyp,isize

c      call fstat(lufl,lent,isuntyp,isize)
      call fstatf(lufl,lent,isuntyp,isize)

cTEST
      print*,"fstat okay",lufl,lent,isuntyp,isize


c      write(0,"('sun file mode (octal): ',8r,i10)") isuntyp
      if(isuntyp.eq.o'20666') then
        call cstexc(jchn(lufl),ires,ierrno)
        jrecl(lufl)=0
        jfile(lufl)=ifile
        jrec(lufl)=irec
        if(ifile.gt.0)  call cmtio(jchn(lufl),1,ifile,ires,ierrno)
        if(irec.gt.0)  call cmtio(jchn(lufl),3,irec,ires,ierrno)
      else
        jrecl(lufl)=iabs(lrec)
        jfile(lufl)=isuntyp
        if(iopt.eq.1) then
          jfile(lufl)=200
        endif
        jrec(lufl)=irec
        lenglu(lufl)=(isize+jrecl(lufl)-1)/jrecl(lufl)
        if(iap.eq.3.or.iap.eq.5.or.iap.eq.6) then
          jrec(lufl)=lenglu(lufl)+irec
        endif
        jrec(lufl)=max0(0,jrec(lufl))
      endif
      return
      end
      subroutine rcreao(path,ichan)
      character*1 null
      parameter (null=char(0))
      character*(*) path
      ichan=-1
      itry=0
      ifail=0
      ilen=istlen(path)
c read-write
      isunop=2
   30 ip2=ilen
      if(path(ip2:ip2).eq.'/'.and.ifail.ne.0) return
      call copen(path(1:ip2)//null,ichan,isunop,iopner,2,z'1a4')
      if(iopner.eq.0) then
c       call cclose(ichan,ires,ier)
        return
      else if(ifail.ne.0) then
        call check('copen in recreat failed')
      else if(iopner.eq.2) then
        itry=itry+1
        ifail=1
      else





        call check('copen in rcreat')
      endif

   10 continue
      do while(ip2.gt.0.and.path(ip2:ip2).ne.'/') 
        ip2=ip2-1
      enddo

   20 call cmkdir(path(1:ip2-1)//null,z'1ed',ires,ierrno)
      if(ierrno.eq.2) then
        itry=1+itry
        ip2=ip2-1
        goto 10
      else if(ierrno.eq.0) then
        itry=itry-1
        if(itry.gt.0) then
          ip2=1+ip2
          do while(path(ip2:ip2).ne.'/')
            ip2=1+ip2
          enddo

          goto 20
        else 




          goto 30
        endif

      else





        call check('rcreat')
      endif

      return
      end
      subroutine rcreat(path)
      character*(*) path
      call rcreao(path,ichan)
      if(ichan.gt.0) call cclose(ichan,ires,ier)
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine rewfl(lufl)
      include "openfile.h"
      if(jrecl(lufl).eq.0) then
        if(jfile(lufl).le.0) then
           call rewtp(lufl)
           jfile(lufl)=0
        else
          call cmtio(jchn(lufl),2,1,ires,ierrno)
          if(ierrno.ne.0) call check('cmtio:backfile in rewfl')
          call cmtio(jchn(lufl),1,1,ires,ierrno)
          if(ierrno.ne.0) call check('cmtio:forwardfile in rewfl')
        endif
      endif
      jrec(lufl)=0
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine rewtp(lufl)
      include "openfile.h"
      if(jrecl(lufl).eq.0) then
        call cmtio(jchn(lufl),5,1,ires,ierrno)
        if(ierrno.ne.0) call check('cmtio:rewind in rewtp')
        jrec(lufl)=0
        jfile(lufl)=0
      else
        call rewfl(lufl)
      endif
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine skipfl(lufl)
      include "openfile.h"
      if(jrecl(lufl).eq.0) then
        call cmtio(jchn(lufl),1,1,ires,ierrno)
        if(ierrno.ne.0) call check('cmtio:forwardfile in skipfl')
        jfile(lufl)=jfile(lufl)+1
      endif
      jrec(lufl)=0
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine upper(string)
      character*(*) string
      ll=len(string)
      do 10 i=1,ll
      ic=ichar_(string(i:i))
      if(ic.gt.95) ic=ic-32
   10 string(i:i)=char(ic)
      return
      end
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine wait(millis)
      isec=millis/1000
      msec=mod(millis,1000)
      nrep=msec/32
      mrem=mod(msec,32)
      if(isec.ne.0) call csleep(isec)
      do i=1,nrep
        call cusleep(32000)
      enddo
      if(mrem.ne.0) call cusleep(mrem*1000)
      return
      end
c--------following added by Lapo 11 Feb 04

      integer function ICHAR_(C)
      character C
      integer tmp
      tmp = ichar(C)
      if(tmp.ge.0) then
         ichar_=tmp
      else
         ichar_=256+tmp
      endif
      return
      end
