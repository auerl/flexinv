      subroutine import(file1,file2,specf)
      character*(*) file1,file2,specf
      integer*4 STOLD,STUNK
      data STOLD/0/,STUNK/3/

      character*1 null
      parameter (null=char(0))
      character*1 spec(1000),line1(80),c
      character*80 line
      equivalence (line,line1)
      parameter (mxlev=20)
      dimension ist(mxlev),nrep(mxlev)

      parameter (mxtypes=5)
      parameter (mxleng=8)

      character*1 types(mxtypes)
      dimension lengs(mxtypes)
      data types/'a','i','j','f','d'/
c     program needs these to be even
      data lengs/ 2 , 4 , 2 , 4 , 8 /
      integer*2 jtemp(4)
      integer*4 itemp(2)
      real ftemp(2)
      double precision dtemp
      equivalence (itemp,jtemp,ftemp,dtemp)

      parameter (mxbyts=8192)

      integer*4 ibuf((mxbyts+mxleng)/4)
      integer*2 jbuf((mxbyts+mxleng)/2)
      equivalence (ibuf,jbuf)
c     write(6,*) ' import ',file1,' to ',file2,' using ',specf

      open(7,file=specf,status='old',access='read')

      ierr=0
      i=0




      do while(ierr.eq.0)
        read(7,"(a80)",iostat=ierr) line
        if(ierr.eq.0) then
          j=1



          do while (line1(j).ne.' '.and.line1(j).ne.null)
            i=i+1

            spec(i)=line1(j)
            j=j+1
          enddo
        endif
      enddo
      nlen=i
c     write(6,"(72a1)") (spec(i),i=1,nlen)

      luin=-1
      call openfc(luin,file1,1,0,0,istat,-1,STOLD)
      luout=-1
      call openfc(luout,file2,4,0,0,istat,-1,STUNK)
      numc=0
      lev=0

      i=1



      ibyt=0
      iwrd=0
      ngot=0
      igot=0
      nproc=0
      irem=0
      iwgot=0
      iwrd=0
  100 continue
      if(i.gt.nlen) then
	write(6,"('pattern exhaused')")
	call exit(1)
      endif
      c=spec(i)
      ic=ichar(c)
      if(ic.ge.z'30'.and.ic.le.z'39') then
        ival=ic-z'30'
        numc=10*numc+ival
        i=i+1


        goto 100
      else if(c.eq.'*') then
        numc=99999999
        i=i+1

        goto 100
      else if(c.eq.'(') then
        lev=lev+1
        if(lev.gt.mxlev) pause 'parens nested too deep'



        ist(lev)=i+1
        nrep(lev)=numc
        numc=0
        i=i+1
        goto 100
      else if(c.eq.')') then
        nrep(lev)=nrep(lev)-1
        numc=0
        if(nrep(lev).eq.0) then
          i=i+1

          lev=lev-1
          if(lev.lt.0) pause 'parens not properly paired'
        else




          i=ist(lev)
        endif
        goto 100
      else




        if(numc.eq.0) numc=1

        do k=1,mxtypes
          if(c.eq.types(k)) then
            ktyp=k
            ll=lengs(k)
            l=ll/2
            goto 80
          endif
        enddo

        pause 'type not found'
  80    continue
c       print*,'after 80 ktyp=',ktyp

        do ii=1,numc
   50     continue
c         print*,'ibyt+ll ',ibyt+ll,'  igot ',igot
          if(ibyt+ll.gt.igot) then
            if(ibyt.ne.0) call bffo(luout,1,ibuf,ibyt,j,0)
            do k=iwrd+1,iwgot
              jbuf(k-iwrd)=jbuf(k)
            enddo

            irem=iwgot-iwrd
            nword=1+(irem+1)/2
            call bffi(luin,1,ibuf(nword),mxbyts,j,mread,0)
c       print*,'after bffi ktyp=',ktyp,'  mread=',mread
            ngot=ngot+mread
c           print*,ngot,' bytes read'
            if(j.eq.3.and.mread.ne.0) pause 'unexpected mread'
            if(j.ne.2.and.j.ne.3) pause 'unexpected status from bffi'
            ishf=(nword-1)*2-irem

            if(ishf.ne.1.and.ishf.ne.0) pause 'unexpected ishf'
            if(ishf.ne.0) then
              kof=irem+ishf
              do k=1,mread/2
                jbuf(k+irem)=jbuf(k+kof)
              enddo
            endif

            iwgot=mread/2+irem
            igot=2*iwgot
            ibyt=0
            iwrd=0
            if(mread.eq.0) goto 99
            goto 50
          endif

          do k=1,l
            jtemp(k)=jbuf(iwrd+k)
          enddo

          goto(1,2,3,4,5),ktyp
          pause 'ktyp out of range'
    1     goto 200
    2     goto 200
    3     goto 200
    4     ftemp(1)=fldg(itemp(1))
          goto 200
    5     dtemp=dfldg(itemp(1))
          goto 200
  200     continue
          do k=1,l
            jbuf(iwrd+k)=jtemp(k)
          enddo

          ibyt=ibyt+ll
          iwrd=iwrd+l

          nproc=nproc+ll
        enddo

        numc=0
        i=i+1

        goto 100
      endif
   99 if(igot.ne.0) then
        call bffo(luout,1,ibuf,igot,j,0)
        pause 'incomplete set at end of file -- not processed'
      endif

c     write(6,"('bytes read',i8,'  bytes pocessed'i8)") ngot,nproc
      close(7)
      call closfl(luin,ista)
      call closfl(luout,ista)
      return
      end
