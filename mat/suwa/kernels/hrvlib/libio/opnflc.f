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
      call fstat(lufl,lent,isuntyp,isize)
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
