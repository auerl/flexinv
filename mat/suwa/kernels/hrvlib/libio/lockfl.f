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
