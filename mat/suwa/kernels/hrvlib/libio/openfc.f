

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
