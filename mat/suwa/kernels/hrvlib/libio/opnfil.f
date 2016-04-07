
c-------------------------------------------------------------------
c-------- Programmed by J. H. Woodhouse ----------------------------
c-------------------------------------------------------------------
      subroutine opnfil(lufl,namein,iap,ifile,irec,istat,lrec)
      character*(*) namein
      call opnflc(lufl,namein,iap,ifile,irec,istat,lrec,0)
      return
      end
