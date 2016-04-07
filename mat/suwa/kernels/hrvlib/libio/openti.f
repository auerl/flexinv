      subroutine openti(lufl,namea,iap,ifile,irec,istat)
      integer*4 namea(5)
      character*26 name
      write(name,"('/opt//',5a4)") (namea(i),i=1,5)
      call openfl(lufl,name,iap,ifile,irec,istat,0)
      return
      end
