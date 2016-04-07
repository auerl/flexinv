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
