      subroutine openfil_easy(ifile, filename, iopt, ierr)
      character*(*) filename
      integer ifile, iopt, ierr
c     iopt: 	1 --- simple read
c     iopt: 	2 --- simple write
c     iopt: 	3 --- append
c     iopt: 	4 --- write unformatted 
c
c
      ierr = 0
      if(iopt.eq.1) then
         open(ifile, name=filename(1:lnblnk(filename)), iostat=ios)
	 if(ios.ne.0) then
		ierr = 1
	 endif
      else if(iopt.eq.2) then
         open(ifile, name=filename(1:lnblnk(filename)), status='unknown')
      else if(iopt.eq.3) then
         open(ifile, name=filename(1:lnblnk(filename)), access="append")	 				   
      else if(iopt.eq.4) then
c unformatted option has been disabled due to operating system problem
         open(ifile, name=filename(1:lnblnk(filename)), 
     &			status='unknown')	 
      else
	 ierr = 1
      endif
      if(ierr.ne.0) print*, 'Error opening file ', filename(1:lnblnk(filename))
      return
      end
