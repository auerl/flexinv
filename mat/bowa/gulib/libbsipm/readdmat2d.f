      subroutine readdmat2d(io,ngpt,dmat)
c
c read 2-D damping matrix 
c Input:  io   --- file id
c        ngpt  --- number of expected horizontal elements
c Output:
c	 dmat  --- damping matrix
c
c
      dimension dmat(1)
      character*256 dampfile	
      write(*, "('enter the 2D damping matrix file:')")
      read(*, "(a)") dampfile
      open(io, file=dampfile, status='old')
      idamp = 1
5     read(io,*, end=11) dmat(idamp)
      idamp = idamp + 1
      goto 5
11    continue
      idamp=idamp-1
      print*, '***damping parameters read:  ', idamp
      iexpected=ngpt*(ngpt+1)/2
      if(idamp.ne.iexpected) then
	print*, 'problem with damping matrix file:'
	print*, 'number of damping elements not consistent!'
	stop
      endif
      close(io)
      return
      end     
     
