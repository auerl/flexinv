      subroutine savearow(io,elat,elon,slat,slon,edep,res,ttime,nn)
c
c... saves a row of the A matrix.
c    Input:   elat  ---  earthquake latitude
c    	      elon  ---  earthquake longitude
c    	      slat  ---  station longitude
c    	      slon  ---  station longitude
c    	      edep  ---  earthquake longitude
c    	      res   ---  residual
c    	      ttime ---  1-D travel time
c    Output:
c	      all the above
c 	      indarow --- index of nonzero element of a row of A matrix
c 	      arow    --- nonzero element of a row of A matrix
c
c
      implicit double precision (a-h,o-z)

c *** basic parameters, naming is obsolete
      parameter(maxfreq=8)
      parameter(maxnode=maxfreq*maxfreq*10+2)
      parameter(maxrknot=50)
      parameter(maxparm=maxnode*maxrknot)
      data pi /3.141592653579/
      data rad/57.2957795130824/
      data reprad /57.295779579/    

      common/amatrix/arow(maxparm),indarow(maxparm),nonzero_a
      real*4 a1(maxparm)
c      integer indarow(maxparm) ! alreadz defined in common - lapo 22.1.2010
      integer indtest(maxparm)

c... find nonzero elements
      nonzero_a=0
      j=1
      do i=1, nn
      	if(abs(arow(i)).gt.0.000001) then
		a1(j)=arow(i)
c		indarow(j)=i
		indtest(j)=i
		j=j+1
	endif
      enddo
      nonzero_a=j-1
c	print*, 'nn=',nn, nonzero_a, indtest(nonzero_a), a1(1)
      write(io) sngl(elat), sngl(elon), sngl(slat), sngl(slon),
     +		sngl(edep), sngl(res), sngl(ttime), nonzero_a
      write(io) (indtest(i),i=1,nonzero_a)
      write(io) (a1(i),i=1,nonzero_a)
      return
      end
