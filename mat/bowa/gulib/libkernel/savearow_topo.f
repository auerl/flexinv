      subroutine savearow_topo(io,elat,elon,slat,slon,edep,
     &			res,ttime,nn,igrade,ifwate)
c
c... saves a row of the A matrix.
c    Input:   elat  ---  earthquake latitude
c    	      elon  ---  earthquake longitude
c    	      slat  ---  station longitude
c    	      slon  ---  station longitude
c    	      edep  ---  earthquake longitude
c    	      res   ---  residual
c    	      ttime ---  1-D travel time
c    	      igrade --- weighting factor*100, need to be devided by 100
c    	      ifwate --- option to weigh data set by the data coverage,
c			 this option is specifically designed for secondary
c			 data like SS-SdS or SdS1-SdS2.
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
c      integer indarow(maxparm) ! already defined in common - lapo 22.1.2010
      integer indnzero(maxparm)

c... find nonzero elements
      nonzero_a=0
      res1=res
      if(ifwate.ne.0) res1=res*real(igrade)*0.01
      j=1
      do i=1, nn
      	if(abs(arow(i)).gt.0.000001) then
		if(ifwate.ne.0) then
			a1(j)=arow(i)*real(igrade)*0.01
		else
			a1(j)=arow(i)
		endif
		indnzero(j)=i
		j=j+1
	endif
      enddo
      nonzero_a=j-1
      write(io) sngl(elat), sngl(elon), sngl(slat), sngl(slon),
     +		sngl(edep), sngl(res1), sngl(ttime), nonzero_a
      write(io) (indnzero(i),i=1,nonzero_a)
      write(io) (a1(i),i=1,nonzero_a)
cTEST
c      write(*,*) "datum info:",sngl(elat), sngl(elon), sngl(slat), sngl(slon),
c     +		sngl(edep), sngl(res1), sngl(ttime), nonzero_a
c	write(*,*) "matrix:"
c      write(*,*) (indnzero(i),i=1,nonzero_a)
c      write(*,*) (a1(i),i=1,nonzero_a)
c	print*,"............"

      return
      end
