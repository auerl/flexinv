      subroutine readrad(ifile,rknot,isplayer,nknot,numsplit)
c
c... read knot and layer information for radial B-spline basis
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


      dimension rknot(maxrknot)
      integer itemp(maxrknot)
      integer nknot,numsplit
      dimension isplayer(maxrknot,1)
      character*256 radfile
	
      write(*, "('enter file containing depth of spline nodes:')")
      read(*,"(a)") radfile
      call openfil_easy(ifile, radfile, 1, ierr)
      read(ifile, *) nknot, numsplit, (itemp(i), i=1, numsplit)
      k=0
      do i=1, nknot
	 read(ifile,*) depth
	 if(depth.gt.6371.d0.or.depth.lt.0.d0) stop 'unreasonable depth for spline nodes!!'
	 rknot(i)=6371.d0-depth
      enddo

      isplayer(1,1)=1      		! beginning radius of layer
      if(numsplit.eq.0) then      
      	 isplayer(1,2)=nknot		! ending radius of layer
      else
      	 do i=1, numsplit
	 	ii=itemp(i)
	 	isplayer(i,2)=ii
	 	isplayer(i+1,1)=ii+1
      	 enddo
      endif
      isplayer(numsplit+1,2)=nknot      	! top of the final layer
      write(*,"('number of splits =',i5)") numsplit
      do i=1, nknot
	 write(*,"(i8, f12.2)") i, rknot(i)
      enddo
      do i=1, numsplit+1
	 write(*, "('layer', i3, '   bottom index=', i3, '  top index=', i3)") i, 
     +				isplayer(i,1), isplayer(i,2)
      enddo
      close(ifile)
      return
      end

