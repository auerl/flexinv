      subroutine find_rsplit(nknot,rknot,isplayer,nsplit)
c
c finds the spliting indices for the radial knots.  Note: this sub 
c does not care about the ordering of the spline nodes (CMB to MOHO
c or  MOHO to CMB)
c  Input: 	nknot  ----  number of radial basis functions
c  	 	rknot  ----  radius for the radial basis
c  Output:
c             isplayer ----  the indices for all the spliting depths
c               nsplit ----  number of splits
c
      dimension  rknot(1)
      dimension  isplayer(1)
		
      nsplit=0
      do i=1, nknot-1
	 if(abs(rknot(i+1)-rknot(i)).lt.0.1) then
		nsplit=nsplit+1
		isplayer(nsplit)=i
	 endif
      enddo
      write(*,"('number of spliting depths =',i4)") nsplit		
      return
      end
       
