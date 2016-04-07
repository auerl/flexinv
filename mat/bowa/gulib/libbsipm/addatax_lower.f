      subroutine addatax_lower(ata,atd,ngpt,nstruc,xmod,iopt)
c
c Add or subtract a model from ATA and ATD.  The idea is to form
c ATA*X and add that to ATD.  This is also useful to save ATA*X
c for resolution tests.
c Note:  This is only for ATA saved in the lower-triangular form.
c Input:	ata,atd  ---  as such
c		ngpt     ---  number of horizontal knots
c		nstruc   ---  number of structural parameters (radial knots+2d topography)
c		xmod     ---  model to be added or subtracted
c		iopt     ---  1=add  
c			      2=subtract
c			      3=resolution test
c Output:
c		ata      ---  updated ATA matrix
c		atd	 ---  updated ATD matrix
c Y.G. 2002.
      real ata(1),atd(1),xmod(1)

      if(iopt.gt.3.or.iopt.lt.1) then
	 print*, 'Option for addatax_lower not reasonable, do nothing...'
	 return
      endif
      print*, 'adding option = ', iopt
      print*, xmod(1), xmod(5068), xmod(5069), xmod(10136)
      print*, xmod(10137), xmod(10498), xmod(10499), xmod(10860)
      ip1=0
      do ik1=1, nstruc
	do ilm1=1, ngpt
	   ip1=ip1+1
	   ip2=0
	   do ik2=1, nstruc
		do ilm2=1, ngpt
			ip2=ip2+1
			if(ip2.le.ip1) then
				iord=((ip1-1)*ip1)/2+ip2
			else
				iord=((ip2-1)*ip2)/2+ip1
			endif
			dot=ata(iord)*xmod(ip2)
			if(iopt.eq.1) then
			   atd(ip1)=atd(ip1)+dot
			else if(iopt.eq.2) then
			   atd(ip1)=atd(ip1)-dot
			else if(iopt.eq.3) then
			   atd(ip1)=dot
			endif
		enddo
           enddo
	enddo
      enddo
      return
      end
