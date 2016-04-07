	subroutine getataf_lower(ata,numatd,ataf,mxparm)
c... fill out a full ata matrix from a lower triangular matrix
c... for anisotropic inversion
c  Input:	ata  ---   ata in lower-triangular form
c  		mxparm ---  max dimension of ataf
c  		numatd ---  number of atd elements
c  Output:
c		ataf ---   full ata matrix
c
	dimension ata(1),ataf(mxparm,mxparm)

	if(mxparm.lt.numatd)then
	   print *,"dimension of ataf in subroutine get_ataf"
	   print *,"exceeds the limit"
	   call exit(1)
	endif
c
c	fill lower triangle
c
	ind=numatd*(numatd+1)/2

	do j=numatd,1,-1
	   do i=j,1,-1
              ataf(i,j)=ata(ind)
              ind=ind-1
	   enddo
	enddo
c
c	fill upper triangle
c
	do i=1,numatd
	   do j=1,i-1
	      ataf(i,j)=ataf(j,i)
	   enddo
	enddo


c	k=1
c	do i=1,numatd
c	   do j=1,i
c	      ataf(i,j)=ata(k)
c	      ataf(j,i)=ata(k)
c	      k=k+1
c	   enddo
c	enddo

	return
	end
