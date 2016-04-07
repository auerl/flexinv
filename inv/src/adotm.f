c====================================================================
c
c // PROGRAM ADOTM, Ludwig Auer 2012
c
c computes synthetic data by computing the A-matrix with a certain
c model in a chosen parameterization. Doesn't do any checks so be
c careful to chose the correct model, assumes that vals are in %/100
c
c original: from Lapo Boschi, 1995-2009
c
c compile with ifort -i8 -132 -shared-intel -O3 -mcmodel=large adotm.f -o adotm
c
c====================================================================



	character*150 matfile,indfile,pntfile,modfile,outfile

	parameter(pi=3.141592653589)
        parameter(m=25000000)
        parameter(nonz=1500000000)
	parameter(nmax=2500000)

	dimension x(nmax)
        dimension indx(nonz),mpoin(0:m)
        dimension values(nonz),t(m),t0vec(m)

	character chnfree*6
	real ran2
        integer*8 jj

	iseed=-1
	xtmp= ran2(iseed)

 
        print*,"------- SYNTHETICS FOR NEXT SUBMATRIX -------"
        
 
c        print*,"matrix file?"
        read(*,*)matfile
c        print*,"index file?"
        read(*,*)indfile
c        print*,"pointer file?"
        read(*,*)pntfile
c        print*,"input model file?"
        read(*,*)modfile
c        print*,"uncertainty for gen data <0 read from sigma.dat for each data"
	read*,sigma
c	if(sigma.lt.0)then
c	   readsfile=1
c	   open(99,file='sigma.dat',status="old")
c	else
c	   readsfile=0
c	endif
	
c	print*,"file for synthetic database?"
	read*,outfile
c	print*,"number of free parameters?"
	read(*,*)nfree
c        print*,"number of observations?"
        read(*,*)ndata


c open and read model file, values expected in %

	open(1,file=modfile,status="old")
	do k=1,nfree
	   read(1,*)kin,x(k)
c	   x(k)=x(k)/100.
	enddo
	close(1)

c open and read A matrix files
	print*,"Reading the matrix..."
	relwei=1.
	icol=1
	jj=1
	mpoin(0)=0

c readmatrix2 is faster than lapos original one      
        call readmatrix2(ndata,matfile,indfile,pntfile,
     &  mpoin,indx,values,icol,jj,nonz,nfree,m,relwei,
     &  nmax)


        print*," ... Matrix read!"


	open(79,file=outfile)
c dot A with input model 
	print*,"A dot m, data:",icol-1
        do l=1,icol-1
           tot=0.
           do ll=mpoin(l-1)+1,mpoin(l)
              tot=tot+(values(ll)*x(indx(ll)) )
c              print*,(values(ll)*x(indx(ll)) ),tot
           enddo
           if(mod(l,1000).eq.0) print*,"Datum #",l
c now add random noise to predicted datum. Box-Muller transformation.
c         y1 = sqrt( - 2 ln(x1) ) cos( 2 pi x2 )
c         y2 = sqrt( - 2 ln(x1) ) sin( 2 pi x2 )
	   x1ran=ran2(iseed)
	   x2ran=ran2(iseed)
	   y1ran=sqrt(-2.*log(x1ran))*cos(2.*pi*x2ran)
	   if(readsfile.eq.1)then
	      read(99,*)sigma
	   endif
	   y1ran = y1ran * sigma ! normal distribution with mean=0 and stand. dev. = sigma
	   tot = tot + y1ran           
	   write(79,*)tot
        enddo

	close(79)
	if(readsfile.eq.1)then
	   close(99)
	endif
	end
c==============================================================



c==============================================================
	subroutine readmatrix2(ndata,namexxx,nameind,namepoi,
     &  mpoin,indx,values,icol,jj,nonz,n,m,relwei,nmax)


        integer*8 mpoin(0:m),nonz,jj
	integer*8, dimension(:) :: indx(nonz)
	integer*4, dimension(:) :: indx0(2000000000)
                                         
	dimension values(nonz),values0(2000000000)
	character*150 namexxx,nameind,namepoi,namerhs,wdir

	open(3,file=namepoi,status='old')
	print*,"icol,jj,n: ",icol,jj,n

	icol0=icol
	nrec0=1
	do icol=icol0,icol0+ndata-1 ! loop over pointer values
	   read(3,*)mptemp 
	   mpoin(icol)=mptemp+mpoin(icol0-1) 
	enddo
	close(3)

	open(1,file=namexxx,status='old',access='direct',form='unformatted',recl=4*mptemp)
	open(4,file=nameind,status='old',access='direct',form='unformatted',recl=4*mptemp)
	read(4,rec=1) (indx0(iii), iii=1,mptemp) ! index in row = index of corresponding pixel
	read(1,rec=1) (values0(iii), iii=1,mptemp) ! values: raypath through this pixel

c        ! apply weighting and store weighted cummulative sensitivities
	do iii=jj,jj+mptemp-1
	   indx(iii)=indx0(iii-jj+1)
	   indx0(iii-jj+1)=0
	   if(indx(iii).gt.n)then
	      print *,"ERROR: undefined voxel index", n, indx(iii)
	      stop
	   endif
	   values(iii)=values0(iii-jj+1)
           values0(iii-jj+1)=0.

	end do

	jj=jj+mptemp 
	close(4)
	close(1)
	close(3)
	return
	end
c==============================================================






c==============================================================
	FUNCTION ran2(idum)
	INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
	REAL ran2,AM,EPS,RNMX
	PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     &   IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     &   NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
	INTEGER idum2,j,k,iv(NTAB),iy
	SAVE iv,iy,idum2
	DATA idum2/123456789/, iv/NTAB*0/, iy/0/
	if (idum.le.0) then
	   idum=max(-idum,1)
	   idum2=idum
	   do 11 j=NTAB+8,1,-1
	      k=idum/IQ1
	      idum=IA1*(idum-k*IQ1)-k*IR1
	      if (idum.lt.0) idum=idum+IM1
	      if (j.le.NTAB) iv(j)=idum
 11	   continue
	   iy=iv(1)
	endif
	k=idum/IQ1
	idum=IA1*(idum-k*IQ1)-k*IR1
	if (idum.lt.0) idum=idum+IM1
	k=idum2/IQ2
	idum2=IA2*(idum2-k*IQ2)-k*IR2
	if (idum2.lt.0) idum2=idum2+IM2
	j=1+iy/NDIV
	iy=iv(j)-idum2
	iv(j)=idum
	if(iy.lt.1)iy=iy+IMM1
	ran2=min(AM*iy,RNMX)
	return
	END
c==============================================================

