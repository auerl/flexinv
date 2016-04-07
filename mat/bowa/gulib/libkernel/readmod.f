      subroutine readmod
c*** read in the model
c    coef(4,8,20) contains the 4 polynomial coefficients for
c    each of 8 parameters for up to 20 layers
c    assume that first layer is inner core and second layer is outer core
c    note: layer number starts from surface to core (opposite to
c    the order that it is read in).
c    edited by Yu Gu, 1998.
      implicit real*8(a-h,o-z)
      character*256 filnam,modnam
      common/coeff/coef(4,8,20)
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)
      data pi/3.14159265358979d0/,rtol/0.d0/
      integer moho

      rad=180.d0/pi
      write(*,"('enter model file name :')")
      read(*,'(a256)') filnam
      open(8,file=filnam)
      read(8,'(a256)') modnam
c ifanis = 1 (anisotropic)
      read(8,*) ifanis,tref,ifdeck
c
c numlyr=number of layers, nic = number of sublayers in Inner Core
c noc = number of sublayers in Outer Core, rnorm is normalization radius
c
      read(8,*) numlyr,nic,noc,rnorm,moho
      npar=5
      if(ifanis.ne.0) npar=8
      do 5 i=1,numlyr
        k=numlyr-i+1	! reversing the indexing for the model
        read(8,*) nl(k),xb(k),xt(k)
c        read(8,905) nl(k),junk,xb(k),xt(k)
c        write(99,805) nl(k),xb(k),xt(k)
c  805 format(i5,2f10.3)
c  905   format(2i5,2d15.9)
c        dum=0.d0
        do j=1,npar
          read(8,810)(coef(jj,j,k),jj=1,4)
c*** this assumes polynomials are in big Q and are constant only
          if(j.eq.4.or.j.eq.5.and.coef(1,j,k).ne.0.d0)
     +        coef(1,j,k)=1.d0/coef(1,j,k)
c          write(99,820)(coef(jj,j,k),jj=1,4),dum
  810 format(5g16.9)
c  820 format(5g16.9)
c  810 format(5g9.5)
c  910   format(4d16.9)
        enddo
c*** i'm not sure if i really need to do this
        if(npar.eq.5) then
          do kk=1,4
            coef(kk,6,k)=coef(kk,2,k)
            coef(kk,7,k)=coef(kk,3,k)
            coef(kk,8,k)=0.d0
          enddo
          coef(1,8,k)=1.d0
        end if
c the line below finds the layer for 670.
c added 220 and 400 by J. Gu.
        if(abs(xt(k)-5700.).lt.20.) llm=k
    5   continue
c**** we are assuming a prem-type parameterization for the derivs
c  ie, cubic polynomials so k=1,4 and we have 5 parameters (vph,
c  vpv,vsh,vsv,eta) so we have 20 coefficients per layer
      nplay=20
      close(8)
c  find slownesses just below and above icb,cmb, and 660
c  here I have added 400 and 220.      
      lic=numlyr
      loc=numlyr-1
      iq=numlyr
c y is the normalized radius, getmod evaluate the velocities
c and eta from cubic polynomial coef at top of each layer.
      y=(xt(iq)-rtol)/rnorm
 	print*, 'y =', y*rnorm
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
      picblp=y*rnorm/vpv/rad
      picbls=y*rnorm/vsv/rad
      iq=numlyr-1
      y=(xb(iq)+rtol)/rnorm
 	print*, 'y =', y*rnorm
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
      picbup=y*rnorm/vpv/rad
      y=(xt(iq)-rtol)/rnorm
 	print*, 'y =', y*rnorm
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
      pcmblp=y*rnorm/vpv/rad
      iq=numlyr-2
      y=(xb(iq)+rtol)/rnorm
 	print*, 'y =', y*rnorm
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
      pcmbup=y*rnorm/vpv/rad
      pcmbus=y*rnorm/vsv/rad
c looks like there is no slowness for 660 upper? Yu, 1998.
      iq=llm
      y=(xt(iq)-rtol)/rnorm
      call getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
      p660lp=y*rnorm/vpv/rad
      p660ls=y*rnorm/vsv/rad
	print*, 'picbup=', picbup, '  picblp=', picblp
	print*, 'picbup=', picbup, '  picblp=', picblp
	print*, 'pcmbup=', pcmbup, '  pcmblp=', pcmblp
      return
      end
