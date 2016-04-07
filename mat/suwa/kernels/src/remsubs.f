      subroutine rdearmod(ilu,icm,ici,nlyr,r4norm)
c-----reads an earth model and fills the blocks /earmod/
c
      implicit double precision (a-h,o-z)
      logical isotrp
      common/earmod/ rnorm,coef(4,8,20),xb(20),xt(20),numlyr,nic,
     # noc,moho,nl(20)
c
c-----imod = 1 ; model is anisotropic
c
      imod = 1
      isotrp = .false.
      read(ilu,1) numlyr,nic,noc,rnorm,moho
    1 format(3i5,d15.9,i5)
      do 5 k=1,numlyr
      read(ilu,2)nl(k),junk,xb(k),xt(k)
    2 format(2i5,2d15.9)
      do 3 j=1,8
    3 read(ilu,4)(coef(jj,j,k),jj=1,4)
    4 format(4d16.9)
    5 continue
      close(ilu)
      nlyr=numlyr
      r4norm=(rnorm)
c-----find mantle and core layers
      do  10 k=1,numlyr
      if(coef(1,3,k).lt.1.d-8.and.coef(1,3,k+1).gt.1.d-8) icm=k
      if(coef(1,3,k).gt.1.d-8.and.coef(1,3,k+1).lt.1.d-8) ici=k
   10 continue
      return
      end
c
      double precision function cubic(x,c)
      implicit double precision (a-h,o-z)
      dimension c(4)
      cubic=c(1)+x*(c(2)+x*(c(3)+x*c(4)))
      return
      end

      subroutine evanisomod(ras8,rho,vpv,vsv,q2,q1,vph,vsh,eta)
      implicit double precision (a-h,o-z)
      common/earmod/ rnorm,coef(4,8,20),xb(20),xt(20),numlyr,nic,
     # noc,moho,nl(20)
      do 11 i=1,numlyr
         if(ras8.lt.xt(i)) go to 12
   11 continue
   12 iq=i
c-----normalize to unit earth radius
      y=ras8/rnorm
      rho=cubic(y,coef(1,1,iq))
      vpv=cubic(y,coef(1,2,iq))
      vsv=cubic(y,coef(1,3,iq))
      q1=cubic(y,coef(1,4,iq))
      q2=cubic(y,coef(1,5,iq))
      vph=cubic(y,coef(1,6,iq))
      vsh=cubic(y,coef(1,7,iq))
      eta=cubic(y,coef(1,8,iq))
      return
      end
c
