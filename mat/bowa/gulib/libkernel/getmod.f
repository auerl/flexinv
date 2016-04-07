      subroutine getmod(y,iq,rho,vpv,vsv,vph,vsh,eta)
c*** returns model vector at normalized radius y which
c*** shell iq. When iq is changed, the routine checks
c*** to see if the layer is isotropic (returned through 
c*** common isot). 
      implicit real*8(a-h,o-z)
      logical isotrp
      dimension test(4)
      common/coeff/coef(4,8,20)
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/isot/isotrp

      data iqs/-1/,test/1.d0,0.d0,0.d0,0.d0/
      if(iq.ne.iqs) then
        iqs=iq
        isotrp=.false.

        if(iso.ne.1) then
          sum=0.d0
          do i=1,4
            sum=sum+coef(i,2,iq)-coef(i,6,iq)
     +        +coef(i,3,iq)-coef(i,7,iq)+coef(i,8,iq)-test(i)
          enddo

          if(sum.eq.0.d0) isotrp=.true.
        end if
      end if

cl--test * *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c	  print*,"la: voxint: radcontrib: calcqvec: qtau: getmod: CRASH!"
c	  print*,"la: voxint: radcontrib: calcqvec: qtau: getmod: iq (max:20) = ",iq

        if(iq.gt.20) then
	  print*,"IQ IS TOO LARGE: ",iq

      end if

      rho=coef(1,1,iq)+y*(coef(2,1,iq)+y*(coef(3,1,iq)+y*coef(4,1,iq)))
      vpv=coef(1,2,iq)+y*(coef(2,2,iq)+y*(coef(3,2,iq)+y*coef(4,2,iq)))
      vsv=coef(1,3,iq)+y*(coef(2,3,iq)+y*(coef(3,3,iq)+y*coef(4,3,iq)))

cl--test * *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c	  print*,"la: voxint: radcontrib: calcqvec: qtau: getmod: SURVIVED!"


      if(.not.isotrp) then
         vph=coef(1,6,iq)+y*(coef(2,6,iq)+y*(coef(3,6,iq)+y*coef(4,6,iq)))
         vsh=coef(1,7,iq)+y*(coef(2,7,iq)+y*(coef(3,7,iq)+y*coef(4,7,iq)))
         eta=coef(1,8,iq)+y*(coef(2,8,iq)+y*(coef(3,8,iq)+y*coef(4,8,iq)))

cl--test * *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c	  print*,"la: voxint: radcontrib: calcqvec: qtau: getmod: anisotropic case"

      else
        vph=vpv
        vsh=vsv
        eta=1.d0
      end if
      return
      end

