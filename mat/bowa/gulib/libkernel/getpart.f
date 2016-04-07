	subroutine getpart(qray,rpath,part)
c This program calculates the kernel G for integration of G*f(k)
c evaluation of G uses cubic polynomials saved in qray
c Input:	qray --- 3-D array of model derivatives
c		rpath -- radius
c Output:	
c		partial derivative of Qtau with Vph,Vpv,Vsh,Vsv,eta
c		saved in part(1:5) in the order above.
c
        implicit real*8(a-h,o-z)
	dimension qray(4,5,20), part(5)
        logical isotrp
        dimension test(4)
        data iqs/-1/,test/1.d0,0.d0,0.d0,0.d0/
	common/layr/nl(20),xb(20),xt(20),ifanis,nplay
	common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
	common/anisopath/theta_path(3200),rad_path(3200),derv_path(3200,10)
     +       ,ntheta, nraytp(20), nrayst(20), ileg
	do ip=1, 5
		part(ip)=0.d0
	enddo
c	write(*,*) '# layers in the model = ', numlyr
	ierror = 1
	i1=1
	iq=0
	do j=1, numlyr
		if(rpath.ge.xb(j).and.rpath.lt.xt(j)) then
			iq = j
			y=rpath/rnorm
			ierror = 0
			goto 10
		endif
	enddo
10	if(ierror.eq.1) then
		write(*,*) 'error in the raypath in getgreen.f'
		stop
	endif
c	write(*,*) 'y=', y, ' iq = ', iq, ' xb(iq)=', xb(iq), ' xt(iq)=', xt(iq), ' rpath=', rpath
c	write(*,*) iq, qray(1,1,iq), qray(1,2,iq), qray(1,3,iq), qray(1,4,iq), qray(1,5,iq)
	do ip=1, 5
           part(ip)=qray(1,ip,iq)+y*(qray(2,ip,iq)+y*(qray(3,ip,iq)+y*qray(4,ip,iq)))
	enddo
c	write(*,*) part(1),part(2),part(3),part(4),part(5)
        return
      end
