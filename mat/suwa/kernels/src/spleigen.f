      subroutine spleigen(imodel,lora,ommin,ommax,omsp0,domsp,
     #                    nspl,splef,maxspl,maxefs)
      include 'gemodes.h'
      parameter (twopi=6.2831853072)
      common /plevel/iprtlv
      dimension splef(maxspl,maxefs,2)
      parameter (maxmodes=200)
      parameter (maxsplines=50)
      parameter (maxfns=12)
      dimension ff(maxmodes,maxfns)
      dimension om(maxmodes)
      dimension amat(maxsplines*maxmodes)
      dimension work(maxsplines*2)
      real*8 dbata(maxsplines,maxsplines)
      logical jtnorm
c
      jtnorm=.true.
c
c---- loop on all modes of desired type and save all wanted ef variables
c
      if(lora.eq.1) then
        ii=2
      else if(lora.eq.2) then
        ii=3
      endif
      nn=0
      nmods=0
      first=0.
      xlast=0.
      do ll=0,1000
        call findmode(imodel,nn,ii,ll,imod)
        if(imod.gt.0) then
          if(omegamod(imod,imodel).gt.omsp0.and.
     #       omegamod(imod,imodel).lt.ommax+domsp*0.5) then
            xk=float(ll)+0.5
            ufactor=omegamod(imod,imodel)/
     #              sqrt(pvelmod(imod,imodel)*gvelmod(imod,imodel))
            vfactor=xk*omegamod(imod,imodel)/
     #              sqrt(pvelmod(imod,imodel)*gvelmod(imod,imodel))
            nmods=nmods+1
            om(nmods)=omegamod(imod,imodel)
            ff(nmods,1)=pvelmod(imod,imodel)
            ff(nmods,2)=gvelmod(imod,imodel)
            ff(nmods,3)=xk
            ff(nmods,4)=smallqmod(imod,imodel)
            ff(nmods,5)=umod(imod,imodel)/radiusmod(imodel)
            ff(nmods,6)=upmod(imod,imodel)
            ff(nmods,7)=vmod(imod,imodel)/radiusmod(imodel)
            ff(nmods,8)=vpmod(imod,imodel)
            ff(nmods,9)=vaccmod(imod,imodel)
            ff(nmods,10)=haccmod(imod,imodel)
            ff(nmods,11)=uppmod(imod,imodel)*radiusmod(imodel)
            ff(nmods,12)=vppmod(imod,imodel)*radiusmod(imodel)
            if(jtnorm) then
              ff(nmods,5)=ff(nmods,5)*ufactor
              ff(nmods,6)=ff(nmods,6)*ufactor
              ff(nmods,7)=ff(nmods,7)*vfactor
              ff(nmods,8)=ff(nmods,8)*vfactor
              ff(nmods,9)=ff(nmods,9)*ufactor
              ff(nmods,10)=ff(nmods,10)*vfactor
              ff(nmods,11)=ff(nmods,11)*ufactor
              ff(nmods,12)=ff(nmods,12)*vfactor
            endif
c
c---- adopt the convention that the acceleration is always negative
c---- i.e. Utt and Wtt are negative
c
            if(lora.eq.1.and.ff(nmods,10).gt.0.) then
              ff(nmods,7)=-ff(nmods,7)
              ff(nmods,8)=-ff(nmods,8)
              ff(nmods,10)=-ff(nmods,10)
              ff(nmods,12)=-ff(nmods,12)
            endif
            if(lora.eq.2.and.ff(nmods,9).gt.0.) then
              ff(nmods,5)=-ff(nmods,5)
              ff(nmods,6)=-ff(nmods,6)
              ff(nmods,7)=-ff(nmods,7)
              ff(nmods,8)=-ff(nmods,8)
              ff(nmods,9)=-ff(nmods,9)
              ff(nmods,10)=-ff(nmods,10)
              ff(nmods,11)=-ff(nmods,11)
              ff(nmods,12)=-ff(nmods,12)
            endif
c
            if(nmods.eq.1) first=omegamod(imod,imodel)
            xlast=omegamod(imod,imodel)
          endif
        endif
      enddo
      nspl=int((ommax-omsp0)/domsp)
c
      if(iprtlv.gt.1) then
        if(lora.eq.1) then
          write(6,"('storing',i3,' knots for',i3,' T modes')") nspl,nmods
        else if(lora.eq.2) then
          write(6,"('storing',i3,' knots for',i3,' S modes')") nspl,nmods
        endif
        write(6,"('range:',2f10.3)") twopi/first,twopi/xlast
      endif
c
      if(nmods.lt.nspl) stop
c
c---- calculate the ATA^-1 matrix for a given number of spline points
c
      call findcbsplmat(om,nmods,omsp0,domsp,nspl,amat,dbata,work)
c
c---- calculate the b-spline coefficients for each ef variable
c
      do ief=1,12
        call findcbsplcoef(ff(1,ief),nmods,amat,dbata,nspl,splef(1,ief,lora),work)
      enddo
      return
      end
c
      subroutine findcbsplmat(xpt,npt,x0sp,dxsp,nkn,aa,ata,work)
c-----this subroutine finds a cubic spline that best fits a set of data
c-----points in a least-squares sense. It may blow up if there are not
c-----enough points to fit!
c
      dimension work(nkn)
      dimension xpt(npt)
      dimension aa(npt,nkn)
      real*8 ata(nkn,nkn)
      real*8 deter
c
      parameter (maxspl=50)
c
      do ikn=1,nkn
        do j=1,nkn
          work(j)=0.
          if(j.eq.ikn) work(j)=1.
        enddo
c
        do ipt=1,npt
          value=ecbspl(xpt(ipt),x0sp,dxsp,nkn,work)
          aa(ipt,ikn)=value
        enddo
      enddo
c
c---- make the inner product matrix
c
      do ikn=1,nkn
        do jkn=1,nkn
          ata(ikn,jkn)=0.d0
          do ipt=1,npt
            ata(ikn,jkn)=ata(ikn,jkn)+dble(aa(ipt,ikn)*aa(ipt,jkn))
          enddo
        enddo
      enddo
c
c---- invert it
c
      call matinv(ata,nkn,nkn,deter)
      if(iprtlv.gt.1 ) then
        write(6,"('determinant:',g15.5)") deter
      endif
c
      return
      end
c
c
      subroutine findcbsplcoef(ypt,npt,aa,ata,nkn,ysp,atb)
      dimension ypt(npt)
      dimension aa(npt,nkn)
      real*8 ata(nkn,nkn)
      dimension ysp(nkn)
      real*8 atb(nkn)
c
c---- make the coefficients
c
      do ikn=1,nkn
        atb(ikn)=0.d0
        do ipt=1,npt
          atb(ikn)=atb(ikn)+dble(aa(ipt,ikn)*ypt(ipt))
        enddo
      enddo
      do ikn=1,nkn
        ysp(ikn)=0.
        do jkn=1,nkn
          ysp(ikn)=ysp(ikn)+sngl(ata(ikn,jkn)*atb(jkn))
        enddo
      enddo
c
      return
      end
