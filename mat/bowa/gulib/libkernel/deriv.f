      subroutine deriv(pp,qvec,xtu)
c   deriv accepts as input the ray parameter pp in seconds/degree
c   and returns array qvec:
c   qvec(1)=distance; qvec(2)=dxdp; qvec(3)=time; qvec(4)= tstar
c   ...rest are derivatives of travel time w.r.t. the model 
c   coefficients for each of vph,vpv,vsh,vsv,and eta - there are
C   5 physical params each described by a cubic (for prem) so there 
C   are 20(=nplay) model parameters per layer
c   xtu is the turning point radius 
c
c   further input must be passed by the calling routine through
c   common blocks:
c         common/layr/nl(20),xb(20),xt(20),ifanis,nplay
c         common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
c   where
c         nl(i) = number of levels in layer i, xb(i) and xt(i) are
c                 bottom and top in kms. of layer i
c         rdep=source radius in kilometers
c         itype=1,2,3 for sh,p,sv waves
c         imod=0 for anisotropic model,1 for isotropic
c         itim: =1 to calculate only the deltas
c               =2 to calculate deltas and dXdp
c               =3 to calculate deltas, travel times, and deriva-
c                  tives of delta w.r.t. p and tstar
c          ider=1 to calculate model derivs
c         idleg=degree of the legendre polynomials to be used in
c                 the integrations
c         rnorm=length scale with respect to which the polynomial
c                 coefficients are normalized
c         numlyr=number of layers in the model - up to 20
c
      implicit double precision (a-h,o-z)
      logical isplit
c*** qints must be dimensioned large enough to be nplay+4
      dimension qvec(*),qints(24)
      common/path$/delt(2,50000,2),nfin(2),kdep !la 1600->50000
      common/ray$/pmin,pmax,nps(20,2),npsrc(2),nn,iup,jref,npleg,nsleg
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
      common/layr/nl(20),xb(20),xt(20),ifanis,nplay
      common/slow/picblp,picbls,picbup,pcmblp,pcmbup,pcmbus,
     +   p660lp,p660ls,lic,loc,llm,lsrc,klay(20)
      common/xnorm/xnorm
c -----the following is for the anisotropic inversion, Yu--------------
      common/kernel/sum(362,20,10),ifsplit,ifdiff
c----------------------------------------------------------------------
      data pi/3.14159265358979d0/
      rad=180.d0/pi
      p=pp*rad
      npar=numlyr*nplay+4
      do j=1,npar
        qvec(j)=0.d0
      enddo
c*** this is loop over ray type
      do 1000 ijk=1,2
      if(ijk.eq.1.and.npleg.eq.0) goto 1000
      if(ijk.eq.2.and.nsleg.eq.0) goto 1000
       
c*** itype=1 for SH, =2 for P, = 3 for SV
      itype=ijk+1
      if(ijk.eq.2.and.ish.ne.0) itype=1
c*** jpath counts increments in delta used in making path
      jpath=0
c ===================================================================
c ----loop over all the layers crossed by the ray
c isplit(logical): true if the source is in the current (n-th) shell
c rdep  : distance of the epicenter from the center of the earth
c iflag : =0(if ray traversed entire layer n); =1(if ray is reflected
c         at the top of the layer n); =2(if ray turned within layer n);
c         its value returns from the call to subroutine integr
c         n ---- PREM layer number

      do 100 n=1,nn
        nlast=n
        if(n.ne.1) jpath=klay(n-1)
        if(nps(n,ijk).eq.0) go to 100
        isplit=.false.
        if(n.eq.lsrc) then
          xbsav=xb(n)
          xtsav=xt(n)
          if(rdep.ne.xt(n))then
c    if the focus is in this shell, the bottom of the shell is
c    set equal to rdep and the upper half-shell is done, the
c    lower half-shell is done by means of the goto 105 at the
C    end of the do 100 block
             xb(n)=rdep
             isplit=.true.
          endif
        end if
  105   continue
        call integr(n,qints,iflag,jpath,ijk)
c -------ray reflected at the top of layer n
        if(iflag.eq.1) then
           nlast=n-1
           xtu=xt(n)
           go to 1000
        endif
        fmult=nps(n,ijk)
        if(n.eq.lsrc.and.isplit) fmult=npsrc(ijk)
        do i=1,4
          qvec(i)=qvec(i)+qints(i)*fmult
        enddo
        if(ider.ne.0) then
          ind=(n-1)*nplay+4
          do i=1,nplay
            qvec(i+ind)=qvec(i+ind)+qints(4+i)*fmult
          enddo
        end if
c*** if this is split source layer -- go back for bottom bit
        if(n.eq.lsrc) then
          xb(n)=xbsav
          if(isplit)then
            xt(n)=rdep
            isplit=.false.
            goto 105 
          else
            xt(n)=xtsav
          endif
       end if
c ---- if flag=2 ray turned in this layer; xtu (turning point
c ----radius) is set equal to xnorm (through common block /xnorm/)
        if(iflag.eq.2)then
          xtu=xnorm
          goto 1000
        endif
  100   continue
c ======================================================================
 1000 nfin(ijk)=jpath
c*** end of loop over wave type
      qvec(1)=qvec(1)*rad
      qvec(2)=qvec(2)*rad*rad

      if(jref.eq.1) xtu=xt(loc)
      if(jref.eq.2) xtu=xt(lic)
      return
      end
