      double precision function f(pp)
c   f is delta as a function of p if itim=1 in calling routine
c   f is d(delta)/d(p) as a function of p if itim=2
      implicit double precision (a-h,o-z)
      dimension qvec(500)
      common/in/p,rdep,itype,imod,itim,idleg,rnorm,numlyr,iso,ider,ish
c -----the following is for the anisotropic inversion, Yu--------------
      common/kernel/sum(362,20,10),ifsplit,ifdiff
c----------------------------------------------------------------------

      call deriv(pp,qvec,xtu)
      f=qvec(itim)
      return
      end

