      parameter (nknts=1000)
      parameter (mxdsc=30)
      real lcon,ncon
      common/modl1/n,nic,noc,moho,nsl,ifanis,r(nknts)
     1            ,rho(nknts),qrho(3,nknts),g(nknts),ell(nknts)
     2            ,eta(nknts)
      common/modl2/acon(nknts),qacon(3,nknts),ccon(nknts),qccon(3,nknts)
     1            ,lcon(nknts),qlcon(3,nknts),ncon(nknts),qncon(3,nknts)
     2            ,fcon(nknts),qfcon(3,nknts)
      common/modl3/qshear(nknts),qkappa(nknts)
      common/modl4/ndisc,ndsc(mxdsc)
      common/nond/rn,wn,vn,gn,rhobar


      common/mode/nord,jcom,lord,wcom,qbar,cgp,avert,ahor,phis,
     #            u(nknts),up(nknts),v(nknts),vp(nknts),
     #            p(nknts),pp(nknts)
      save
