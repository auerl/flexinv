      parameter (mxmodels=2)
      parameter (mxlevels=1000)
      real vpvlayer(mxlevels,mxmodels)
      real vphlayer(mxlevels,mxmodels)
      real vsvlayer(mxlevels,mxmodels)
      real vshlayer(mxlevels,mxmodels)
      real etalayer(mxlevels,mxmodels)
      real radlayer(mxlevels,mxmodels)
      real rholayer(mxlevels,mxmodels)
      real qmulayer(mxlevels,mxmodels)
      real qkalayer(mxlevels,mxmodels)
      character*80 modelname(2)
      integer ifanismodel(mxmodels)
      integer ifdeckmodel(mxmodels)
      real trefmodel(mxmodels)
      integer nmodels
      integer nlevmodel(mxmodels)
      integer nlevoutmodel(mxmodels)
      integer ninnercore(mxmodels)
      integer noutercore(mxmodels)
      integer nmohorovic(mxmodels)
      integer ioutlayer(mxlevels,mxmodels)
c
      integer lminmodel(mxmodels)
      integer lmaxmodel(mxmodels)
      integer nbranmodel(mxmodels)
      integer fminmodel(mxmodels)
      integer fmaxmodel(mxmodels)
      real epsmodel(mxmodels)
      real wgravmodel(mxmodels)
c
      common /modelpara/ nmodels,
     #        nlevmodel,nlevoutmodel,
     #        ninnercore,noutercore,nmohorovic,ifanismodel,
     #        ifdeckmodel,trefmodel,
     #        lminmodel,lmaxmodel,nbranmodel,
     #        fminmodel,fmaxmodel,epsmodel,wgravmodel,
     #        vpvlayer,vphlayer,vsvlayer,vshlayer,
     #        etalayer,radlayer,rholayer,
     #        qmulayer,qkalayer,
     #        ioutlayer,
     #        modelname
c
      parameter (mxmodes=12000)
      integer ifrstrad(mxmodels)
      integer ifrstsph(mxmodels)
      integer ifrsttor(mxmodels)
      integer nbytesmo(mxmodels)
      integer ibytfrst(mxmodels)
      integer nmodemod(mxmodels)
      integer imodefil(mxmodes,mxmodels)
      integer iovermod(mxmodes,mxmodels)
      integer itypemod(mxmodes,mxmodels)
      integer ilordmod(mxmodes,mxmodels)
      real    omegamod(mxmodes,mxmodels)
      real    smallqmod(mxmodes,mxmodels)
      real    vaccmod(mxmodes,mxmodels)
      real    haccmod(mxmodes,mxmodels)
      real    vdismod(mxmodes,mxmodels)
      real    hdismod(mxmodes,mxmodels)
      real    gvelmod(mxmodes,mxmodels)
      real    pvelmod(mxmodes,mxmodels)
      real    potmod(mxmodes,mxmodels)
      real    umod(mxmodes,mxmodels)
      real    upmod(mxmodes,mxmodels)
      real    uppmod(mxmodes,mxmodels)
      real    vmod(mxmodes,mxmodels)
      real    vpmod(mxmodes,mxmodels)
      real    vppmod(mxmodes,mxmodels)
c
      real    depthmod(mxmodels)
      real    radiusmod(mxmodels)
c      
      common /modepara/ depthmod,radiusmod,
     #        ifrstrad,ifrstsph,ifrsttor,nbytesmo,ibytfrst,
     #        nmodemod,imodefil,iovermod,itypemod,ilordmod,omegamod,
     #        smallqmod,vaccmod,haccmod,vdismod,hdismod,
     #        gvelmod,pvelmod,potmod,
     #        umod,upmod,uppmod,vmod,vpmod,vppmod
      
      save
