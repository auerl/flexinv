c Inversion parameters:
c maxfreq --->  maximum tesselation frequency
c maxnode --->  maximum number of horizontal nodes
c maxrknot --->  maximum number of radial parameters
c maxparm --->  maximum number of unknown parameters in model
c Y.G., 2002.
c
c
      parameter(maxfreq=8)
      parameter(maxnode=maxfreq*maxfreq*10+2)
      parameter(maxrknot=1)
      parameter(maxparm=maxnode*maxrknot)

      data pi /3.141592653579/
      data rad/57.2957795130824/
      data reprad /57.295779579/    
