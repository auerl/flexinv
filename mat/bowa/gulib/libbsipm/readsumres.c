#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>
#define resrcl 72      /* Xianfeng's residual file record length */

int readsumres_(fd, elat,elon,edep,slat,slon,res,ttime,delta)
float *elat,*elon,*edep,*slat,*slon;
float *res,*ttime,*delta;
int fd;
{  /* read summary data written in C */
int byterd, ig, ierr;
float stelv, elipc, crustcor, evcor, coric;
char eqname[9], stname[5];

typedef struct TEMP{		/* structure defined for reading only */
	char	eqname[8];	/* earthquake name */
	float elat;		/* summary earthquake latitude (in degree) */
	float elon;		/* longitude */
	float edep;		/* earthquake depth */
	char	stname[4];	/* station name */
	float slat;		/* summary station latitude (in degree) */
	float slon;		/* summary station longitude (in degree) */
	float stelv;		/* station elevation (in km) */
	float res;		/* observed residual */
	int	ig;		/* grade */
	float ttime;		/* total travel time */
	float elipc;		/* ellipcity correction */
	float crustcor;		/* crustal correction */
	float evcor;		/* station elevation correction */
	float delta;		/* epicentral distance */
	float coric;		/* innor core correction */
	float spare;
} TEMP;

TEMP temp;

  if((byterd=read(fd,&temp,resrcl))!=resrcl)return 0;
  strncpy(eqname,temp.eqname,8);
  eqname[8]='\0';
  *elat=temp.elat;
  *elon=temp.elon;
  *edep=temp.edep;
  strncpy(stname,temp.stname,4);
  stname[4]='\0';
  *slat=temp.slat;
  *slon=temp.slon;
  stelv=temp.stelv;
  *res=temp.res;
  *ttime=temp.ttime;
  elipc=temp.elipc;
  crustcor=temp.crustcor;
  evcor=temp.evcor;
  *delta=temp.delta;
  coric=temp.coric;
  ig=temp.ig;
  ierr = 1;

  printf("%s %f %f %f %f %f\n", eqname,*elat,*elon,*edep,*delta,*res);          return ierr;

}
