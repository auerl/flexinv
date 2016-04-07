#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#define  PERM  0644   /* read write for owner */

void writeipm_c_(leny,nodelat, nodelon,nrad,radknot,nsorce,nstruc,maskst,maskso,numata,numatd,ata,atd)
int *leny, *nrad, *nsorce, *nstruc,*maskst,*maskso,*numata,*numatd;
float ata[], atd[], radknot[];
float nodelon[], nodelat[];
{
  int i,k,fd,nwritten;
  int format, ntopo;
  float theta_phi[2];
  char  outfile1[80];
  float RAD=0.01745329251994;
  float diagavg;

   printf("output inner product matrix name:\n");
   scanf("%s", outfile1);
   fd = creat(outfile1, PERM);
   if(fd == -1) {
      printf(" error in opening file %s\n", outfile1);
   }
   format = 6 ; /* b-spline format with topography  */
   nwritten=write(fd,&format,sizeof(int));
/* write out model knots (horizontal) */
   nwritten=write(fd,leny,sizeof(int));
   printf("output leny = %d\n", *leny);
   for(i=0;i<*leny;i++){
/*    
      printf("nodelat=%f, nodelon=%f\n", nodelat[i], nodelon[i]);
*/
      theta_phi[0]=(90.0-nodelat[i])*RAD;
      theta_phi[1]=nodelon[i];
      if(theta_phi[1]<0.0) theta_phi[1]+=360.0;
      theta_phi[1]=theta_phi[1]*RAD;
      nwritten=write(fd,theta_phi,2*sizeof(float));
   }
/* write out model knots (radial) */
   printf("number of radial knots in writeipm_cpp= %d\n", *nrad);
   printf("number of atd element=%d\n, number of ata=%d\n", *numatd, *numata);
   nwritten=write(fd, nrad, sizeof(int));
   nwritten=write(fd, radknot, (*nrad)*sizeof(float));
/* write out matrices */
   nwritten=write(fd,numatd, sizeof(int));
   k=0;
   diagavg=0.0;
   for(i=0;i<*numatd; i++){
	  if(ata[k]<0.0) {
	     printf("\n Error, ATA diagnal less than 0.0!!!\n\n");
	     exit(-1);
	  }
	  diagavg += ata[k];
	  k+=i+2;
   }
   printf("average of diagonal elements= %f\n", diagavg/(*numatd));   
   ntopo=0;
/*
 *  below finds the number of surface b-spline layers in the inversion
 */
   for(i=0;i<2; i++){
      if(maskst[i] != 0) ntopo++;
   }
   printf("number of topography parameters = %d\n", ntopo);
   nwritten=write(fd,&ntopo, sizeof(int));
   nwritten=write(fd,ata,(*numata)*sizeof(float));
   nwritten=write(fd,atd,(*numatd)*sizeof(float));
   printf( "output ata file %s", outfile1);
   close(fd);
   return;
}
