#include <sys/types.h>
/* #include <rmt.h> */  /* omit this line if remote tape access library not installed*/
#include <sys/stat.h>
#include <sys/times.h>
#include <sys/ioctl.h>
#include <sys/mtio.h>
#include <sys/file.h>
#include <sys/dir.h>
/* #include <sys/ttold.h> */
#include <stdio.h>
#include <unistd.h>
/* #include <sys/fcntlcom.h> */

/* struct sgttyb   pack; */
struct stat     buf;

cmkdir_(pname,mode,ires,ierrno)
        long *ires,*ierrno,*mode;
        char *pname;
{

        extern int errno, mkdir();
        errno=0;
        *ires=(long) mkdir(pname,(int)(*mode));
        *ierrno=(long)errno;
}
 
cstexc_(ichan, ires, ierrno)
	long           *ichan, *ires, *ierrno;
{
	extern int      errno, ioctl();
	char           *dummy;
	errno = 0;
	*ires = (long) ioctl((int) (*ichan), TIOCEXCL, dummy);
	*ierrno = (long) errno;
}

copen_(pname, ichan, iopt, ierrno, inew, mode)
	char           *pname;
	long           *ichan, *iopt, *ierrno, *inew, *mode;
{
	extern int      errno;
	extern int      open();
	int             kopt, knew;
	kopt = (int) (*iopt);
	knew = (int) (*inew);
	switch (kopt) {
	case 0:
		kopt = O_RDONLY;
		break;
	case 1:
		kopt = O_WRONLY;
		break;
	case 2:
		kopt = O_RDWR;
		break;
	case 8:
		kopt = O_APPEND;
		break;
	default:
		fprintf(stderr, "copen: unknown option %d", *iopt);
		exit(2);
		break;
	}
	switch (knew) {
	case 0:
		break;
	case 1:
		kopt = kopt | O_EXCL | O_CREAT;
		break;
	case 2:
		kopt = kopt | O_TRUNC | O_CREAT;
		break;
	default:
		fprintf(stderr, "copen: unknown status %d", *inew);
		exit(2);
		break;
	}
	errno = 0;
	*ichan = (long) open(pname, kopt, (int)(*mode) );
	*ierrno = (long) errno;
}

cclose_(ichan, ires, ierrno)
	long           *ichan, *ires, *ierrno;
{
	extern int      errno;
	extern int      close();
	errno = 0;
	*ires = (long) close((int) (*ichan));
	*ierrno = (long) errno;
}

cread_(ichan, pbuf, pnbyt, ires, ierrno)
	char           *pbuf;
	long           *ichan, *pnbyt, *ires, *ierrno;
{
	extern int      errno, read();
	errno = 0;
	*ires = (long) read((int) (*ichan), pbuf, (int) (*pnbyt));
	*ierrno = (long) errno;
}

clseek_(ichan, offst, iopt, ires, ierrno)
	long           *ichan, *offst, *iopt, *ires, *ierrno;
{
/*	extern int      errno, lseek(); */
	extern int      errno, llseek();
	int             kopt;
	kopt = (int) (*iopt);
	switch (kopt) {
	case 0:
		kopt = L_SET;
		break;
	case 1:
		kopt = L_INCR;
		break;
	case 2:
		kopt = L_XTND;
		break;
	default:
		fprintf(stderr, "clseek: unknown optiion %d", *iopt);
		exit(2);
		break;
	}
	errno = 0;
/*	*ires = (long) lseek((int) (*ichan), (*offst), kopt); */
	*ires = (long) llseek((int) (*ichan), (*offst), kopt);
	*ierrno = (long) errno;
}

cmtio_(ichan, iop, icnt, ires, ierrno)
	long           *ichan, *iop, *icnt, *ires, *ierrno;
{
	extern int      errno, ioctl();
	int             kopt;
	struct mtop     magop;
	kopt = (int) (*iop);
	switch (kopt) {
	case 0:
		kopt = MTWEOF;
		break;
	case 1:
		kopt = MTFSF;
		break;
	case 2:
		kopt = MTBSF;
		break;
	case 3:
		kopt = MTFSR;
		break;
	case 4:
		kopt = MTBSR;
		break;
	case 5:
		kopt = MTREW;
		break;
	case 6:
		kopt = MTOFFL;
		break;
	case 7:
		kopt = MTNOP;
		break;
	default:
		fprintf(stderr, "cmtio: unknown optiion %d", *iop);
		exit(2);
		break;
	}
	magop.mt_op = kopt;
	magop.mt_count = (long) (*icnt);
	errno = 0;
	*ires = (long) ioctl((int) (*ichan), MTIOCTOP, &magop);
	*ierrno = (long) errno;
}

csleep_(psec)
	long           *psec;
{
	int             idum;
	idum = sleep((unsigned) (*psec));
}

cusleep_(psec)
	long           *psec;
{
	int             idum;
	idum = usleep((unsigned) (*psec));
}

