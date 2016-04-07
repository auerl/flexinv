#include <sys/types.h>
/*#include <rmt.h>*/  /* omit this line if remote tape access library not installed*/
#include <sys/stat.h>
#include <sys/ioctl.h>
#include <sys/mtio.h>
#include <sys/file.h>
#include <stdio.h>


#define TIOCSETP        0x7409
#define TIOCGETP        0x7408

struct sgttyb   pack;
struct stat     buf;

cgtflag_(ichan, iflag, ires, ierrno)
	long           *ichan, *iflag, *ires, *ierrno;
{
    extern int      errno;//, ioctl();
	errno = 0;
	*ires = (long) ioctl((int) (*ichan), TIOCGETP, &pack);
	*ierrno = (long) errno;
	*iflag = (long) pack.sg_flags;
}

cstflag_(ichan, iflag, ires, ierrno)
	long           *ichan, *iflag, *ires, *ierrno;
{
    extern int      errno;//, ioctl();
	errno = 0;
	pack.sg_flags = (int) (*iflag);
	*ires = (long) ioctl((int) (*ichan), TIOCSETP, &pack);
	*ierrno = (long) errno;
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

cnoblock_(ichan, iopt, ires, ierrno)
	long           *ichan, *iopt, *ires, *ierrno;
{
	extern int      errno, ioctl();
	int             tiopt;
	tiopt = (int) (*iopt);
	errno = 0;
	*ires = (long) ioctl((int) (*ichan), FIONBIO, &tiopt);
	*ierrno = (long) errno;
}

copen_(pname, ichan, iopt, ierrno)
	char           *pname;
	long           *ichan, *iopt, *ierrno;
{
	extern int      errno;
	extern int      open();
	int             kopt;
	kopt = (int) (*iopt);
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
	errno = 0;
	*ichan = (long) open(pname, kopt);
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

cwrite_(ichan, pbuf, pnbyt, ires, ierrno)
	char           *pbuf;
	long           *ichan, *pnbyt, *ires, *ierrno;
{
	extern int      errno, write();
	errno = 0;
	*ires = (long) write((int) (*ichan), pbuf, (int) (*pnbyt));
	*ierrno = (long) errno;
}

clseek_(ichan, offst, iopt, ires, ierrno)
	long           *ichan, *offst, *iopt, *ires, *ierrno;
{
	extern int      errno, lseek();
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
	*ires = (long) lseek((int) (*ichan), (*offst), kopt);
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


cperror_(s)
	char           *s;
{
	int             idum;
	idum = perror(s);
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

cfstat_(ichan, size, istat, ierrno)
	long           *ichan, *size, *istat, *ierrno;
{
	extern int      errno;
	extern int      fstat();
	int             dummy;
	errno = 0;
	/* printf("calling fstat in c\n"); */
	dummy = fstat((int) (*ichan), &buf);
	/* printf("fstat through in c\n"); */
	*ierrno = (long) errno;
	*size = (long) (buf.st_size);
	*istat = (long) (buf.st_mode);
}

cstat_(file, size, istat, ierrno)
	char           *file;
	long           *size, *istat, *ierrno;
{
	extern int      errno;
	extern int      stat();
	int             dummy;
	errno = 0;
	dummy = stat(file, &buf);
	*ierrno = (long) errno;
	*size = (long) (buf.st_size);
	*istat = (long) (buf.st_mode);
}

ctrun_(ichan, leng, ierrno)
	long           *ichan, *leng, *ierrno;
{
	extern int      errno, ftruncate();
	int             dummy;
	errno = 0;
	dummy = ftruncate((int) (*ichan), (unsigned long) (*leng));
	*ierrno = (long) errno;
}

cgetppid_(pid)
	int		*pid;
{
	*pid = getppid();
}

creadlink_(path, file, namlen, ierrno)
	char		*path, *file;
	int		*namlen, *ierrno;
{
	extern int	errno, readlink();
	int		dummy;
	errno = 0;
	dummy = readlink(path, file, namlen);
	if (dummy = -1) { *ierrno = errno; };
	if (dummy!= -1) { *ierrno = dummy; };
}
