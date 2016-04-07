      subroutine readata(io,natd,atd,ata)
c
c... read ATA matrix
c... Input:
c		natd  --- number of ATD elements
c... Output:	
c		atd   --- ATD matrix
c		ata   --- ATA matrix

      dimension atd(1), ata(1)

      write(*,"('Reading ATA matrix...')")
      nata=natd*(natd+1)/2
      write(*,"('number of ATD: ',i14)") natd
      write(*,"('number of ATA: ',i14)") nata

      read(io) natd
      read(io) (atd(i),i=1,natd)
      read(io) (ata(i),i=1,nata)
      close(io)
      return
      end
            
