      subroutine writeata(io,natd,atd,ata)
c
c... output ATA matrix
c
      dimension atd(1), ata(1)

      write(*,"('Output ATA matrix...')")
      nata=natd*(natd+1)/2
      write(*,"('number of ATD: ',i14)") natd
      write(*,"('number of ATA: ',i14)") nata

      write(io) natd
      write(io) (atd(i),i=1,natd)
      write(io) (ata(i),i=1,nata)
      close(io)
      return
      end
            
