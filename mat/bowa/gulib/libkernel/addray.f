      subroutine addray (nps,ic1,ic2,i)
c Keeps sum of number of P and S legs in each layer of the earth model.
      dimension nps(*)
      do 1 j=ic1,ic2
    1 nps(j)=nps(j)+i
      return
      end 
