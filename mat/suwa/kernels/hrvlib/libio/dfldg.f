      DOUBLE PRECISION FUNCTION DFLDG(X)
      DOUBLE PRECISION X,Y,SMALL,D1,D2,D3,D4
      DATA SMALL/1.52587890625D-05/
      INTEGER*2 I(4),IMANT,IEXP,MINUS8,JJ(2),FIFTEEN
      EQUIVALENCE (Y,I(1)),(JJ(1),J)
      DATA IMANT/Z'00FF'/,IEXP/Z'7F00'/,MINUS8/-8/,JJ/0,0/
     1 ,FIFTEEN/15/
      Y=X                                        



      D1=AND(I(1),IMANT)
      JJ(2)=I(2)

      D2=J                              


      JJ(2)=I(3)
      D3=J                            


      JJ(2)=I(4)
      D4=J                           


      IPOW=ISHFT(AND(IEXP,I(1)),MINUS8)-66
      DFLDG=(16.D0**IPOW)*(D1+SMALL*(D2+SMALL*(D3+SMALL*D4)))
      IF(BTEST(I(1),FIFTEEN)) DFLDG=-DFLDG
      RETURN
      END
