C$PROG WIER
C     *** WIER ******** VERSION 1, MODIFICATION LEVEL 0 *** DKO00001 ***
C     *                                                                *
C     *   WRITE ERROR MESSAGES FOR SL-MATH SUBROUTINES                 *
C     *                                                                *
C     *   5736-XM7 COPYRIGHT IBM CORP. 1971                            *
C     *   REFER TO INSTRUCTIONS ON COPYRIGHT NOTICE FORM NO. 120-2083  *
C     *   FE SERVICE NO. 200281                                        *
C     *                                                                *
C     ******************************************************************
C
      SUBROUTINE WIER(IER,NO)
      WRITE (6,1) NO,IER
    1 FORMAT (//' ***** DKO',I5,' RAISED ERROR INDICATOR TO ',I4,
     1     ' *****'///)
      RETURN
      END
