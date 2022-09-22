      SUBROUTINE EQUTOGAL(RA,DEC,L,B,CODE)

C     Converts equatorial coordinates (right ascension, declination)
C     to (new) galactic coordinates (l,b).

C     CODE = 0: Input in radians, output in radians
C            1: Input in radians, output in degrees
C            2: Input in hours and degrees, output in radians
C            3: Input in hours and degrees, output in degrees

      REAL*8 PI, ONEDEG, PIOV2, TWOPI, SINOBLIQ, COSOBLIQ, RA0,
     :       L0
      PARAMETER (PI=3.141592653589793238D0)
      PARAMETER (ONEDEG=PI/180., PIOV2=PI/2, TWOPI=2*PI)
      PARAMETER (COSOBLIQ=0.4601997847838517D0,
     :           SINOBLIQ=0.8878153851364013D0)
  
C                62.6 degrees

      PARAMETER (RA0=ONEDEG*282.25, L0=ONEDEG*33.0)
      INTEGER CODE
      REAL RA, DEC, L, B
      REAL*8 LONG, LAT, SINLONG, COSLONG, SINLAT, COSLAT, SINB,
     :       COSB, LONG2, SINLONG2, COSLONG2

      IF (CODE .EQ. 0 .OR. CODE .EQ. 1) THEN
            LONG = RA - RA0
            LAT = DEC
         ELSE IF (CODE .EQ. 2 .OR. CODE .EQ. 3) THEN
            LONG = ONEDEG*15*RA - RA0
            LAT = ONEDEG*DEC
         ELSE
            STOP 'STOP: Bad code passed to EQUTOGAL'
      END IF

      SINLAT = SIN(LAT)
      COSLAT = COS(LAT)
      SINLONG = SIN(LONG)
      COSLONG = COS(LONG)

      SINB = MIN(SINLAT*COSOBLIQ-COSLAT*SINLONG*SINOBLIQ,1.0D0)
      B = ASIN( MAX(SINB,-1.0D0) )
      COSB = MAX( SQRT(1 - SINB**2) , 1.D-7 )
      COSLONG2 = MIN(COSLAT*COSLONG/COSB,1.0D0)
      LONG2 = ACOS( MAX(COSLONG2,-1.0D0) )
      SINLONG2 = (COSLAT*SINLONG*COSOBLIQ + SINLAT*SINOBLIQ)
     :               / COSB
      IF (SINLONG2 .GE. 0) THEN
            L = LONG2 + L0
         ELSE
            L = (TWOPI-LONG2) + L0
      END IF
      IF (L .GE. TWOPI) L = L - TWOPI

      IF (CODE .EQ. 1 .OR. CODE .EQ. 3) THEN
         L = L/ONEDEG
         B = B/ONEDEG
      END IF

      END

C     ******************************************************************

      SUBROUTINE GALTOEQU(L,B,RA,DEC,CODE)

C     Converts (new) galactic coordinates (l,b) to equatorial
C     coordinates (right ascension, declination).

C     CODE = 0: Input and output are in radians
C            1: Input in radians, output in hours and degrees
C            2: Input in degrees, output in radians
C            3: Input in degrees, output in hours and degrees

      REAL*8 PI, ONEDEG, PIOV2, TWOPI, SINOBLIQ, COSOBLIQ, RA0,
     :       L0
      PARAMETER (PI=3.141592653589793238D0)
      PARAMETER (ONEDEG=PI/180., PIOV2=PI/2, TWOPI=2*PI)
      PARAMETER (COSOBLIQ=0.4601997847838517D0,
     :           SINOBLIQ=0.8878153851364013D0)  
c 62.6 degrees
      PARAMETER (RA0=ONEDEG*282.25, L0=ONEDEG*33.0)
      INTEGER CODE
      REAL RA, DEC, L, B
      REAL*8 LONG, LAT, SINLONG, COSLONG, SINLAT, COSLAT, SINDEC,
     :       COSDEC, LONG2, SINLONG2, COSLONG2

      IF (CODE .EQ. 0 .OR. CODE .EQ. 1) THEN
            LONG = L - L0
            LAT = B
         ELSE IF (CODE .EQ. 2 .OR. CODE .EQ. 3) THEN
            LONG = ONEDEG*L - L0
            LAT = ONEDEG*B
         ELSE
            STOP 'STOP: Bad code passed to GALTOEQU'
      END IF

      SINLAT = SIN(LAT)
      COSLAT = COS(LAT)
      SINLONG = SIN(LONG)
      COSLONG = COS(LONG)

      SINDEC = MIN(SINLAT*COSOBLIQ+COSLAT*SINLONG*SINOBLIQ,1.0D0)
      DEC = ASIN( MAX(SINDEC,-1.0D0) )
      COSDEC = MAX( SQRT(1 - SINDEC**2) , 1.D-7 )
      COSLONG2 = MIN(COSLAT*COSLONG/COSDEC,1.0D0)
      LONG2 = ACOS( MAX(COSLONG2,-1.0D0) )
      SINLONG2 = (COSLAT*SINLONG*COSOBLIQ - SINLAT*SINOBLIQ)
     :               / COSDEC
      IF (SINLONG2 .GE. 0) THEN
            RA = LONG2 + RA0
         ELSE
            RA = (TWOPI-LONG2) + RA0
      END IF
      IF (RA .GE. TWOPI) RA = RA - TWOPI

      IF (CODE .EQ. 1 .OR. CODE .EQ. 3) THEN
         RA = RA/(15*ONEDEG)
         DEC = DEC/ONEDEG
      END IF

      END

C     ******************************************************************

      SUBROUTINE EQUTOECL(RA,DEC,L,B,CODE)

C     Converts 1950 equatorial coordinates (right ascension,declination)
C     to 1950 ecliptic coordinates (lambda,beta).

C     CODE = 0: Input in radians, output in radians
C            1: Input in radians, output in degrees
C            2: Input in hours and degrees, output in radians
C            3: Input in hours and degrees, output in degrees

      REAL PI, ONEDEG, PIOV2, TWOPI, SINOBLIQ, COSOBLIQ
      PARAMETER (PI=3.141592653589793238)
      PARAMETER (ONEDEG=PI/180., PIOV2=PI/2, TWOPI=2*PI)
      PARAMETER (SINOBLIQ=0.39788119, COSOBLIQ=0.91743695) 
c! 23.45 deg
      INTEGER CODE
      REAL RA, DEC, L, B, LONG, LAT, SINLONG, COSLONG, SINLAT,
     :     COSLAT, SINB, COSB, LONG2, SINLONG2, COSLONG2

      IF (CODE .EQ. 0 .OR. CODE .EQ. 1) THEN
            LONG = RA
            LAT = DEC
         ELSE IF (CODE .EQ. 2 .OR. CODE .EQ. 3) THEN
            LONG = ONEDEG*15*RA
            LAT = ONEDEG*DEC
         ELSE
            STOP 'STOP: Bad code passed to EQUTOECL'
      END IF

      SINLAT = SIN(LAT)
      COSLAT = COS(LAT)
      SINLONG = SIN(LONG)
      COSLONG = COS(LONG)

      SINB = MIN(SINLAT*COSOBLIQ-COSLAT*SINLONG*SINOBLIQ,1.0)
      B = ASIN( MAX(SINB,-1.0) )
      COSB = MAX(COS(B),1.E-7)
      COSLONG2 = MIN(COSLAT*COSLONG/COSB,1.0)
      LONG2 = ACOS( MAX(COSLONG2,-1.0) )
      SINLONG2 = (COSLAT*SINLONG*COSOBLIQ + SINLAT*SINOBLIQ)
     :               / COSB
      IF (SINLONG2 .GE. 0) THEN
            L = LONG2
         ELSE
            L = TWOPI - LONG2
      END IF

      IF (CODE .EQ. 1 .OR. CODE .EQ. 3) THEN
         L = L/ONEDEG
         B = B/ONEDEG
      END IF

      END

C     ******************************************************************

      SUBROUTINE ECLTOEQU(L,B,RA,DEC,CODE)

C     Converts 1950 ecliptic coordinates (lambda,beta) to equatorial
C     coordinates (right ascension, declination).

C     CODE = 0: Input and output are in radians
C            1: Input in radians, output in hours and degrees
C            2: Input in degrees, output in radians
C            3: Input in degrees, output in hours and degrees

      REAL PI, ONEDEG, PIOV2, TWOPI, SINOBLIQ, COSOBLIQ
      PARAMETER (PI=3.141592653589793238)
      PARAMETER (ONEDEG=PI/180., PIOV2=PI/2, TWOPI=2*PI)
      PARAMETER (SINOBLIQ=0.39788119, COSOBLIQ=0.91743695) 
c! 23.45 deg
      INTEGER CODE
      REAL RA, DEC, L, B, LONG, LAT, SINLONG, COSLONG, SINLAT,
     :     COSLAT, SINDEC, COSDEC, LONG2, SINLONG2, COSLONG2

      IF (CODE .EQ. 0 .OR. CODE .EQ. 1) THEN
            LONG = L
            LAT = B
         ELSE IF (CODE .EQ. 2 .OR. CODE .EQ. 3) THEN
            LONG = ONEDEG*L
            LAT = ONEDEG*B
         ELSE
            STOP 'STOP: Bad code passed to GALTOEQU'
      END IF

      SINLAT = SIN(LAT)
      COSLAT = COS(LAT)
      SINLONG = SIN(LONG)
      COSLONG = COS(LONG)

      SINDEC = MIN(SINLAT*COSOBLIQ+COSLAT*SINLONG*SINOBLIQ,1.0)
      DEC = ASIN( MAX(SINDEC,-1.0) )
      COSDEC = MAX(COS(DEC),1.E-7)
      COSLONG2 = MIN(COSLAT*COSLONG/COSDEC,1.0)
      LONG2 = ACOS( MAX(COSLONG2,-1.0) )
      SINLONG2 = (COSLAT*SINLONG*COSOBLIQ - SINLAT*SINOBLIQ)
     :               / COSDEC
      IF (SINLONG2 .GE. 0) THEN
            RA = LONG2
         ELSE
            RA = TWOPI - LONG2
      END IF

      IF (CODE .EQ. 1 .OR. CODE .EQ. 3) THEN
         RA = RA/(15*ONEDEG)
         DEC = DEC/ONEDEG
      END IF

      END

C***********************************************************************

      SUBROUTINE EQUTOSGAL(RA,DEC,SGL,SGB,CODE)

C     Converts equatorial coordinates (right ascension, declination)
C     to (new) supergalactic coordinates (SGL,SGB), as given in RC2
C     (de Vaucouleurs et al. 1976).

C     The SG north pole is defined in old galactic coordinates
C     (l=15,b=+5) in RC1 and RC2, but we assume here that the
C     approximate value in new galactic coordinates (l=47.37,b=+6.32)
C     given in RC2 is exactly correct. Note that the origin of
C     coordinates (new galactic coordinates) given in RC2 is incorrect;
C     the printed value l=137.29 is the new galactic latitude of the
C     OLD supergalactic origin. The correct origin is at l=137.37, b=0. 

C     In converting to and from equatorial coordinates, I use rotation
C     angles which were determined on a TI-30 pocket calculator, and
C     which are probably accurate to a few thousandths of a degree.

C     CODE = 0: Input in radians, output in radians
C            1: Input in radians, output in degrees
C            2: Input in hours and degrees, output in radians
C            3: Input in hours and degrees, output in degrees

      REAL PI, ONEDEG, PIOV2, TWOPI, SINOBLIQ, COSOBLIQ, RA0,
     :     SGL0
      PARAMETER (PI=3.1415927)
      PARAMETER (ONEDEG=PI/180., PIOV2=PI/2, TWOPI=2*PI)
      PARAMETER (COSOBLIQ=0.269661, SINOBLIQ=0.962955)  
c! 74.356 deg
      PARAMETER (RA0=ONEDEG*13.189, SGL0=ONEDEG*-63.269)
      INTEGER CODE
      REAL RA, DEC, SGL, SGB
      REAL LONG, LAT, SINLONG, COSLONG, SINLAT, COSLAT, SINSGB,
     :     COSSGB, LONG2, SINLONG2, COSLONG2

      IF (CODE .EQ. 0 .OR. CODE .EQ. 1) THEN
            LONG = RA - RA0
            LAT = DEC
         ELSE IF (CODE .EQ. 2 .OR. CODE .EQ. 3) THEN
            LONG = ONEDEG*15*RA - RA0
            LAT = ONEDEG*DEC
         ELSE
            STOP 'STOP: Bad code passed to EQUTOGAL'
      END IF

      SINLAT = SIN(LAT)
      COSLAT = COS(LAT)
      SINLONG = SIN(LONG)
      COSLONG = COS(LONG)

      SINSGB = MIN(SINLAT*COSOBLIQ-COSLAT*SINLONG*SINOBLIQ,1.0)
      SGB = ASIN( MAX(SINSGB,-1.0) )
      COSSGB = MAX( SQRT(1 - SINSGB**2) , 1.E-7 )
      COSLONG2 = MIN(COSLAT*COSLONG/COSSGB,1.0)
      LONG2 = ACOS( MAX(COSLONG2,-1.0) )
      SINLONG2 = (COSLAT*SINLONG*COSOBLIQ + SINLAT*SINOBLIQ)
     :               / COSSGB
      IF (SINLONG2 .GE. 0) THEN
            SGL = LONG2 + SGL0
         ELSE
            SGL = (TWOPI-LONG2) + SGL0
      END IF
      IF (SGL .GE. TWOPI) SGL = SGL - TWOPI

      IF (CODE .EQ. 1 .OR. CODE .EQ. 3) THEN
         SGL = SGL/ONEDEG
         SGB = SGB/ONEDEG
      END IF

      END

C     ******************************************************************

      SUBROUTINE SGALTOEQU(SGL,SGB,RA,DEC,CODE)

C     Converts (new) supergalactic coordinates (SGL,SGB), as given in
C     RC2 (de Vaucouleurs et al. 1976), to equatorial coordinates
C     (right ascension, declination).

C     The SG north pole is defined in old galactic coordinates
C     (l=15,b=+5) in RC1 and RC2, but we assume here that the
C     approximate value in new galactic coordinates (l=47.37,b=+6.32)
C     given in RC2 is exactly correct. Note that the origin of
C     coordinates (new galactic coordinates) given in RC2 is incorrect;
C     the printed value l=137.29 is the new galactic latitude of the
C     OLD supergalactic origin. The correct origin is at l=137.37, b=0. 

C     In converting to and from equatorial coordinates, I use rotation
C     angles which were determined on a TI-30 pocket calculator, and
C     which are probably accurate to a few thousandths of a degree.

C     CODE = 0: Input and output are in radians
C            1: Input in radians, output in hours and degrees
C            2: Input in degrees, output in radians
C            3: Input in degrees, output in hours and degrees

      REAL PI, ONEDEG, PIOV2, TWOPI, SINOBLIQ, COSOBLIQ, RA0,
     :     SGL0
      PARAMETER (PI=3.1415927)
      PARAMETER (ONEDEG=PI/180., PIOV2=PI/2, TWOPI=2*PI)
      PARAMETER (COSOBLIQ=0.269661, SINOBLIQ=0.962955)  
c 74.256 degs
      PARAMETER (RA0=ONEDEG*13.189, SGL0=ONEDEG*-63.269)
      INTEGER CODE
      REAL RA, DEC, SGL, SGB
      REAL LONG, LAT, SINLONG, COSLONG, SINLAT, COSLAT, SINDEC,
     :     COSDEC, LONG2, SINLONG2, COSLONG2

      IF (CODE .EQ. 0 .OR. CODE .EQ. 1) THEN
            LONG = SGL - SGL0
            LAT = SGB
         ELSE IF (CODE .EQ. 2 .OR. CODE .EQ. 3) THEN
            LONG = ONEDEG*SGL - SGL0
            LAT = ONEDEG*SGB
         ELSE
            STOP 'STOP: Bad code passed to GALTOEQU'
      END IF

      SINLAT = SIN(LAT)
      COSLAT = COS(LAT)
      SINLONG = SIN(LONG)
      COSLONG = COS(LONG)

      SINDEC = MIN(SINLAT*COSOBLIQ+COSLAT*SINLONG*SINOBLIQ,1.0)
      DEC = ASIN( MAX(SINDEC,-1.0) )
      COSDEC = MAX( SQRT(1 - SINDEC**2) , 1.E-7 )
      COSLONG2 = MIN(COSLAT*COSLONG/COSDEC,1.0)
      LONG2 = ACOS( MAX(COSLONG2,-1.0) )
      SINLONG2 = (COSLAT*SINLONG*COSOBLIQ - SINLAT*SINOBLIQ)
     :               / COSDEC
      IF (SINLONG2 .GE. 0) THEN
            RA = LONG2 + RA0
         ELSE
            RA = (TWOPI-LONG2) + RA0
      END IF
      IF (RA .GE. TWOPI) RA = RA - TWOPI

      IF (CODE .EQ. 1 .OR. CODE .EQ. 3) THEN
         RA = RA/(15*ONEDEG)
         DEC = DEC/ONEDEG
      END IF

      END

C***********************************************************************

      SUBROUTINE GALTOSGAL(L,B,SGL,SGB,CODE)

C     Converts (new) galactic coordinates to (new) supergalactic
C     coordinates (SGL,SGB), as given in RC2 (de Vaucouleurs et al.
C     1976).

C     The SG north pole is defined in old galactic coordinates
C     (l=15,b=+5) in RC1 and RC2, but we assume here that the
C     approximate value in new galactic coordinates (l=47.37,b=+6.32)
C     given in RC2 is exactly correct. Note that the origin of
C     coordinates (new galactic coordinates) given in RC2 is incorrect;
C     the printed value l=137.29 is the new galactic latitude of the
C     OLD supergalactic origin. The correct origin is at l=137.37, b=0. 

C     CODE = 0: Input in radians, output in radians
C            1: Input in radians, output in degrees
C            2: Input in hours and degrees, output in radians
C            3: Input in hours and degrees, output in degrees

      REAL PI, ONEDEG, PIOV2, TWOPI, SINOBLIQ, COSOBLIQ, L0
      PARAMETER (PI=3.1415927)
      PARAMETER (ONEDEG=PI/180., PIOV2=PI/2, TWOPI=2*PI)
      PARAMETER (COSOBLIQ=0.1100813, SINOBLIQ=0.9939226)  
c! 83.68 deg
      PARAMETER (L0=ONEDEG*137.37)
      INTEGER CODE
      REAL L, B, SGL, SGB
      REAL LONG, LAT, SINLONG, COSLONG, SINLAT, COSLAT, SINSGB,
     :     COSSGB, LONG2, SINLONG2, COSLONG2

      IF (CODE .EQ. 0 .OR. CODE .EQ. 1) THEN
            LONG = L - L0
            LAT = B
         ELSE IF (CODE .EQ. 2 .OR. CODE .EQ. 3) THEN
            LONG = ONEDEG*L - L0
            LAT = ONEDEG*B
         ELSE
            STOP 'STOP: Bad code passed to EQUTOGAL'
      END IF

      SINLAT = SIN(LAT)
      COSLAT = COS(LAT)
      SINLONG = SIN(LONG)
      COSLONG = COS(LONG)

      SINSGB = MIN(SINLAT*COSOBLIQ-COSLAT*SINLONG*SINOBLIQ,1.0)
      SGB = ASIN( MAX(SINSGB,-1.0) )
      COSSGB = MAX( SQRT(1 - SINSGB**2) , 1.E-7 )
      COSLONG2 = MIN(COSLAT*COSLONG/COSSGB,1.0)
      LONG2 = ACOS( MAX(COSLONG2,-1.0) )
      SINLONG2 = (COSLAT*SINLONG*COSOBLIQ + SINLAT*SINOBLIQ)
     :               / COSSGB
      IF (SINLONG2 .GE. 0) THEN
            SGL = LONG2
         ELSE
            SGL = TWOPI - LONG2
      END IF
      IF (SGL .GE. TWOPI) SGL = SGL - TWOPI

      IF (CODE .EQ. 1 .OR. CODE .EQ. 3) THEN
         SGL = SGL/ONEDEG
         SGB = SGB/ONEDEG
      END IF

      END

C     ******************************************************************

      SUBROUTINE SGALTOGAL(SGL,SGB,L,B,CODE)

C     Converts (new) supergalactic coordinates (SGL,SGB), as given in
C     RC2 (de Vaucouleurs et al. 1976), to (new) galactic coordinates.

C     The SG north pole is defined in old galactic coordinates
C     (l=15,b=+5) in RC1 and RC2, but we assume here that the
C     approximate value in new galactic coordinates (l=47.37,b=+6.32)
C     given in RC2 is exactly correct. Note that the origin of
C     coordinates (new galactic coordinates) given in RC2 is incorrect;
C     the printed value l=137.29 is the new galactic latitude of the
C     OLD supergalactic origin. The correct origin is at l=137.37, b=0. 

C     CODE = 0: Input and output are in radians
C            1: Input in radians, output in hours and degrees
C            2: Input in degrees, output in radians
C            3: Input in degrees, output in hours and degrees

      REAL PI, ONEDEG, PIOV2, TWOPI, SINOBLIQ, COSOBLIQ, L0
      PARAMETER (PI=3.1415927)
      PARAMETER (ONEDEG=PI/180., PIOV2=PI/2, TWOPI=2*PI)
      PARAMETER (COSOBLIQ=0.1100813, SINOBLIQ=0.9939226)  

c                   83.68 deg

      PARAMETER (L0=ONEDEG*137.37)
      INTEGER CODE
      REAL L, B, SGL, SGB
      REAL LONG, LAT, SINLONG, COSLONG, SINLAT, COSLAT, SINB,
     :     COSB, LONG2, SINLONG2, COSLONG2

      IF (CODE .EQ. 0 .OR. CODE .EQ. 1) THEN
            LONG = SGL
            LAT = SGB
         ELSE IF (CODE .EQ. 2 .OR. CODE .EQ. 3) THEN
            LONG = ONEDEG*SGL
            LAT = ONEDEG*SGB
         ELSE
            STOP 'STOP: Bad code passed to GALTOEQU'
      END IF

      SINLAT = SIN(LAT)
      COSLAT = COS(LAT)
      SINLONG = SIN(LONG)
      COSLONG = COS(LONG)

      SINB = MIN(SINLAT*COSOBLIQ+COSLAT*SINLONG*SINOBLIQ,1.0)
      B = ASIN( MAX(SINB,-1.0) )
      COSB = MAX( SQRT(1 - SINB**2) , 1.E-7 )
      COSLONG2 = MIN(COSLAT*COSLONG/COSB,1.0)
      LONG2 = ACOS( MAX(COSLONG2,-1.0) )
      SINLONG2 = (COSLAT*SINLONG*COSOBLIQ - SINLAT*SINOBLIQ)
     :               / COSB
      IF (SINLONG2 .GE. 0) THEN
            L = LONG2 + L0
         ELSE
            L = (TWOPI-LONG2) + L0
      END IF
      IF (L .GE. TWOPI) L = L - TWOPI

      IF (CODE .EQ. 1 .OR. CODE .EQ. 3) THEN
         L = L/ONEDEG
         B = B/ONEDEG
      END IF

      END
