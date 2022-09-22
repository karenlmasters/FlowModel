c Useful parameters

      real onedeg
      parameter (onedeg=3.1415927/180.)
      real pi
      parameter (pi=3.1415927)

c DECLARATIONS
      CHARACTER name*8, infile*60
      INTEGER i
      REAL x

      PRINT *,' Enter name of input file to check'
      READ(*,'(a60)') INFILE
      OPEN(10,FILE=INFILE,STATUS='OLD')




C **************************************************
      SUBROUTINE CARTES(R1,R2,D1,D2,S)
C
C  Calculates separation between two objects
C     Input and output in degrees
C

      REAL COSSEP,R1,R2,D1,D2,S
      real onedeg
      parameter (onedeg = 0.017453292)

      r1 = r1*onedeg
      r2 = r2*onedeg
      d1 = d1*onedeg
      d2 = d2*onedeg

      COSSEP = SIN(D1)*SIN(D2) + 
     *         COS(D1)*COS(D2)*(COS(R1)*COS(R2) 
     *         + SIN(R1)*SIN(R2))
      COSSEP = DMIN1(COSSEP,1.0)
      COSSEP = DMAX1(COSSEP,-1.0)
      S = ACOS(COSSEP)
      s = s/onedeg

      r1 = r1/onedeg
      r2 = r2/onedeg
      d1 = d1/onedeg
      d2 = d2/onedeg

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c SUPERGALACTIC CO-ORDS
      sgx = r*cos(SGL)*cos(SGB)
      sgy = r*sin(SGL)*cos(SGB)
      sgz = r*sin(SGB)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c Subroutine writen May 15th 2003. KLM
      subroutine carttolbr(x,y,z,l,b,r)
c
c Input is x,y,z, output is r in same units and l and b in radians.
c SUPERGALACTIC x,y,z, l=SGL, b=SGB.      

      real x,y,z,l,b,r
      real onedeg
      parameter (onedeg=3.1415927/180.)
c Find r      
      r = sqrt(x**2 + y**2 + z**2)

c Find b
            if (z .ge. 0) then
               b = asin(z/r)
            elseif (z .eq. 0) then
               b = 0
            else 
               b = -asin(-z/r)
            endif
c Find l
            if (x .gt. 0 .and. y .ge. 0) then
                l = atan(y/x) 
            elseif (x .eq. 0 .and. y .ge. 0) then
               l = 90.*onedeg
            elseif (x .gt. 0 .and. y .lt. 0) then
               l = 360.*onedeg - atan(-y/x)
            elseif (x .eq. 0 .and. y .lt. 0) then
               l = 270.*onedeg
            elseif (x .lt. 0 .and. y .ge. 0) then
               l = 180.*onedeg - atan(-y/x)
            else
               l = 180.*onedeg + atan(y/x)
            endif

      return
      end

c Change RA and DEC to radians from hh mm ss, dd mm ss format.
c
      RA = 1.0*(irh+(irm + irs10/600.)/60.)     

      RAd = 15.0*RA

      if (SGN .eq. '+') then
         DEC = idd + (idm + ids/60.)/60. 
      elseif (SGN .eq. '-') then
         DEC = -(idd + (idm + ids/60.)/60.) 
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function chisq(x,y,sig,ndata)
c
c Subroutine to work out traditional chisq of a vector of data. 
c chisq = sum [(x - y)**2/sig**2]
c
c
 
      integer NMAX, NDATA
      parameter (NMAX=10000)
      real x(NMAX), y(NMAX), sig(NMAX)
      real chisq

      do 10 i=1,NDATA

         chisq = chisq + (x(i) - y(i))**2/sig(i)**2

 10   continue
      return
      end
