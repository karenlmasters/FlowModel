CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine myflowmodelcall(CO1,CO2,COOPT,VEL,VELOPT,DISTANCE,
     *                         EDIST,FLAG)
c
c Subroutine to call flowmodel subroutine and deal with plausible errors on 
c distances and picking out of multiple versions. 
c
c Written by Karen L. Masters, based on flow model described in Masters (2005) - PhD 
c thesis. Please cite appropriately. 
c 
c In this version you can pick your favourite value of H0 at the end as all
c distances are dealt with in km/s units. Just divide the output distance 
c (which is in km/s) by H0.
c
c
c       FLAG    (INT)   Flag to indicate possible problems
c                       0: No distance found. Usually close to D=0Mpc
c                       1: Everything is fine, single valued distance
c                       2: Double valued distance
c                       3: Triple valued distance
c                       10: Assigned to Virgo Core
c                       11: Near Virgo (within 6Mpc which is region where SB00 
c                           claim that model is "uncertain"), single valued 
c                       12: Near Virgo (within 6Mpc which is region where SB00 
c                           claim that model is "uncertain"), double valued 
c                       13: Near Virgo (within 6Mpc which is region where SB00 
c                           claim that model is "uncertain"), triple valued 
c                       20: Assigned to GA Core
c                       21: Near GA (within 10Mpc which is region where SB00 
c                           claim that model is "uncertain"), single valued
c                       22: Near GA (within 10Mpc which is region where SB00 
c                           claim that model is "uncertain"), double valued
c                       23: Near GA (within 10Mpc which is region where SB00 
c                           claim that model is "uncertain"), triple valued
c                       30: Assigned to Core of Fornax cluster
c 
c Other outputs are as described in flowmodel.f subroutine.
c 
c Don't allow any variables which are no explicitly defined here.
c      implicit undefined(a-z)
c Parameters
      real onedeg
      parameter (onedeg=3.1415927/180.)
      real pi
      parameter (pi=3.1415927)
c NBRACK is the number of solutions the program will search for (input into 
c ZBRAK)
      integer NBRACK
      parameter (NBRACK=3)
c Variables
      integer velopt, coopt, NSOLNS, FLAG
      real co1, co2, ra, dec, sgl, sgb, l, b, r, vel, dist(NBRACK), 
     * vs(NBRACK), distance, edist, x, y, z
c
c Parameters of clusters
c
      real SBFh0, theta, veldiff, xv, yv, zv, dv, sigv, sglv, sgbv, vv, 
     * thetav, xg, yg, zg, sigg, sglg, sgbg, dg, vg, thetag, xf, yf, zf, 
     * sigf, sglf, sgbf, df, vf, thetaf, ratv, ratg, vcmb
c
c Different estimates of error on distance picked.
c       
      real diffdist, diffdist1, diffdist3
c
c Look for galaxies in the Virgo Cluster, the GA or Fornax in the model before
c bothering to find a distance for them. 
c NOTE: These are positions in model fit by T00.
c Need supergalactic co-ords for galaxy
c
c     Convert equatorial to sgal
      IF (COOPT .EQ. 0) THEN
         RA=CO1
         DEC=CO2
         CALL EQUTOSGAL(RA,DEC,SGL,SGB,0)
C     Convert gal to sgal
      ELSE IF (COOPT .EQ. 1) THEN
         L=CO1
         B=CO2
         CALL GALTOSGAL(L,B,SGL,SGB,0)
C     Or have sgal
      ELSE IF (COOPT .EQ. 2) THEN
         SGL=CO1
         SGB=CO2
      ELSE 
         pause 'Incorrect entry for COOPT'
      ENDIF      

c
c Convert velocity to CMB frame (using Lineweaver (1996) result of 
c V = 369.0+/-2.5 km/s in direction of l = 264.31+/-0.19 degrees and 
c b = 48.05 +/- 0.1 degrees: 
c deltaV = va{cos(ba)cos(la)cos(b)cos(l) + cos(ba)sin(la)cos(b)sin(l)+
c sin(ba)sin(b)} 
c      
      if (velopt .eq. 0) then
         vcmb = vel - 24.46*cos(b)*cos(l)-
     : 245.45*cos(b)*sin(l) + 274.44*sin(b)
      elseif (velopt .eq. 1) then
         Vcmb = vel
      else
         pause 'Incorrect entry for VELOPT'
      endif

c
c Search for galaxies within an angular core radius and with 
c velocities no more different from the cluster velocity than the velocity 
c dispersion.
 
c VIRGO CLUSTER
* VA parameters
      h0 = 70.
      xv = -4.7*h0
      yv = 15.9*h0
      zv = -0.3*h0
      call carttolbr(xv,yv,zv,sglv,sgbv,dv)
c CMB velocity of Virgo (assume at rest in CMB)
      vv = dv
      sigv = 650
c Angular size of 2Mpc core.
      thetav = 2.0*h0/dv
c
c If within angular distance of 2Mpc of core and velocity dispersion of 
c cluster velocity then set to distance of cluster. Goto the end of the 
c subroutine.
c
      theta1 = sqrt((SGL-sglv)**2 + (SGB-sgbv)**2)      
      veldiff = abs(vv - Vcmb)

c      write(*,*) theta, thetav, veldiff


      if (theta1 .le. thetav .and. veldiff .le. sigv) then
            NSOLNS = 1
            distance = dv
            EDIST = 2.0*h0
            FLAG=10
            goto 100
      endif

* GA parameters
      xg = -30.5*h0
      yg =  26.2*h0
      zg = -30.8*h0
      sigg = 500
      call carttolbr(xg,yg,zg,sglg,sgbg,dg)
c CMB velocity of GA (assume at rest in CMB)
      vg = dg
c Angular size of 2Mpc core.
      thetag = 2.0*h0/dg
c
c If within angular distance of 2Mpc of core and velocity dispersion of 
c cluster velocity then set to distance of cluster. Goto the end of the 
c subroutine.
c
      theta2 = sqrt((SGL-sglg)**2 + (SGB-sgbg)**2)      
      veldiff = abs(vg - Vcmb)

      if (theta2 .le. thetag .and. veldiff .le. sigg) then
            NSOLNS = 1
            distance = dg
            EDIST = 2.0*h0
            FLAG=20
            goto 100
      endif

c
c First pass at distance.
c 
      call myflowmodel(CO1,CO2,COOPT,VEL,VELOPT,DIST,VS,NSOLNS,FLAG)
c
c Galaxies which the model just has one possible distance
c
      if (NSOLNS .eq. 1) then
         DISTANCE = DIST(1)
         EDIST = VS(1)
c     
c Galaxies with NSOLNS = 0 are usually at d=0Mpc (or very close by). Set 
c error equal to thermal velocity dispersion divided by h0.
c
      elseif (NSOLNS .eq. 0) then
         DISTANCE = 0.0
         Edist = 163.
c      
c Two solns. Pick Dist(1). Error is max of difference between the two distances
c or error on DIST(1) using above method
c
      elseif (NSOLNS .eq. 2) then
         Diffdist = abs(DIST(1) - DIST(2))

         EDIST = VS(1)
         EDIST = amax1(EDIST,diffdist)
c
c Three solns. Pick DIST(2). Error is max of difference between the three 
c distances or error on DIST(2) using above method.
c
      elseif (NSOLNS .eq. 3) then
        distance = DIST(2)
        Diffdist1 = abs(DIST(1) - DIST(2))         
        Diffdist3 = abs(DIST(2) - DIST(3))
        Diffdist = (Diffdist1 + Diffdist3)/2.0
        
        EDIST = VS(2)
        EDIST = amax1(EDIST,diffdist)
c
c Something wierd happened.
c
      else
         write(*,*) 'Don"t know how to deal with NSOLNS=',Nsolns
      endif

c
c Search for galaxies near Virgo or GA. 
c
c Need supergalactic x, y ,z.
 
      r=distance
      x = r*cos(SGL)*cos(SGB)
      y = r*sin(SGL)*cos(SGB)
      z = r*sin(SGB)

      ratv = sqrt((x-xv)**2+(y-yv)**2+(z-zv)**2)
      if (ratv .le. 6.0*h0 .or. theta1 .le. thetav) then
         FLAG = FLAG+10
      endif

      write(*,*) ratv

* GA parameters
  
      ratg = sqrt((x-xg)**2+(y-yg)**2+(z-zg)**2)
      if (ratg .le. 10.0*h0 .or. theta2 .le. thetag) then
         FLAG = FLAG+20
      endif      

c Sent to here if assigned to a cluster core.
c
 100  continue
 50   return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine myflowmodel(CO1,CO2,COOPT,VEL,VELOPT,DIST,VS,NSOLNS,
     *                     FLAG)
c
c Compute distances from galaxy velocities and co-ordinates in the non-linear
c flow model of Tonry et al (2000) (uses Yahill (1985) approximation which is
c good to an overdensity of 30). This model includes the spherical infall 
c on the Virgo Attractor and the Great Attractor, as well as a quadrupole 
c correction which could be due either to non-spherical infall or the 
c influence of external masses. There is also a velocity dispersion component 
c which is a combination of thermal velocity dispersion and virial motions in 
c Virgo, the GA and the Fornax Cluster. The parameters of the model were fit 
c to a sample of about 300 early-type galaxies with SBF distances (see Tonry 
c et al 2000). 
c The output of this subroutine consists of up to three possible best 
c solutions, a file from which you can make a plot of the unnormalized 
c probability given the model that the galaxy is at a given distance (which 
c allows an estimate of the error in the model on the distance from the 
c velocity dispersion) and also a file from which you can plot the function
c mydeltav(dist) for the galaxy (to check the roots found are sensible).
c
c Created 24th Jan 2003, Karen L. Masters (KLM)
c Last update 4th Feb 2003, KLM 
c
c This version provides distances independent of h0 (ie. in units of km/s).
c May 16th 2003, KLM.
c J2000 co-ords, Feb 2005, KLM
c Uses flow model fit to SFI++. June 2005, KLM
c
c INPUT:
c
C       CO1,CO2 (REAL)  Coordinates of galaxy
C       COOPT   (INT)   0: Input RA, DEC in radians (J2000)
C                       1: Input L,B in radians
C                       2: Input SGL,SGB in radians
C
c
C       VEL     (REAL)  Radial velocity of source (km/s)
C       VELOPT  (INT)   0: Input heliocentric velocity
C                       1: Input CMB velocity
C
c OUTPUT:
c
C       DIST    (REAL)  3-dimensional vector containing NB distance
C                       solutions in km/s
C       NSOLNS  (INT)   Number of solutions found
c
c       FLAG    (INT)   Flag to indicate if in/near a cluster
c                       0: Everything is fine
c                       10: Assigned to Virgo Core
c                       11: Near Virgo (within 6*78.4km/s which is region 
c                           where SB00 
c                           claim that model is "uncertain")
c                       20: Assigned to GA Core
c                       21: Near GA (within 10*78.4km/s which is region where 
c                           SB00 
c                           claim that model is "uncertain")
c                       30: Assigned to Core of Fornax cluster
c
c SUBROUTINES CALLED:
c
c     sbf2flow(x,y,z,vx,vy,vz,vr,vsig) 
c              Model for velocity field
c              From www.ifa.hawaii.edu/~jt/SOFT/sbf2flow.f 
c              Discussed in Tonry et al (2000) ApJ 530, 625
c              (and included below)
c
c     equtogal, equtosgal, galtosgal, galtoequ, sgaltogal, sgaltoequ   
c              Co-ordinate transformations. 
c              In /home/esperanza3/galaxy/gallib/source/coords.f
c
c     zbrak,zbrent    
c              Root finding from Numerical Recipes (zbrent is a function)
c
c     mydeltav   Function giving difference between observed and predicted 
c              velocity at a given distance. 
c 
c     carttolbr 
c              Find lattitude, longitude and distance from cartesian 
c              co-ords.
c
c PARAMETERS IN SBF2FLOW:
c
c MODEL:
c
c Cosmology:
c       Ho             Hubble's constant
c       (wx,wy,wz)     Peculiar velocity of Sun in CMB frame
c       omega          Mass density parameter
c       thermal        Thermal velocity dispersion in universe
c       
c Quadrupole:
c       rquad          Cut-off radius
c       qxx, qxy, qxz, qyz, qzz
c
c Attractors:
c       xa,ya,za       Position in Supergalactic cartesian co-ords relative
c                      to Local Group
c       deltaa         Overdensity at position of LG
c       gammaa         Power law fall off
c       rcore          Core radius
c       rcuta          Cut-off radius
c       siga           Virial velocity dispersion which falls off like a 
c                      Gaussian centred on the attractor with width rcore
c     
c           Note: model also has non-attracting overdensities which are
c                 simply described by a position and virial velocity dispersion
c
c Input to model
c 
c       sgx,sgy,sgz     Supergalactic cartesian co-ords of galaxy
c       r               Trial distance for galaxy      
c
c Output from model
c
c       vx,vy,vz        Velocity in supergalactic cartesian co-ords in 
c                       CMB frame
c       vr              Radial velocity in CMB frame
c       vsig            Velocity dispersion at position of galaxy
c

c Don't allow any variables which are no explicitly defined here.
c      implicit undefined(a-z)
c Parameters
      real onedeg
      parameter (onedeg=3.1415927/180.)
      real pi
      parameter (pi=3.1415927)
c NBRACK is the number of solutions the program will search for (input into 
c ZBRAK)
      integer NBRACK
      parameter (NBRACK=3)
c TOL is the accuracy to which the solutions will be found (input into ZBRENT)
      real TOL
      parameter (TOL=1.0E-4)
c Variables
      integer velopt, coopt, NSOLNS, i, FLAG
      real co1, co2, vel, ra, dec, l, b, dist(NBRACK), vs(NBRACK),
     * r, sgx, sgy, sgz, vx,vy,vz,vr,vsig,prob,rmin,rmax,xb1(NBRACK),
     * xb2(NBRACK), vcmb, del
c Variables and parameters common to FLOWMODEL, MYDELTAV and ZBRAK
      common /galaxy/ sgl, sgb, vgal
      real sgl, sgb, vgal
c Functions
      real mydeltav,zbrent
      external mydeltav,zbrent

c
c Convert co-ordinates into galactic and supergalactic, for corrections and
c input into SBF2FLOW
c
c     Convert equatorial to sgal and gal
      IF (COOPT .EQ. 0) THEN
         RA=CO1
         DEC=CO2
         CALL EQU2000TOSGAL(RA,DEC,SGL,SGB,0)
         CALL EQU2000TOGAL(RA,DEC,L,B,0)
C     Convert gal to sgal and equ
      ELSE IF (COOPT .EQ. 1) THEN
         L=CO1
         B=CO2
         CALL GALTOSGAL(L,B,SGL,SGB,0)
         CALL GALTOEQU2000(L,B,RA,DEC,0)
C     Or convert sgal to gal and equ
      ELSE IF (COOPT .EQ. 2) THEN
         SGL=CO1
         SGB=CO2
         CALL SGALTOGAL(SGL,SGB,L,B,0)
         CALL SGALTOEQU2000(SGL,SGB,RA,DEC,0)
      ELSE 
         pause 'Incorrect entry for COOPT'
      ENDIF
c
c Convert velocity to CMB frame (using Lineweaver (1996) result of 
c V = 369.0+/-2.5 km/s in direction of l = 264.31+/-0.19 degrees and 
c b = 48.05 +/- 0.1 degrees: 
c deltaV = va{cos(ba)cos(la)cos(b)cos(l) + cos(ba)sin(la)cos(b)sin(l)+
c sin(ba)sin(b)} 
c      
      if (velopt .eq. 0) then
         vgal = vel - 24.46*cos(b)*cos(l)-
     : 245.45*cos(b)*sin(l) + 274.44*sin(b)
         Vcmb = vgal
      elseif (velopt .eq. 1) then
         vgal = vel 
         Vcmb = vel
      else
         pause 'Incorrect entry for VELOPT'
      endif
c
c
c Plausible range of distances. Peculiar velocity in excess of 2000 would be 
c extremely unusual. I 
c made this large so as to not miss any possible distances, but there is no 
c point in it being the whole universe.
c
      rmin = MAX((vgal-2000),0.0)
      rmax = vgal + 2000
c
c Calculate probability distribution for distances given the velocity and 
c position. P = (1.0/sqrt(2*pi)*sigv)*exp(-0.5*((vgal-vmod)/vsig)**2), from 
c Tonry et al. (2000)
c 
c We do this between dist=rmin and rmax Mpc in intervals of 0.5/h Mpc
c
c File to write out to
c      write(*,*) 'Probability distribution in distprob.out'
c      open(1,file='distprob.out',status='unknown')
c      write(1,*) 'Prob. of given distance (Mpc) for galaxy at SGL=',SGL,
c     * ', SGB=',SGB,', Vcmb=',vgal
c      write(*,*) 'Function to find roots of output to deltav.out'
c      open(2,file='deltav.out',status='unknown')
c      write(2,*) 'Observed minus predicted velocity at given distances f
c     *or galaxy at SGL=',SGL,', SGB=',SGB,', Vcmb=',vgal

      r = rmin
      do while (r .le. rmax) 
c Calculate supergalactic cartesian co-ords of galaxy which are the input to
c the SBF2FLOW subroutine 
         sgx = r*cos(SGL)*cos(SGB)
         sgy = r*sin(SGL)*cos(SGB)
         sgz = r*sin(SGB)
c Call SBF2FLOW subroutine 
         call klmflow(sgx,sgy,sgz,vx,vy,vz,vr,vsig)
c Calcluate probability of distance given observed velocity vgal
         prob = (1000.0/(sqrt(2*pi)*vsig))*exp(-0.5*((vgal-vr)/vsig)**2)
c Write out to a file
c         write(1,'(2F10.3)') r, prob
c Calculate function DELTAV at the same intervals
         del = mydeltav(r)
c Write out to a file (since it is good practise to look at the function
c you are finding roots of
c         write(2,'(2F10.3)') r, del
c Step size of 0.5 Mpc
         r = r+50.
      enddo
c
c Search for distances at which model predicts observed CMB velocities
c
c Call ZBRAK (Numerical Recipes) to look for NBRAK=3 zero's in deltav
      NSOLNS = NBRACK
c NSOLNS will be reset to the actual number of solutions found. deltav is 
c a function which computes the difference between the observed and predicted 
c velocity at a given distance. Interval is divided into 50 segments to 
c search for zero crossings. xb1(NSOLNS), and xb2(NSOLNS) are output as the 
c NSOLNS bracketing pairs.
      call zbrak(mydeltav,rmin,rmax,50,xb1,xb2,NSOLNS)
      FLAG = NSOLNS
c Use function ZBRENT (Numerical Recipes) which searches for the root in the 
c intervals found by ZBRAK. Root is returned as zbrent with accuracy tol=1.E-4
c Algorithm is 'guaranteed' by Brent to converge!
      do 10 i=1,NSOLNS
         rmin = xb1(i)
         rmax = xb2(i)
         dist(i) = zbrent(mydeltav,rmin,rmax,tol)
 10      continue

c
c Want velocity dispersion at possible distances to get an estimate of the 
c errors
c

       do 11 i=1,NSOLNS
         sgx = dist(i)*cos(SGL)*cos(SGB)
         sgy = dist(i)*sin(SGL)*cos(SGB)
         sgz = dist(i)*sin(SGB)         
         call klmflow(sgx,sgy,sgz,vx,vy,vz,vr,vs(i))
 11      continue

      return
      end

c      include 'zbrack.f'
c      include 'zbrent.f'
c      include 'carttolbr.f'
c      include 'coords2000.f'
c      include 'coord.f'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      real function mydeltav(dist)
c
c Difference of observed velocity and that predicted by SBF2FLOW at a given 
c distance. This will be minimized by ZBRAK and ZBRENT in subroutine 
c FLOWMODEL above. KLM Jan 2003.
c
c      implicit undefined(a-z)
      real dist, sgx, sgy, sgz, vr, vx, vy, vz, vsig
c Variables and parameters common to FLOWMODEL and DELTAV
      common /galaxy/ sgl, sgb, vgal
      real sgl, sgb, vgal
c Calculate supergalactic cartesian co-ords of galaxy at distance=dist
      sgx = dist*cos(SGL)*cos(SGB)
      sgy = dist*sin(SGL)*cos(SGB)
      sgz = dist*sin(SGB)
c Call model subroutine 
      call klmflow(sgx,sgy,sgz,vx,vy,vz,vr,vsig)      

      mydeltav = vgal - vr

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine klmflow(x,y,z,vx,vy,vz,vr,vsig)
C
C Based on subroutine written by John Tonry (available on-line).
C Edited to deal with distances in km/s (independent of Hubble's constant).
C This version also uses NFW profiles for the attractors, and does not include
C cutoff radii, except for the quadrupole.
C

      real x,y,z,vx,vy,vz,vr,vsig
c Parameters
      real h0, omega, wx, wy, wz
      real rcore, thermal
      real gammav, sigv, rcutv, xv, yv, zv
      real gammag, sigg, rcutg, xg, yg, zg
      real sigf, xf, yf, zf
      real rquad, qxx, qxy, qxz, qyz, qzz
c Stuff in program
      integer ivrh, ivrv, ivrg, ivrq
      real d3, delta0, delta, uinfall
      real distv, ratv
      real distg, ratg
      real ratq, cut, qyy
      real ratf

      ezp(x) = exp(amax1(-20.0,x))

      if(x.eq.0.and.y.eq.0.and.z.eq.0) x = 0.001
c FIXED PARAMETERS
* Cosmology parameters
      h0 = 70.
      omega = 0.3
* Common parameters
c NFW concentration parameter.
      c=5.
      deltac = 200.*c**3/(3.*(log(1+c) - c/(1.+c)))
c Thermal velocity dispersion
      rcore = 2.*h0
      thermal = 163.
* VA parameters
      sigv = 650
* GA parameters
      sigg = 500
* Quadrupole parameters
      rquad = 50*h0

c VARYING PARAMETERS
* Cosmology parameters
      wx = -204
      wy =  160
      wz = -261
* VA parameters
      xv = -4.73*h0
      yv = 15.87*h0
      zv = 0.32*h0
      r200v = 1.7*h0
* GA parameters
      xg=-30.51*h0
      yg=26.18*h0
      zg=-30.79*h0
      r200g = 2.17*h0
* Quadrupole parameters
      qxx =  8.52
      qxy =  0.51
      qxz = -2.71
      qyz =  0.36
      qzz = -8.81

**************************************************

* Peculiar velocity
      vx = wx
      vy = wy
      vz = wz

      ivrw = nint((vx*x + vy*y + vz*z) / sqrt(x*x+y*y+z*z))

* Hubble flow contribution
      vx = vx + x
      vy = vy + y
      vz = vz + z

      ivrh = nint(sqrt(x*x+y*y+z*z))

* Virgo Attractor
      distv = amax1(0.001, sqrt(xv**2+yv**2+zv**2))
c      r200v = 0.0207*(Mv)**(1./3.)/1000.
      rcv0 = c*distv/r200v
      delta0 = deltac*(1./(1+rcv0) + log(1+rcv0) - 1)/(omega*rcv0**3)

      ratv = amax1(0.001, sqrt((x-xv)**2+(y-yv)**2+(z-zv)**2))
      rcv = c*ratv/r200v
      delta = deltac*(1./(1+rcv) + log(1+rcv) - 1)/(omega*rcv**3)

      uinfall = 1./3. * ratv * omega**0.6 * delta * (1+delta)**-0.25
      
      ux = -uinfall * (x-xv)/ratv
      uy = -uinfall * (y-yv)/ratv
      uz = -uinfall * (z-zv)/ratv

      vx = vx + ux
      vy = vy + uy
      vz = vz + uz

      ivrv = nint((ux*x + uy*y + uz*z) / sqrt(x*x+y*y+z*z))

* Great Attractor
      distg = amax1(0.001, sqrt(xg**2+yg**2+zg**2))
c      r200g = 0.0207*(Mg)**(1./3.)/1000.
      rcg0 = c*distg/r200g
      delta0 = deltac*(1./(1+rcg0) + log(1+rcg0) - 1)/(omega*rcg0**3)

      ratg = amax1(0.001, sqrt((x-xg)**2+(y-yg)**2+(z-zg)**2))
      rcg = c*ratg/r200g
      delta = deltac*(1./(1+rcg) + log(1+rcg) - 1)/(omega*rcg**3)

      uinfall = 1./3. * ratg * omega**0.6 * delta * (1+delta)**-0.25
      
      ux = -uinfall * (x-xg)/ratg
      uy = -uinfall * (y-yg)/ratg
      uz = -uinfall * (z-zg)/ratg

      vx = vx + ux
      vy = vy + uy
      vz = vz + uz

      ivrg = nint((ux*x + uy*y + uz*z) / sqrt(x*x+y*y+z*z))

* Quadrupole
      ratq = sqrt(x*x+y*y+z*z)
      cut = ezp(-0.5*ratq*ratq/(rquad*rquad))
      qyy = -qxx - qzz
      ux = cut*(x*qxx + y*qxy + z*qxz)/h0
      uy = cut*(x*qxy + y*qyy + z*qyz)/h0
      uz = cut*(x*qxz + y*qyz + z*qzz)/h0
c
      vx = vx + ux
      vy = vy + uy
      vz = vz + uz
c
      ivrq = nint((ux*x + uy*y + uz*z) / sqrt(x*x+y*y+z*z))

* Radial component
      vr = (vx*x + vy*y + vz*z) / sqrt(x*x+y*y+z*z)

* Thermal velocity
      vsig = thermal*thermal
      vsig = vsig + sigv*sigv*ezp(-ratv*ratv/(rcore*rcore))
      vsig = vsig + sigg*sigg*ezp(-ratg*ratg/(rcore*rcore))
      vsig = sqrt(vsig)

* Details (if desired)
c      write(*,*) ivrh, ivrw, ivrv, ivrg, ivrq, nint(vr), nint(vsig)

 99   return
      end


