CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      real function deltav(dist)
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
      call wb01flow(sgx,sgy,sgz,vx,vy,vz,vr,vsig)      

      deltav = vgal - vr

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
