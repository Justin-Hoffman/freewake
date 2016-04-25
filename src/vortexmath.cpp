#include <cmath>
#include <limits>

#include "stdio.h"
#include "vortexmath.h"


Vec3D BiotSavart( Vec3D &r1, Vec3D &r2, double rc ){
    /*Vec3D ds = (r2-r1);
     Vec3D r = (r2+r1)/2.0;
    double mag_r = r.magnitude();
    
    return ( ds.cross(r) / ( -4.0 * M_PI * ( mag_r * mag_r * mag_r + rc*rc ) ) );
    */
    Vec3D r0 = r2-r1;
    r1 = -1.0*r1;
    r2 = -1.0*r2;
    
    Vec3D r1r2 = r1.cross(r2);
    double mr1r2 = r1r2.magnitude();

    return r1r2/( 4.0 * M_PI * mr1r2 * mr1r2 + rc*rc ) * ( r0.dot( r1.norm()-r2.norm() ) );
    
}

