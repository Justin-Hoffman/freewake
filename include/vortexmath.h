#ifndef VORTEXMATH_H
#define VORTEXMATH_H

#include "vec3d.h"

Vec3D BiotSavart( Vec3D &r1, Vec3D &r2, double rc );

double VortexCoreGrowth( double rc_old, double dt );

#endif
