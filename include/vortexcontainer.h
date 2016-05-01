#ifndef VORTEXCONTAINER_H
#define VORTEXCONTAINER_H
#include "vec3d.h"

class VortexContainer {
    public:
        VortexContainer(){};
        virtual ~VortexContainer(){};
 
        virtual Vec3D calcInfluenceCoefficient( Vec3D p, int n) = 0;
        virtual Vec3D calcInducedVelocity( Vec3D, int jStart ) = 0;
};

#endif
