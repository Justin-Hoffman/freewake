#ifndef VORTEXLATTICE_H
#define VORTEXLATTICE_H

#include <vector>
#include <utility>

#include "vec3d.h"
#include "vortexcontainer.h"

class VortexLattice : public VortexContainer{
    public: 
        VortexLattice();
        VortexLattice( int ni, int nj );
        VortexLattice( const VortexLattice &vl );
        ~VortexLattice();

        int ni();
        int nj();
        
        std::vector<std::vector<Vec3D>>& getEndPoints();
        std::vector<std::vector<double>>& getGammaI();
        std::vector<std::vector<double>>& getGammaJ();
        
        virtual Vec3D calcInfluenceCoefficient( Vec3D p, int n);
        virtual Vec3D calcInducedVelocity( Vec3D );
        std::pair<int, int> ijFromN( int n );
        
        
    private:

        int ni_;
        int nj_;

        std::vector<std::vector<Vec3D>> endPoints_;
        std::vector<std::vector<Vec3D>> endPointV_;
        std::vector<std::vector<double>> gammaI;
        std::vector<std::vector<double>> gammaJ;

};

#endif
