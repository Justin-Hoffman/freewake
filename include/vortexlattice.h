#ifndef VORTEXLATTICE_H
#define VORTEXLATTICE_H

#include <vector>
#include <utility>

#include "horseshoelattice.h"
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
        
        std::vector< std::vector<Vec3D> >& endPoints();
        std::vector< std::vector<Vec3D> >& endPointVelocity();
        std::vector< std::vector<double> >& gammaI();
        std::vector< std::vector<double> >& gammaJ();
        std::vector< std::vector<double> >& rcI();
        std::vector< std::vector<double> >& rcJ();

        void fixToTrailingEdge( HorseshoeLattice &h );
        void advect( double dt );
        void initializeToHelix( Vec3D axis, double dTheta, double dZ );
        void printState();
        
        virtual Vec3D calcInfluenceCoefficient( Vec3D p, int n);
        virtual Vec3D calcInducedVelocity( Vec3D );
        std::pair<int, int> ijFromN( int n );
        
        
    private:

        int ni_;
        int nj_;

        std::vector< std::vector<Vec3D> > endPoints_;
        std::vector< std::vector<Vec3D> > endPointV_;
        std::vector< std::vector<double> > gammaI_;
        std::vector< std::vector<double> > rcI_;
        std::vector< std::vector<double> > gammaJ_;
        std::vector< std::vector<double> > rcJ_;

};

#endif
