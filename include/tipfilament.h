#ifndef TIPFILAMENT_H
#define TIPFILAMENT_H

#include <vector>
#include <utility>

#include "horseshoelattice.h"
#include "vec3d.h"
#include "vortexcontainer.h"
#include "vortexlattice.h"

class TipFilament : public VortexContainer{
    public: 
        TipFilament();
        TipFilament( int ni, int nj );
        TipFilament( const TipFilament &tf );
        ~TipFilament();

        int ni();
        int nj();
        
        std::vector< std::vector<Vec3D> >& endPoints();
        std::vector< std::vector<Vec3D> >& endPointVelocity();
        std::vector< std::vector<double> >& gamma();
        std::vector< std::vector<double> >& rc();

        void fixToWake( VortexLattice &h );
        void advect( double dt );
        void advectAndRotate( double dt, Vec3D axis, double omega );
        void advectPCC( double dt, Vec3D axis, double omega );
        void initializeToHelix( Vec3D axis, double dTheta, double dZ );
        void printState();
        
        virtual Vec3D calcInfluenceCoefficient( Vec3D p, int n);
        virtual Vec3D calcInducedVelocity( Vec3D, int jStart = 0 );
        std::pair<int, int> ijFromN( int n );
        
        
    private:

        int ni_;
        int nj_;

        std::vector< std::vector<Vec3D> > endPoints_;
        std::vector< std::vector<Vec3D> > endPointV_;
        std::vector< std::vector<double> > gamma_;
        std::vector< std::vector<double> > rc_;

};
#endif
