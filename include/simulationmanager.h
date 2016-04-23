#ifndef SIMULATIONMANAGER_H
#define SIMULATIONMANAGER_H

#include <vector>
#include <utility>

#include "liftingsurface.h"

class SimulationManager {
    public:
        SimulationManager();
        SimulationManager( SimulationManager& ); 
        ~SimulationManager();
        SimulationManager & operator=( const SimulationManager& ); 
        
        void addSurface( LiftingSurface* s);
        void setGlobalLinearVelocity( Vec3D v );
        void solve();
        void integrateForceAndMoment();
        Vec3D netMoment(); 
        Vec3D netForce();
        double netLift();
        double netLiftDirect();
        double netDrag();

        Vec3D getGlobalLinearVelocity();
        int hijToN( int, int, int);
        int nInJToSuperN( int, int);
        std::tuple<int, int, int> hijFromN( int n);
        
    private: 
        std::vector<LiftingSurface*> surfaces_;
        std::vector<int>   nOffset_;
        double* a;
        double* b;
        Vec3D globalLinearVelocity_;
        Vec3D globalRotationAxis_;
        double globalRotationRate_;
        Vec3D netForce_;
        Vec3D netMoment_;
        Vec3D vInfinity( Vec3D p );
};
#endif
