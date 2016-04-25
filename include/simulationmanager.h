#ifndef SIMULATIONMANAGER_H
#define SIMULATIONMANAGER_H

#include <vector>
#include <utility>

#include "liftingsurface.h"

struct ForcesAndMoments { //Container forces and moments in body coordinates (X,Y,Z) and aero coordinates(Drag,Sideforce,Lift)
    Vec3D bodyForce;
    Vec3D bodyMoment;
    Vec3D bodyForceCoeff;
    Vec3D bodyMomentCoeff;
    Vec3D aeroForce;
    Vec3D aeroMoment;
    Vec3D aeroForceCoeff;
    Vec3D aeroMomentCoeff;
    ForcesAndMoments() : bodyForce(), bodyMoment(), bodyForceCoeff(), bodyMomentCoeff(), aeroForce(), aeroMoment(), aeroForceCoeff(), aeroMomentCoeff(){}
};


class SimulationManager {
    public:
        SimulationManager();
        SimulationManager( SimulationManager& ); 
        ~SimulationManager();
        SimulationManager & operator=( const SimulationManager& ); 
        
        void addSurface( LiftingSurface* s);
        void setGlobalLinearVelocity( Vec3D v );
        void solve();
        void setReferenceSurface( ReferenceSurface );
        void setReferenceVelocity( double );
        void integrateForceAndMoment();
        void printState();

        Vec3D netMoment(); 
        Vec3D netForce();
        ForcesAndMoments forcesAndMoments();
        double netLift();
        double netLiftDirect();
        double netDrag();
        ReferenceSurface referenceSurface();
        double referenceVelocity();
        Vec3D getGlobalLinearVelocity();
        int hijToN( int, int, int);
        int nInJToSuperN( int, int);
        std::tuple<int, int, int> hijFromN( int );
        
    private:
        bool needsSolve_; 
        std::vector<LiftingSurface*> surfaces_;
        std::vector<int>   nOffset_;
        double* a;
        double* b;
        ReferenceSurface refSurf_;
        double refV_;
        Vec3D globalLinearVelocity_;
        Vec3D globalRotationAxis_;
        double globalRotationRate_;
        ForcesAndMoments fomo_; 

        Vec3D vInfinity( Vec3D p );
};

#endif
