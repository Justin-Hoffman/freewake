#ifndef SIMULATIONMANAGER_H
#define SIMULATIONMANAGER_H

#include <vector>
#include <utility>

#include "liftingsurface.h"
#include "vortexlattice.h"


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
        SimulationManager( const SimulationManager& ); 
        ~SimulationManager();
        
        void addSurface( LiftingSurface* s);
        void advectWake();
        void calculateWakeVelocities();
        void fillWakeBC();
        void setGlobalLinearVelocity( Vec3D v );
        void setGlobalRotationAxis( Vec3D v );
        void setGlobalRotationRate( double );
        void solve();
        void setReferenceSurface( ReferenceSurface );
        void setReferenceVelocity( double );
        void setDt( double );
        void step();
        void integrateForceAndMoment();
        void printState();

        Vec3D netMoment(); 
        Vec3D netForce();
        ForcesAndMoments forcesAndMoments();
        double netLift();
        double netLiftDirect();
        double netDrag();
        double dt();
        ReferenceSurface referenceSurface();
        double referenceVelocity();
        LiftingSurface& getSurface( int );
        int getNSurfaces();
        Vec3D getGlobalLinearVelocity();
        Vec3D getGlobalRotationAxis();
        double getGlobalRotationRate();
        int hijToN( int, int, int);
        int nInJToSuperN( int, int);
        std::tuple<int, int, int> hijFromN( int );
        
    private:
        bool needsSolve_; 
        std::vector<LiftingSurface*> surfaces_;
        std::vector<int>   nOffset_;
        std::vector< std::vector<double> > lastGamma_;
        std::vector< std::vector<double> > thisGamma_;
        ReferenceSurface refSurf_;
        double refV_;
        double dt_;
        Vec3D globalLinearVelocity_;
        Vec3D globalRotationAxis_;
        double globalRotationRate_;
        ForcesAndMoments fomo_; 
        
        void fillRHS( double* );
        void fillLHS( double* );
        
        
        Vec3D vInfinity( Vec3D p );
};

#endif
