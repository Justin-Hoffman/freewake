#ifndef MATLABINTERFACE_H
#define MATLABINTERFACE_H

#include <vector>
#include <utility>

#include "liftingsurface.h"
#include "simulationmanager.h"
#include "vec3d.h"

#include "matrix.h"
#include "mex.h"

extern "C" 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

// A structure containing all data from the mxArray passed in as the RHS from MATLAB
struct MatlabInterfaceStruct {
    IntegrationScheme integrationScheme;

    double omega;                     /*!< Global rotation rate */ 
    double dt; 
    int nt;
    Vec3D globalLinearVelocity;
    Vec3D globalRotationAxis;

    double refL;
    double refC; 
    double refA;
    double refV;
    double vMach;

    int nSurfaces;
    int nChord; 
    int nSpan;
    int nNearWake; 
    int nFarWake; 


    double surfaceAR;
    double surfacePitch;
    double surfaceTipDihedral;
    double surfaceTipDihedralBreak;

    bool isFreeWake; 
    bool isFreeTipVortex;
    bool hasFixedTrailers;
    bool doPrandtlGlauert;
};

MatlabInterfaceStruct validateArgs( int nrhs, const mxArray *prhs[] );
#endif
