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

#endif
