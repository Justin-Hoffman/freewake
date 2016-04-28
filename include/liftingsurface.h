#ifndef LIFTINGSURFACE_H
#define LIFTINGSURFACE_H

#include <array>
#include <stdint.h>
#include <tuple>

#include "horseshoelattice.h"
#include "vortexlattice.h"

class VortexLattice;
class LiftingSurface{
    public:
        LiftingSurface();         //!< 
        LiftingSurface( const LiftingSurface& );         //!< 
        LiftingSurface( int nSpan, int nChord );         //!< 
        LiftingSurface( int nSpan, int nChord, int nWake );         //!< 
        ~LiftingSurface();        //!< 
        
        HorseshoeLattice& getHorseshoeLattice();
        VortexLattice& getVortexLattice();
        
        bool freeWake();   
     
        double getAspectRatio();
        double getSpan();
        double getSweep();
        double getPitch();
        double getRootChord();
        double getTipChord();
        double getTaperRatio();
        double getTipDihedral();
        double getTipDihedralBreak();

        std::vector< Vec3D >& spanwiseForce();
        
        Vec3D calcInducedVelocity( Vec3D );
    
        int nSpan();
        int nChord();
        int nWake();
        
        void setFreeWake( bool );   

        void setAspectRatio( double );
        void setSpan( double );
        void setSweep( double );
        void setPitch( double );
        void setRootChord( double );
        void setTipChord( double );
        void setTaperRatio( double );
        void setTipDihedral( double );
        void setTipDihedralBreak( double );
        void setVortexLattice( VortexLattice &v );

        void updateLattice();
               
    private:
        bool freeWake_;
        double span_;
        double sweep_;
        double pitch_;
        double rootChord_;
        double tipChord_;
        double tipDihedral_;
        double tipDihedralBreak_;
        double taperRatio_;
        
        int nSpan_;
        int nChord_;
        int nWake_;
        
        std::vector< Vec3D > spanwiseForce_;
        
        HorseshoeLattice horseshoeLattice_;    
        VortexLattice vortexLattice_;    
};

struct ReferenceSurface {
    double S;
    double b;
    double cbar;
    ReferenceSurface() : S( 1.0 ), b( 1.0 ), cbar(1.0) {}   
    ReferenceSurface( double s_in, double b_in, double c_in ) : S(s_in), b(b_in), cbar(c_in){}
};
#endif
