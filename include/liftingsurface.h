#ifndef LIFTINGSURFACE_H
#define LIFTINGSURFACE_H

#include <array>
#include <stdint.h>
#include <tuple>

#include "horseshoelattice.h"
#include "tipfilament.h"
#include "vortexlattice.h"

class VortexLattice;
class LiftingSurface{
    public:
        LiftingSurface();         //!< 
        LiftingSurface( const LiftingSurface& );         //!< 
        LiftingSurface( int nSpan, int nChord );         //!< 
        LiftingSurface( int nSpan, int nChord, int nWake, int nFilament );         //!< 
        ~LiftingSurface();        //!< 
        
        HorseshoeLattice& getHorseshoeLattice();
        VortexLattice& getVortexLattice();
        TipFilament& getTipFilament();
        
        bool freeWake();   
        bool freeTipVortex();   
     
        double getAspectRatio();
        double getCoreRadius();
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
        int nFilament();
        
        void setFreeWake( bool );   
        void setFreeTipVortex( bool );   

        void setAspectRatio( double );
        void setCoreRadius( double );
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
        bool freeTipVortex_;
        double span_;
        double sweep_;
        double pitch_;
        double rc_;
        double rootChord_;
        double tipChord_;
        double tipDihedral_;
        double tipDihedralBreak_;
        double taperRatio_;
        
        int nSpan_;
        int nChord_;
        int nWake_;
        int nFilament_;
        
        std::vector< Vec3D > spanwiseForce_;
        
        HorseshoeLattice horseshoeLattice_;    
        VortexLattice vortexLattice_;    
        TipFilament tipFilament_;    
};

struct ReferenceSurface {
    double S;
    double b;
    double cbar;
    double vMach; 
    bool pgCorrection;
    ReferenceSurface() : S( 1.0 ), b( 1.0 ), cbar(1.0), vMach(999.999), pgCorrection(false) {}   
    ReferenceSurface( double s_in, double b_in, double c_in ) : S(s_in), b(b_in), cbar(c_in), vMach(999.99), pgCorrection(false){}
};
#endif
