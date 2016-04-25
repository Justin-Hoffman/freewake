#ifndef LIFTINGSURFACE_H
#define LIFTINGSURFACE_H

#include <array>
#include <stdint.h>
#include <tuple>

#include "horseshoelattice.h"

class LiftingSurface{
    public:
        LiftingSurface();         //!< 
        LiftingSurface( const LiftingSurface& );         //!< 
        LiftingSurface( int nSpan, int nChord );         //!< 
        ~LiftingSurface();        //!< 
        
        HorseshoeLattice& getLattice();
        
        double getAspectRatio();
        double getSpan();
        double getSweep();
        double getPitch();
        double getRootChord();
        double getTipChord();
        double getTaperRatio();
        double getTipDihedral();
        double getTipDihedralBreak();
        
        int nSpan();
        int nChord();
        
        void setAspectRatio( double );
        void setSpan( double );
        void setSweep( double );
        void setPitch( double );
        void setRootChord( double );
        void setTipChord( double );
        void setTaperRatio( double );
        void setTipDihedral( double );
        void setTipDihedralBreak( double );

        void updateLattice();
               
    private:
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

        HorseshoeLattice lattice_;    
};

struct ReferenceSurface {
    double S;
    double b;
    double cbar;
    ReferenceSurface() : S( 1.0 ), b( 1.0 ), cbar(1.0) {}   
    ReferenceSurface( double s_in, double b_in, double c_in ) : S(s_in), b(b_in), cbar(c_in){}
};
#endif
