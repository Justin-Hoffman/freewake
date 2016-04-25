#include "liftingsurface.h"

LiftingSurface::LiftingSurface() : 
                span_( 1.0 ), sweep_( 0.0 ), pitch_(0.0), rootChord_( 1.0 ), tipChord_(1.0), tipDihedral_(0.0), tipDihedralBreak_(0.0), taperRatio_(1.0), nSpan_( 10 ), nChord_( 2 ),
                lattice_( nSpan_, nChord_) {
    lattice_.snapToUnit();
}

LiftingSurface::LiftingSurface( const LiftingSurface &l) : 
                span_( l.span_ ), sweep_( l.sweep_ ), pitch_( l.pitch_ ), rootChord_( l.rootChord_ ), tipChord_( l.tipChord_ ), tipDihedral_( l.tipDihedral_ ), 
                tipDihedralBreak_( l.tipDihedralBreak_ ), taperRatio_( l.taperRatio_ ), nSpan_( l.nSpan_ ), nChord_( l.nChord_ ),
                lattice_( l.lattice_ ) {
}

LiftingSurface::LiftingSurface( int nSpan, int nChord) :
                span_( 1.0 ), sweep_( 0.0 ), pitch_(0.0), rootChord_( 1.0 ), tipChord_(1.0), tipDihedral_(0.0), tipDihedralBreak_(0.0), taperRatio_(1.0), nSpan_( nSpan ), nChord_( nChord ),
                lattice_( nSpan_, nChord_ ) {
    lattice_.snapToUnit();
}        

LiftingSurface::~LiftingSurface(){
}

HorseshoeLattice& LiftingSurface::getLattice(){
    return lattice_;
}

int LiftingSurface::nSpan(){
    return nSpan_;
}

int LiftingSurface::nChord(){
    return nChord_;
}

double LiftingSurface::getAspectRatio(){
    return 2.0*span_ / ( rootChord_ + tipChord_ );
}

void LiftingSurface::setAspectRatio( double AR ){
    span_ =  AR*( rootChord_ + tipChord_ ) / 2.0;
}

double LiftingSurface::getSpan(){
    return span_;
}
void LiftingSurface::setSpan( double span ){
    span_ = span;
}

double LiftingSurface::getSweep(){
    return sweep_;
}
void LiftingSurface::setSweep( double sweep ){
    sweep_ = sweep;
}

double LiftingSurface::getPitch(){
    return pitch_;
}
void LiftingSurface::setPitch( double pitch ){
    pitch_ = pitch;
}


double LiftingSurface::getRootChord(){
    return rootChord_;
}

void LiftingSurface::setRootChord( double rootChord ){
    rootChord_ = rootChord;
    taperRatio_ = tipChord_/rootChord_;
}

double LiftingSurface::getTipChord(){
    return tipChord_;
}

void LiftingSurface::setTipChord( double tipChord ){
    tipChord_ = tipChord;
    taperRatio_ = tipChord_/rootChord_;
}

double LiftingSurface::getTipDihedral(){
    return tipDihedral_;
}

void LiftingSurface::setTipDihedral( double tipDihedral ){
    tipDihedral_ = tipDihedral;
}

double LiftingSurface::getTipDihedralBreak(){
    return tipDihedralBreak_;
}

void LiftingSurface::setTipDihedralBreak( double tipDihedralBreak ){
    tipDihedralBreak_ = tipDihedralBreak;
}

double LiftingSurface::getTaperRatio(){
    return taperRatio_;
}

void LiftingSurface::setTaperRatio( double taperRatio ){
    taperRatio_ = taperRatio;
    tipChord_  = rootChord_*taperRatio_;
}

void LiftingSurface::updateLattice(){
    lattice_.snapToAspectTaperSweep( getAspectRatio(), getTaperRatio(), getSweep());
    lattice_.flipTip( tipDihedralBreak_, tipDihedral_ );
    lattice_.scale( rootChord_ );
    lattice_.rotate( Vec3D(0.0, 0.0, 0.0), Vec3D(1.0, 0.0, 0.0), pitch_);    
}

