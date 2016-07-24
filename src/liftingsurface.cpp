#include "liftingsurface.h"

LiftingSurface::LiftingSurface() : 
                freeWake_(false), freeTipVortex_(false), span_( 1.0 ), sweep_( 0.0 ), pitch_(0.0), rc_( 1E-6 ), rootChord_( 1.0 ), tipChord_(1.0), tipDihedral_(0.0), tipDihedralBreak_(0.0), taperRatio_(1.0), 
                nSpan_( 10 ), nChord_( 2 ), nWake_( 2 ), nFilament_( 2 ),  spanwiseForce_( 10, Vec3D() ),
                horseshoeLattice_( nSpan_, nChord_), vortexLattice_( nSpan_+1, nWake_ ), tipFilament_( 2, nFilament_) {
    horseshoeLattice_.snapToUnit(); 
}

LiftingSurface::LiftingSurface( const LiftingSurface &l) : 
                freeWake_( l.freeWake_ ), freeTipVortex_( l.freeTipVortex_ ), span_( l.span_ ), sweep_( l.sweep_ ), pitch_( l.pitch_ ), rc_( l.rc_), rootChord_( l.rootChord_ ), tipChord_( l.tipChord_ ), tipDihedral_( l.tipDihedral_ ), 
                tipDihedralBreak_( l.tipDihedralBreak_ ), taperRatio_( l.taperRatio_ ), 
                nSpan_( l.nSpan_ ), nChord_( l.nChord_ ), nWake_( l.nWake_ ), nFilament_( l.nFilament_ ), spanwiseForce_( l.spanwiseForce_ ),
                horseshoeLattice_( l.horseshoeLattice_ ), vortexLattice_( l.vortexLattice_ ), tipFilament_( l.tipFilament_ ){
}

LiftingSurface::LiftingSurface( int nSpan, int nChord ) :
                freeWake_( false ), freeTipVortex_( false ), span_( 1.0 ), sweep_( 0.0 ), pitch_(0.0), rc_( 1E-6 ), rootChord_( 1.0 ), tipChord_(1.0), tipDihedral_(0.0), tipDihedralBreak_(0.0), taperRatio_(1.0), 
                nSpan_( nSpan ), nChord_( nChord ), nWake_( 2 ), nFilament_( 2 ), spanwiseForce_( nSpan, Vec3D() ),
                horseshoeLattice_( nSpan_, nChord_ ), vortexLattice_( nSpan_+1, nWake_ ), tipFilament_( 2, 2 ) {
    horseshoeLattice_.snapToUnit();
}        

LiftingSurface::LiftingSurface( int nSpan, int nChord, int nWake, int nFilament ) :
                freeWake_( false ), freeTipVortex_( false ), span_( 1.0 ), sweep_( 0.0 ), pitch_(0.0), rc_( 1E-6), rootChord_( 1.0 ), tipChord_(1.0), tipDihedral_(0.0), tipDihedralBreak_(0.0), taperRatio_(1.0), 
                nSpan_( nSpan ), nChord_( nChord ), nWake_( nWake ), nFilament_( nFilament ), spanwiseForce_( nSpan, Vec3D() ),
                horseshoeLattice_( nSpan_, nChord_ ), vortexLattice_( nSpan_+1, nWake_ ), tipFilament_( 2, nFilament_ ) {
    horseshoeLattice_.snapToUnit();
}        

LiftingSurface::~LiftingSurface(){
}

bool LiftingSurface::freeWake(){
    return freeWake_;
}

bool LiftingSurface::freeTipVortex(){
    return freeTipVortex_;
}

void LiftingSurface::setFreeWake( bool freeWake ){
    freeWake_ = freeWake;
}

void LiftingSurface::setFreeTipVortex( bool freeTipVortex ){
    freeTipVortex_ = freeTipVortex;
}

HorseshoeLattice& LiftingSurface::getHorseshoeLattice(){
    return horseshoeLattice_;
}

VortexLattice& LiftingSurface::getVortexLattice(){
    return vortexLattice_;
}

TipFilament& LiftingSurface::getTipFilament(){
    return tipFilament_;
}

std::vector<Vec3D>& LiftingSurface::spanwiseForce(){
    return spanwiseForce_;
}

int LiftingSurface::nSpan(){
    return nSpan_;
}

int LiftingSurface::nChord(){
    return nChord_;
}

int LiftingSurface::nWake(){
    return nWake_;
}

int LiftingSurface::nFilament(){
    return nFilament_;
}

double LiftingSurface::getAspectRatio(){
    return 2.0*span_ / ( rootChord_ + tipChord_ );
}

void LiftingSurface::setAspectRatio( double AR ){
    span_ =  AR*( rootChord_ + tipChord_ ) / 2.0;
}

void LiftingSurface::setCoreRadius( double rc ){
    rc_ = rc;
    horseshoeLattice_.setRc( rc*10.0 );
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

void LiftingSurface::setVortexLattice( VortexLattice &v ){
    vortexLattice_ = VortexLattice( v );
}

void LiftingSurface::updateLattice(){
    horseshoeLattice_.snapToAspectTaperSweep( getAspectRatio(), getTaperRatio(), getSweep());
    horseshoeLattice_.flipTip( tipDihedralBreak_, tipDihedral_ );
    horseshoeLattice_.scale( rootChord_ );
    horseshoeLattice_.rotate( Vec3D(0.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0.0), pitch_);
    vortexLattice_.fixToTrailingEdge( horseshoeLattice_ ); 
}

Vec3D LiftingSurface::calcInducedVelocity( Vec3D p ){
    Vec3D vInduced = horseshoeLattice_.calcInducedVelocity( p );
    if ( freeWake_ ){
        vInduced += vortexLattice_.calcInducedVelocity( p );
    }
    if ( freeTipVortex_ ){
        vInduced += tipFilament_.calcInducedVelocity( p );
    }
    return vInduced;
}
