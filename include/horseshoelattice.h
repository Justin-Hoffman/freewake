#ifndef HORSESHOELATTICE_H
#define HORSESHOELATTICE_H

#include <vector>
#include <utility>

#include "vec3d.h"
#include "vortexcontainer.h"

enum PointSpacing { Linear, Cosine, HalfCosine };

class HorseshoeLattice : public VortexContainer{
    public: 
        HorseshoeLattice();
        HorseshoeLattice( int ni, int nj );
        HorseshoeLattice( const HorseshoeLattice &vl );
        ~HorseshoeLattice();

        int ni();
        int nj();
       
        bool hasTrailers();
        Vec3D trailerVec();
       
        PointSpacing chordwiseSpacing();
        PointSpacing spanwiseSpacing();
        std::vector< std::vector<Vec3D> >& getEndPoints();
        std::vector< std::vector<Vec3D> >& getControlPoints();
        std::vector< std::vector<Vec3D> >& getControlPointNormals();
        std::vector< std::vector<double> >& getGamma();
        std::vector< double >& getNetGamma();
        
        virtual Vec3D calcInfluenceCoefficient( Vec3D p, int n);
        virtual Vec3D calcInducedVelocity( Vec3D );
        
        void calcControlPointNormals(); 
        
        std::pair<int, int> ijFromN( int n );
        int ijToN( int i, int j);
        int maxN();
        double dSpan(int i, int j);
        double dChord(int i);
        double getRc( );
        Vec3D gammaVector(int i, int j);
        Vec3D gammaCenterPoint(int i, int j);
       
        void centerControlPoints(); 
        void chordwiseSpacing( PointSpacing );
        void spanwiseSpacing( PointSpacing );
        void flipTip( double dihedralBreak, double dihedral);  
        void setHasTrailers( bool b );
        void setTrailerVec( Vec3D v );
        void snapToUnit();
        void snapToAspectTaper( double ar, double taper ); 
        void snapToAspectTaperSweep( double ar, double taper, double Sweep ); 
        void printState();

        void rotate( Vec3D point, Vec3D dir, double theta); 
        void scale( double scale ); 
        void translate( Vec3D dir ); 
        void setRc( double );
        
    private:
        double rc_;
        int ni_;
        int nj_;
        PointSpacing chordwiseSpacing_;
        PointSpacing spanwiseSpacing_;
        bool hasTrailers_;
        Vec3D trailerVec_;


        std::vector< std::vector<Vec3D> > endPoints;
        std::vector< std::vector<Vec3D> > controlPoints;
        std::vector< std::vector<Vec3D> > controlPointNormals;
        std::vector< std::vector<double> > gamma;
        
        std::vector<double> netGamma_;

};

#endif
