#ifndef VEC3D_H
#define VEC3D_H

class Vec3D {
    public:
        double x;
        double y; 
        double z; 
    
        Vec3D();
        Vec3D( double argx, double argy, double argz);
        Vec3D( const Vec3D &v ); 
        ~Vec3D();
        
        double dot( const Vec3D &v );
        Vec3D cross( const Vec3D &v );
        double magnitude();
        Vec3D norm();
        Vec3D rotate( Vec3D point, Vec3D dir, double theta );

        Vec3D operator* ( double );
        Vec3D operator/ ( double );
        Vec3D operator+ ( double );
        Vec3D operator- ( double );
        Vec3D operator+ ( Vec3D v );
        Vec3D operator- ( Vec3D v );
        
        Vec3D operator+=( double );
        Vec3D operator+=( Vec3D v );
            
        void printState();
    private:
 
};

Vec3D operator*(double, Vec3D);
Vec3D operator+(double, Vec3D);
Vec3D operator-(double, Vec3D);

#endif
