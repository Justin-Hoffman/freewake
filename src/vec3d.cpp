#include <cmath>

#include "stdio.h"
#include "vec3d.h"


Vec3D::Vec3D() : x(0.0), y(0.0), z(0.0) {
};

Vec3D::Vec3D( double argx, double argy, double argz) : x(argx), y(argy), z(argz) {
};

Vec3D::Vec3D( const Vec3D &v ) : x(v.x), y(v.y), z(v.z) {
};

Vec3D::~Vec3D(){};

double Vec3D::dot( const Vec3D &v) {
    Vec3D vOut = Vec3D();
    
    vOut.x = this->x*v.x;
    vOut.y = this->y*v.y; 
    vOut.z = this->z*v.z;

    return vOut.x+vOut.y+vOut.z;
};

Vec3D Vec3D::cross( const Vec3D &v) {
    Vec3D vOut = Vec3D();
    
    vOut.x = this->y*v.z - this->z*v.y;
    vOut.y = this->z*v.x - this->x*v.z;
    vOut.z = this->x*v.y - this->y*v.x;

    return vOut;
};

double Vec3D::magnitude(){
    return sqrt( this->x*this->x + this->y*this->y + this->z*this->z);
};

Vec3D Vec3D::norm() {
    Vec3D vOut = Vec3D(*this);
    double m = vOut.magnitude();
    
    return (vOut/m);
};

Vec3D Vec3D::rotate( Vec3D point, Vec3D dir, double theta ){
    Vec3D vOut = Vec3D();
    Vec3D vRelP = *this - point;
    vRelP = vRelP*cos(theta) + dir.cross(vRelP)*sin(theta) + dir * ( dir.dot(vRelP) ) * ( 1.0 - cos(theta) );
    return (point + vRelP);
}

Vec3D Vec3D::operator*( double d ){
    Vec3D vOut = Vec3D(*this);
    vOut.x *= d;
    vOut.y *= d;
    vOut.z *= d;

    return vOut;
};

Vec3D Vec3D::operator/( double d ){
    Vec3D vOut = Vec3D(*this);
    vOut.x /= d;
    vOut.y /= d;
    vOut.z /= d;

    return vOut;
};

Vec3D operator*( double d, Vec3D v ){
    Vec3D vOut = Vec3D(v);
    vOut.x *= d;
    vOut.y *= d;
    vOut.z *= d;
    return vOut;
};

Vec3D Vec3D::operator+( double d ){
    Vec3D vOut = Vec3D(*this);
    vOut.x += d;
    vOut.y += d;
    vOut.z += d;
    return vOut;
};

Vec3D Vec3D::operator+=( double d ){
    this->x += d;
    this->y += d;
    this->z += d;
    return *this;
};

Vec3D operator+( double d, Vec3D v ){
    Vec3D vOut = Vec3D(v);
    vOut.x += d;
    vOut.y += d;
    vOut.z += d;
    return vOut;
};

Vec3D Vec3D::operator+( Vec3D v ){
    Vec3D vOut = Vec3D(*this);
    vOut.x += v.x;
    vOut.y += v.y;
    vOut.z += v.z;
    return vOut;
};

Vec3D Vec3D::operator+=( Vec3D v ){
    this->x += v.x;
    this->y += v.y;
    this->z += v.z;
    return *this;
};


Vec3D Vec3D::operator-( double d ){
    Vec3D vOut = Vec3D(*this);
    vOut.x -= d;
    vOut.y -= d;
    vOut.z -= d;
    return vOut;
};

Vec3D operator-( double d, Vec3D v ){
    Vec3D vOut = Vec3D();
    vOut.x = d-v.x;
    vOut.y = d-v.y;
    vOut.z = d-v.z;
    return vOut;
};

Vec3D Vec3D::operator-( Vec3D v ){
    Vec3D vOut = Vec3D(*this);
    vOut.x -= v.x;
    vOut.y -= v.y;
    vOut.z -= v.z;
    return vOut;
};

void Vec3D::printState(){
    printf("[%7.7f %7.7f %7.7f]\n", this->x, this->y, this->z);
}

