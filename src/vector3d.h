#pragma once
#include<array>
#include<cmath>
#include <vector>
class Vector3D;
double distance(const Vector3D &n0, const Vector3D &n1);


class Vector3D{    
public:      
    std::array<double, 3> xyz;

    explicit Vector3D(double value=0){
        xyz.fill(value);
    }

    Vector3D(double *vec){
        xyz[0] = vec[0];
        xyz[1] = vec[1];
        xyz[2] = vec[2];
    }

    Vector3D(double x, double y, double z){
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
    }

    Vector3D(const Vector3D &v){
        this->xyz = v.xyz;
    }

    Vector3D &operator=(const Vector3D &v){
        this->xyz = v.xyz;
        return *this;
    }

    double &operator[](const int &i){
        return this->xyz[i];
    }
    double operator[](const int &i) const{
        return this->xyz[i];
    }

    Vector3D operator + (const Vector3D &v) const{
        Vector3D rst(*this);
        for(int i=0; i<3; i++){
            rst[i] += v[i];
        }
        return rst;
    }
    Vector3D operator + (const double &value) const{
        Vector3D rst(*this);
        for(int i=0; i<3; i++){
            rst[i] += value;
        }
        return rst;
    }
    void  operator +=(const Vector3D &v){

        for(int i=0; i<3; i++){
            xyz[i] += v[i];
        }
    }
    void  operator +=(const double &value){
        for(int i=0; i<3; i++){
            xyz[i] += value;
        }
    }

    void  operator -=(const Vector3D &v){
        for(int i=0; i<3; i++){
            xyz[i] -= v[i];
        }
    }

    void  operator -=(const double &value){
        for(int i=0; i<3; i++){
            xyz[i] -= value;
        }
    }

    void  operator /=(const double value){
        for(int i=0; i<3; i++){
            xyz[i] /= value;
        }
    }

    void  operator *=(const double value){
        for(int i=0; i<3; i++){
            xyz[i] *= value;
        }
    }

    Vector3D operator - (const Vector3D &v) const{
        Vector3D rst(*this);
        for(int i=0; i<3; i++){
            rst[i] -= v[i];
        }
        return rst;
    }
    Vector3D operator - (const double &value) const{
        Vector3D rst(*this);
        for(int i=0; i<3; i++){
            rst[i] -= value;
        }
        return rst;
    }

    Vector3D operator / (const double value) const{
        Vector3D rst(*this);
        for(int i=0; i<3; i++){
            rst[i] /= value;
        }
        return rst;
    }
    friend Vector3D operator* (const double &k, const Vector3D &v){
        Vector3D rst(v);
        for(int i=0; i<3; i++){
            rst[i] *= k;
        }
        return rst;
    }
    friend Vector3D operator* (const Vector3D &v, const double &k){
        Vector3D rst(v);
        for(int i=0; i<3; i++){
            rst[i] *= k;
        }
        return rst;
    }

    double dot(const Vector3D &v) const{
        double rst=0;
        for(int i=0; i<3; i++){
            rst += xyz[i] * v[i];
        }
        return rst;
    }

    double *data(){
        return xyz.data();
    }

    std::array<double, 3> XYZ(){
        return xyz;
    }

    double norm(){
        double rst=0;
        for(auto v: xyz){
            rst += v*v;
        }
        return sqrt(rst);        
    }

    void normalize(){
        double len = this->norm();
        for(int i=0; i<3; i++){
            xyz[i] /= len;
        }
    }

    Vector3D normalized(){
        Vector3D rst(*this);
        double len = this->norm();
        for(int i=0; i<3; i++){
            rst.xyz[i] /= len;
        }
        return rst;
    }

    void initialize(double *vec){
        for(int i=0; i<3; i++){
            xyz[i] = vec[i];
        }
    }

    std::vector<double> toSTDVector(){
        std::vector<double> rst(3);
        rst[0] = xyz[0];
        rst[1] = xyz[1];
        rst[2] = xyz[2];
        return rst;
    }
};





