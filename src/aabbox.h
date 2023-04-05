#pragma once
#include <array>
#include "vector3d.h"

class AABBox
{
public:

    std::array<double, 3> maximum;
    std::array<double, 3> minimum;
    AABBox(){
        maximum.fill(std::numeric_limits<double>::min());
        minimum.fill(std::numeric_limits<double>::max());
    }

    AABBox(const std::array<double,3> &min, const std::array<double, 3> &max):maximum(max), minimum(min){
    }

    void initialize(const std::array<double,3> &min, const std::array<double, 3> &max){
        minimum = min;
        maximum = max;
    }

    void insert(const std::array<double, 3> &n){
        for(int i=0; i<3; i++){
            if (n[i]>maximum[i]){
                maximum[i] = n[i];
            }
            
            if(n[i]<minimum[i]){
                minimum[i] = n[i];
            }
        }
    }

    void insert(Vector3D &n){
        for(int i=0; i<3; i++){
            if (n[i]>maximum[i]){
                maximum[i] = n[i];
            }
            
            if(n[i]<minimum[i]){
                minimum[i] = n[i];
            }
        }
    }

    bool isValid(){
        bool rst = true;
        for(int i=0; i<3; i++){
            if(minimum[i]>maximum[i]){
                rst = false;
                break;
            }
        }

        return rst;
    }


    void reset(){
        maximum.fill(std::numeric_limits<double>::min());
        minimum.fill(std::numeric_limits<double>::max());        
    }

    bool contain(Vector3D pos, double epsilon = std::numeric_limits<double>::epsilon()){
        for(int i=0; i<3; i++){
            if(pos[i]>(maximum[i]+epsilon) || pos[i]<(minimum[i]-epsilon)){
                return false;
            }
        }
        return true;
    }


    bool intersects(AABBox &another,  double tolerance = std::numeric_limits<double>::epsilon()){
       
        return (another.maximum[0] >= this->minimum[0] - tolerance 
        && another.minimum[0] <= this->maximum[0] + tolerance
        && another.maximum[1] >= this->minimum[1] - tolerance 
        && another.minimum[1] <= this->maximum[1] + tolerance
        && another.maximum[2] >= this->minimum[2] - tolerance 
        && another.minimum[2] <= this->maximum[2] + tolerance);
    }
};
