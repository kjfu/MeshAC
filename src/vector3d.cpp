#include "vector3d.h"

double distance(const Vector3D &n0, const Vector3D &n1){
    double rst=0;
    for(int i=0; i<3; i++){
        double d = n0[i] - n1[i];
        rst += d*d;
    }

    return sqrt(rst);
}

// double squareDistance(const Vector3D &v0, const Vector3D &v1){
//     double rst=0;
//     for(int i=0; i<3; i++){
//         double d = v0[i] - v1[i];
//         rst += d*d;
//     }

//     return sqrt(rst);
// }