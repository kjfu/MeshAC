#pragma once
#include "vector3d.h"


class Sphere{
public:
    Vector3D center;
    double radius;
    Sphere(){

    }

    Sphere(Vector3D &v0, Vector3D &v1, Vector3D &v2, Vector3D &v3);
    void initialize(Vector3D &v0, Vector3D &v1, Vector3D &v2, Vector3D &v3);
    bool contain(const Vector3D &pos);


};















#define REAL double
bool circumsphere(REAL* pa, REAL* pb, REAL* pc, REAL* pd, REAL* cent, REAL* radius);
void LU_solve(REAL lu[4][4], int n, int* ps, REAL* b, int N);
bool LU_decmp(REAL lu[4][4], int n, int* ps, REAL* d, int N);

// dot() returns the dot product: v1 dot v2.
inline REAL Dot(REAL* v1, REAL* v2) 
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

// cross() computes the cross product: n = v1 cross v2.
inline void Cross(REAL* v1, REAL* v2, REAL* n) 
{
  n[0] =   v1[1] * v2[2] - v2[1] * v1[2];
  n[1] = -(v1[0] * v2[2] - v2[0] * v1[2]);
  n[2] =   v1[0] * v2[1] - v2[0] * v1[1];
}

// distance() computes the Euclidean distance between two points.
inline REAL Distance(REAL* p1, REAL* p2)
{
  return sqrt((p2[0] - p1[0]) * (p2[0] - p1[0]) +
              (p2[1] - p1[1]) * (p2[1] - p1[1]) +
              (p2[2] - p1[2]) * (p2[2] - p1[2]));
}
