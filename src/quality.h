/*
 * @Author: Kejie Fu
 * @Date: 2023-03-28 03:09:37
 * @LastEditTime: 2023-04-01 00:19:55
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /Mesher3DForSJTU/src/quality.h
 */
#pragma once

template<typename T>
double calculateTetrahedronQualityWith4Points_ISO(const T &a, const T &b, const T &c, const T &d){
    double       abx,aby,abz,acx,acy,acz,adx,ady,adz,bcx,bcy,bcz,bdx,bdy,bdz;
    double       cdx,cdy,cdz;
    double       vol,v1,v2,v3,rap;

    /* volume: (ac^ad).ab/6. Note that here we compute 6*vol. */
    abx = b[0] - a[0];
    aby = b[1] - a[1];
    abz = b[2] - a[2];
    rap = abx*abx + aby*aby + abz*abz;

    acx = c[0] - a[0];
    acy = c[1] - a[1];
    acz = c[2] - a[2];
    rap += acx*acx + acy*acy + acz*acz;

    adx = d[0] - a[0];
    ady = d[1] - a[1];
    adz = d[2] - a[2];
    rap += adx*adx + ady*ady + adz*adz;

    v1  = acy*adz - acz*ady;
    v2  = acz*adx - acx*adz;
    v3  = acx*ady - acy*adx;
    vol = abx * v1 + aby * v2 + abz * v3;
    if ( vol < 1.0e-200)  return vol;

    bcx = c[0] - b[0];
    bcy = c[1] - b[1];
    bcz = c[2] - b[2];
    rap += bcx*bcx + bcy*bcy + bcz*bcz;

    bdx = d[0] - b[0];
    bdy = d[1] - b[1];
    bdz = d[2] - b[2];
    rap += bdx*bdx + bdy*bdy + bdz*bdz;

    cdx = d[0] - c[0];
    cdy = d[1] - c[1];
    cdz = d[2] - c[2];
    rap += cdx*cdx + cdy*cdy + cdz*cdz;
    if ( rap < 1.0e-200 )  return 0.0;

    /* quality = 6*vol / len^3/2. Q = 1/(12 sqrt(3)) for the regular tetra of length 1. */
    rap = rap * sqrt(rap);        
    


    return vol / rap;


}


template<typename T>
double calculateTetrahedronScaleQualityWith4Points_ISO(const T &a, const T &b, const T &c, const T &d){
    double       abx,aby,abz,acx,acy,acz,adx,ady,adz,bcx,bcy,bcz,bdx,bdy,bdz;
    double       cdx,cdy,cdz;
    double       vol,v1,v2,v3,rap;

    /* volume: (ac^ad).ab/6. Note that here we compute 6*vol. */
    abx = b[0] - a[0];
    aby = b[1] - a[1];
    abz = b[2] - a[2];
    rap = abx*abx + aby*aby + abz*abz;

    acx = c[0] - a[0];
    acy = c[1] - a[1];
    acz = c[2] - a[2];
    rap += acx*acx + acy*acy + acz*acz;

    adx = d[0] - a[0];
    ady = d[1] - a[1];
    adz = d[2] - a[2];
    rap += adx*adx + ady*ady + adz*adz;

    v1  = acy*adz - acz*ady;
    v2  = acz*adx - acx*adz;
    v3  = acx*ady - acy*adx;
    vol = abx * v1 + aby * v2 + abz * v3;
    if ( vol < 1.0e-200)  return vol;

    bcx = c[0] - b[0];
    bcy = c[1] - b[1];
    bcz = c[2] - b[2];
    rap += bcx*bcx + bcy*bcy + bcz*bcz;

    bdx = d[0] - b[0];
    bdy = d[1] - b[1];
    bdz = d[2] - b[2];
    rap += bdx*bdx + bdy*bdy + bdz*bdz;

    cdx = d[0] - c[0];
    cdy = d[1] - c[1];
    cdz = d[2] - c[2];
    rap += cdx*cdx + cdy*cdy + cdz*cdz;
    if ( rap < 1.0e-200 )  return 0.0;

    /* quality = 36*vol^{2/3} / (3^{1/3}*rap). Q = 1 for the regular tetra of length 1. */
    //rap = rap * sqrt(rap);        
    
     
    return 36*pow(vol/6.0,2.0/3.0) / (pow(3, 1.0/3.0)* rap);


}