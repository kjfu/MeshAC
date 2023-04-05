#include "sphere.h"
Sphere::Sphere(Vector3D &v0, Vector3D &v1, Vector3D &v2, Vector3D &v3){
    initialize(v0, v1, v2, v3);
}



void Sphere::initialize(Vector3D &v0, Vector3D &v1, Vector3D &v2, Vector3D &v3){
    double cent[3];
    double r;
    circumsphere(v0.data(), v1.data(), v2.data(), v3.data(), cent, &r);
    center.initialize(cent);
    radius=r;
}

bool Sphere::contain(const Vector3D &pos){
    double d = distance(center, pos);
    bool rst = false;
    if (d<=radius){
        rst = true;
    }
    return rst;
}

bool circumsphere(REAL* pa, REAL* pb, REAL* pc, REAL* pd, REAL* cent, REAL* radius){
    REAL A[4][4], rhs[4], D;
    int indx[4];

    // Compute the coefficient matrix A (3x3).
    A[0][0] = pb[0] - pa[0];
    A[0][1] = pb[1] - pa[1];
    A[0][2] = pb[2] - pa[2];
    A[1][0] = pc[0] - pa[0];
    A[1][1] = pc[1] - pa[1];
    A[1][2] = pc[2] - pa[2];
    if (pd != NULL) {
        A[2][0] = pd[0] - pa[0];
        A[2][1] = pd[1] - pa[1]; 
        A[2][2] = pd[2] - pa[2];
    } else {
        Cross(A[0], A[1], A[2]);
    }

    // Compute the right hand side vector b (3x1).
    rhs[0] = 0.5 * Dot(A[0], A[0]);
    rhs[1] = 0.5 * Dot(A[1], A[1]);
    if (pd != NULL) {
        rhs[2] = 0.5 * Dot(A[2], A[2]);
    } else {
        rhs[2] = 0.0;
    }

    // Solve the 3 by 3 equations use LU decomposition with partial pivoting
    //   and backward and forward substitute..
    if (!LU_decmp(A, 3, indx, &D, 0)) {
        if (radius != (REAL *) NULL) *radius = 0.0;
        return false;
    }    
    LU_solve(A, 3, indx, rhs, 0);
    if (cent != (REAL *) NULL) {
        cent[0] = pa[0] + rhs[0];
        cent[1] = pa[1] + rhs[1];
        cent[2] = pa[2] + rhs[2];
    }
    if (radius != (REAL *) NULL) {
        *radius = sqrt(rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2]);
    }
    return true;
}

// lu_decmp()    Compute the LU decomposition of a matrix.                   //
//                                                                           //
// Compute the LU decomposition of a (non-singular) square matrix A using    //
// partial pivoting and implicit row exchanges.  The result is:              //
//     A = P * L * U,                                                        //
// where P is a permutation matrix, L is unit lower triangular, and U is     //
// upper triangular.  The factored form of A is used in combination with     //
// 'lu_solve()' to solve linear equations: Ax = b, or invert a matrix.       //
//                                                                           //
// The inputs are a square matrix 'lu[N..n+N-1][N..n+N-1]', it's size is 'n'.//
// On output, 'lu' is replaced by the LU decomposition of a rowwise permuta- //
// tion of itself, 'ps[N..n+N-1]' is an output vector that records the row   //
// permutation effected by the partial pivoting, effectively,  'ps' array    //
// tells the user what the permutation matrix P is; 'd' is output as +1/-1   //
// depending on whether the number of row interchanges was even or odd,      //
// respectively.                                                             //
//                                                                           //
// Return true if the LU decomposition is successfully computed, otherwise,  //
// return false in case that A is a singular matrix.                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

bool LU_decmp(REAL lu[4][4], int n, int* ps, REAL* d, int N)
{
  REAL scales[4];
  REAL pivot, biggest, mult, tempf;
  int pivotindex = 0;
  int i, j, k;

  *d = 1.0;                                      // No row interchanges yet.

  for (i = N; i < n + N; i++) {                             // For each row.
    // Find the largest element in each row for row equilibration
    biggest = 0.0;
    for (j = N; j < n + N; j++)
      if (biggest < (tempf = fabs(lu[i][j])))
        biggest  = tempf;
    if (biggest != 0.0)
      scales[i] = 1.0 / biggest;
    else {
      scales[i] = 0.0;
      return false;                            // Zero row: singular matrix.
    }
    ps[i] = i;                                 // Initialize pivot sequence.
  }

  for (k = N; k < n + N - 1; k++) {                      // For each column.
    // Find the largest element in each column to pivot around.
    biggest = 0.0;
    for (i = k; i < n + N; i++) {
      if (biggest < (tempf = fabs(lu[ps[i]][k]) * scales[ps[i]])) {
        biggest = tempf;
        pivotindex = i;
      }
    }
    if (biggest == 0.0) {
      return false;                         // Zero column: singular matrix.
    }
    if (pivotindex != k) {                         // Update pivot sequence.
      j = ps[k];
      ps[k] = ps[pivotindex];
      ps[pivotindex] = j;
      *d = -(*d);                          // ...and change the parity of d.
    }

    // Pivot, eliminating an extra variable  each time
    pivot = lu[ps[k]][k];
    for (i = k + 1; i < n + N; i++) {
      lu[ps[i]][k] = mult = lu[ps[i]][k] / pivot;
      if (mult != 0.0) {
        for (j = k + 1; j < n + N; j++)
          lu[ps[i]][j] -= mult * lu[ps[k]][j];
      }
    }
  }

  // (lu[ps[n + N - 1]][n + N - 1] == 0.0) ==> A is singular.
  return lu[ps[n + N - 1]][n + N - 1] != 0.0;
}


// lu_solve()    Solves the linear equation:  Ax = b,  after the matrix A    //
//               has been decomposed into the lower and upper triangular     //
//               matrices L and U, where A = LU.                             //
//                                                                           //
// 'lu[N..n+N-1][N..n+N-1]' is input, not as the matrix 'A' but rather as    //
// its LU decomposition, computed by the routine 'lu_decmp'; 'ps[N..n+N-1]'  //
// is input as the permutation vector returned by 'lu_decmp';  'b[N..n+N-1]' //
// is input as the right-hand side vector, and returns with the solution     //
// vector. 'lu', 'n', and 'ps' are not modified by this routine and can be   //
// left in place for successive calls with different right-hand sides 'b'.   //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void LU_solve(REAL lu[4][4], int n, int* ps, REAL* b, int N)
{
    int i, j;
    REAL X[4], dot;

    for (i = N; i < n + N; i++) X[i] = 0.0;

    // Vector reduction using U triangular matrix.
    for (i = N; i < n + N; i++) {
        dot = 0.0;
        for (j = N; j < i + N; j++)
        dot += lu[ps[i]][j] * X[j];
        X[i] = b[ps[i]] - dot;
    }

    // Back substitution, in L triangular matrix.
    for (i = n + N - 1; i >= N; i--) {
        dot = 0.0;
        for (j = i + 1; j < n + N; j++)
        dot += lu[ps[i]][j] * X[j];
        X[i] = (X[i] - dot) / lu[ps[i]][i];
    }

    for (i = N; i < n + N; i++) b[i] = X[i];
}
