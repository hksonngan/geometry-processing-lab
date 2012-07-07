//=============================================================================
//
//  CLASS Mat3x3 - IMPLEMENTATION
//
//=============================================================================

//== INCLUDES =================================================================


#include "Mat3x3.hh"
#include <math.h>


//== IMPLEMENTATION ========================================================== 


#define MAX_ITER 100
#define EPS 0.00001


int
Mat3x3::symm_eigenv(Vec3d& eigenvals,
		    Vec3d& eigenvec1, Vec3d& eigenvec2, Vec3d& eigenvec3)
{
  unsigned short int i, j;
  int num_iter=0;
  double theta, t, c, s;
  Mat3x3 V(1.0, 0.0, 0.0,
	   0.0, 1.0, 0.0,
	   0.0, 0.0, 1.0),
         R,
         A=*this;
  
  while (num_iter < MAX_ITER) {
    
    // find largest off-diagonal m
    if ( fabs(A.elem_[0][1]) < fabs(A.elem_[0][2]) ) {
      if ( fabs(A.elem_[0][2]) < fabs(A.elem_[1][2]) ) i = 1, j = 2;
      else i = 0, j = 2;
    }
    else {
      if ( fabs(A.elem_[0][1]) < fabs(A.elem_[1][2]) ) i = 1, j = 2;
      else i = 0, j = 1;
    }
    
    if ( fabs(A.elem_[i][j]) < EPS ) break;
    

    // compute Jacobi-Rotation
    theta = 0.5 * (A.elem_[j][j] - A.elem_[i][i]) / A.elem_[i][j];
    t = 1.0 / ( fabs( theta ) + sqrt( 1.0 + theta*theta ) );
    if (theta < 0.0) t = -t;

    c = 1.0 / sqrt(1 + t*t);
    s = t*c;

    R = Mat3x3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    R.elem_[i][i] = R.elem_[j][j] = c;
    R.elem_[i][j] = s;
    R.elem_[j][i] = -s;

    A = R.transpose() * A * R;
    V *= R;
    
    num_iter++;
  }


  if (num_iter < MAX_ITER) {

    // sort and return
    int sorted[3];
    double d[3]={A.elem_[0][0], A.elem_[1][1], A.elem_[2][2]};

    if (d[0] > d[1]) {
      if (d[1] > d[2]) {
	sorted[0] = 0, sorted[1] = 1, sorted[2] = 2;
      }
      else {
	if (d[0] > d[2]) {
	  sorted[0] = 0, sorted[1] = 2, sorted[2] = 1;
	}
	else {
	  sorted[0] = 2, sorted[1] = 0, sorted[2] = 1;
	}
      }
    }
    else {
      if (d[0] > d[2]) {
	sorted[0] = 1, sorted[1] = 0, sorted[2] = 2;
      }
      else {
	if (d[1] > d[2]) {
	  sorted[0] = 1, sorted[1] = 2, sorted[2] = 0;
	}
	else {
	  sorted[0] = 2, sorted[1] = 1, sorted[2] = 0;
	}
      }
    }
            	
    eigenvals = Vec3d( d[sorted[0]],
		       d[sorted[1]],
		       d[sorted[2]] );
    eigenvec1 = Vec3d( V.elem_[0][sorted[0]],
		       V.elem_[1][sorted[0]],
		       V.elem_[2][sorted[0]] );
    eigenvec2 = Vec3d( V.elem_[0][sorted[1]],
		       V.elem_[1][sorted[1]],
		       V.elem_[2][sorted[1]] );

    eigenvec1.normalize();
    eigenvec2.normalize();

    eigenvec3 = eigenvec1 % eigenvec2;
    eigenvec3.normalize();
    
    return(1);
  }
  else return(0);
}

#undef MAX_ITER
#undef EPS


//=============================================================================
