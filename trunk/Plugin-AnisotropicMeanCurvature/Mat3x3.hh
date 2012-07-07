//=============================================================================
//
//  CLASS Mat3x3
//
//=============================================================================

#ifndef MAT3X3_HH
#define MAT3X3_HH


//== INCLUDES =================================================================

// #include <CG2/common/VectorT.hh>
// #include <OpenMesh/Core/Geometry/VectorT.hh>
#include <ACG/Math/VectorT.hh>

//== CLASS DEFINITION =========================================================

using namespace ACG;

class Mat3x3
{
public:
   
  
  /// Empty constructor - object remains uninitialized.
  Mat3x3( void ) {};

  /// Constructor 
  Mat3x3( double a00, double a01, double a02, 
          double a10, double a11, double a12, 
          double a20, double a21, double a22)
  {
    elem_[0][0] = a00; elem_[0][1] = a01; elem_[0][2] = a02;
    elem_[1][0] = a10; elem_[1][1] = a11; elem_[1][2] = a12;
    elem_[2][0] = a20; elem_[2][1] = a21; elem_[2][2] = a22;
  };
  
  /// Copy-Constructor
  Mat3x3( const Mat3x3 &p )
  {
    elem_[0][0] = p.elem_[0][0];
    elem_[0][1] = p.elem_[0][1];
    elem_[0][2] = p.elem_[0][2];
    elem_[1][0] = p.elem_[1][0];
    elem_[1][1] = p.elem_[1][1];
    elem_[1][2] = p.elem_[1][2];
    elem_[2][0] = p.elem_[2][0];
    elem_[2][1] = p.elem_[2][1];
    elem_[2][2] = p.elem_[2][2];
  };

  
  /// Destructor
  ~Mat3x3( void ) {};

  
  //
  // Operators
  //
  
  /// +=-Operator
  Mat3x3& operator+=( const Mat3x3& p )
  {
    elem_[0][0] += p.elem_[0][0];
    elem_[0][1] += p.elem_[0][1];
    elem_[0][2] += p.elem_[0][2];
    elem_[1][0] += p.elem_[1][0];
    elem_[1][1] += p.elem_[1][1];
    elem_[1][2] += p.elem_[1][2];
    elem_[2][0] += p.elem_[2][0];
    elem_[2][1] += p.elem_[2][1];
    elem_[2][2] += p.elem_[2][2];
    return *this;
  };

  /// -=-Operator
  Mat3x3& operator-=( const Mat3x3& p )
  {
    elem_[0][0] -= p.elem_[0][0];
    elem_[0][1] -= p.elem_[0][1];
    elem_[0][2] -= p.elem_[0][2];
    elem_[1][0] -= p.elem_[1][0];
    elem_[1][1] -= p.elem_[1][1];
    elem_[1][2] -= p.elem_[1][2];
    elem_[2][0] -= p.elem_[2][0];
    elem_[2][1] -= p.elem_[2][1];
    elem_[2][2] -= p.elem_[2][2];
    return *this;
  };

  /// /=-Operator
  Mat3x3& operator/=( double s ) 
  {
    elem_[0][0] /= s;  elem_[0][1] /= s;  elem_[0][2] /= s;
    elem_[1][0] /= s;  elem_[1][1] /= s;  elem_[1][2] /= s;
    elem_[2][0] /= s;  elem_[2][1] /= s;  elem_[2][2] /= s;
    return *this;
  };

  /// *=-Operator : Matrix * Scalar
  Mat3x3& operator*=( double s ) 
  {
    elem_[0][0] *= s;  elem_[0][1] *= s;  elem_[0][2] *= s;
    elem_[1][0] *= s;  elem_[1][1] *= s;  elem_[1][2] *= s;
    elem_[2][0] *= s;  elem_[2][1] *= s;  elem_[2][2] *= s;
    return *this;
  };

  /// *=-Operator : Matrix * Matrix
  Mat3x3& operator*=( const Mat3x3& p )
  {
    return ( *this = *this * p );
  };

  /// +-Operator 
  Mat3x3  operator+( const Mat3x3& p ) const
  {
    return Mat3x3( elem_[0][0]+p.elem_[0][0],
                   elem_[0][1]+p.elem_[0][1],
                   elem_[0][2]+p.elem_[0][2],
                   elem_[1][0]+p.elem_[1][0],
                   elem_[1][1]+p.elem_[1][1],
                   elem_[1][2]+p.elem_[1][2],
                   elem_[2][0]+p.elem_[2][0],
                   elem_[2][1]+p.elem_[2][1],
                   elem_[2][2]+p.elem_[2][2] );
  };

  /// --Operator
  Mat3x3  operator-( const Mat3x3& p ) const 
  {
    return Mat3x3( elem_[0][0]-p.elem_[0][0],
                   elem_[0][1]-p.elem_[0][1],
                   elem_[0][2]-p.elem_[0][2],
                   elem_[1][0]-p.elem_[1][0],
                   elem_[1][1]-p.elem_[1][1],
                   elem_[1][2]-p.elem_[1][2],
                   elem_[2][0]-p.elem_[2][0],
                   elem_[2][1]-p.elem_[2][1],
                   elem_[2][2]-p.elem_[2][2] );
  };
 
  ///-Operator 
  Mat3x3  operator/( double s ) const 
  {
    return Mat3x3( elem_[0][0]/s, elem_[0][1]/s, elem_[0][2]/s,
                   elem_[1][0]/s, elem_[1][1]/s, elem_[1][2]/s,
                   elem_[2][0]/s, elem_[2][1]/s, elem_[2][2]/s );
  };
  
  /// *-Operator : Matrix * Scalar
  Mat3x3  operator*( double s ) const
  {
    return Mat3x3( elem_[0][0]*s, elem_[0][1]*s, elem_[0][2]*s,
                   elem_[1][0]*s, elem_[1][1]*s, elem_[1][2]*s,
                   elem_[2][0]*s, elem_[2][1]*s, elem_[2][2]*s );
  };

  /// friend *-Operator : Scalar * Matrix
  friend Mat3x3 operator*( double s, Mat3x3& p )
  {
    return Mat3x3( p.elem_[0][0]*s, p.elem_[0][1]*s, p.elem_[0][2]*s,
                   p.elem_[1][0]*s, p.elem_[1][1]*s, p.elem_[1][2]*s,
                   p.elem_[2][0]*s, p.elem_[2][1]*s, p.elem_[2][2]*s );  
  };
  
  /// *-Operator : Matrix * Vector
  Vec3d operator*( const Vec3d& vec ) const 
  {
    return Vec3d(elem_[0][0]*vec[0] + elem_[0][1]*vec[1] + elem_[0][2]*vec[2], 
                 elem_[1][0]*vec[0] + elem_[1][1]*vec[1] + elem_[1][2]*vec[2],
                 elem_[2][0]*vec[0] + elem_[2][1]*vec[1] + elem_[2][2]*vec[2]);
  };

  /// *-Operator : Matrix * Matrix
  Mat3x3  operator*( const Mat3x3& p ) const
  {
    Mat3x3 result;
    result.elem_[0][0] = ( elem_[0][0]*p.elem_[0][0] +
                           elem_[0][1]*p.elem_[1][0] +
                           elem_[0][2]*p.elem_[2][0] );
    result.elem_[0][1] = ( elem_[0][0]*p.elem_[0][1] +
                           elem_[0][1]*p.elem_[1][1] +
                           elem_[0][2]*p.elem_[2][1] );
    result.elem_[0][2] = ( elem_[0][0]*p.elem_[0][2] +
                           elem_[0][1]*p.elem_[1][2] +
                           elem_[0][2]*p.elem_[2][2] );
    result.elem_[1][0] = ( elem_[1][0]*p.elem_[0][0] +
                           elem_[1][1]*p.elem_[1][0] +
                           elem_[1][2]*p.elem_[2][0] );
    result.elem_[1][1] = ( elem_[1][0]*p.elem_[0][1] +
                           elem_[1][1]*p.elem_[1][1] +
                           elem_[1][2]*p.elem_[2][1] );
    result.elem_[1][2] = ( elem_[1][0]*p.elem_[0][2] +
                           elem_[1][1]*p.elem_[1][2] +
                           elem_[1][2]*p.elem_[2][2] );
    result.elem_[2][0] = ( elem_[2][0]*p.elem_[0][0] +
                           elem_[2][1]*p.elem_[1][0] +
                           elem_[2][2]*p.elem_[2][0] );
    result.elem_[2][1] = ( elem_[2][0]*p.elem_[0][1] +
                           elem_[2][1]*p.elem_[1][1] +
                           elem_[2][2]*p.elem_[2][1] );
    result.elem_[2][2] = ( elem_[2][0]*p.elem_[0][2] +
                           elem_[2][1]*p.elem_[1][2] +
                           elem_[2][2]*p.elem_[2][2] );
    return result;
  };


  /// read access for matrix elements
  double operator() (int _row, int _col) const { return elem_[_row][_col]; };

  /// write access for matrix elements
  double& operator() (int _row, int _col) { return elem_[_row][_col]; };


  //
  // Methods
  //

   /// determinant of 3x3 Matrix
  double  det()
  {
    return ((elem_[0][1]*elem_[1][2] - elem_[0][2]*elem_[1][1]) * elem_[2][0]
          + (elem_[0][2]*elem_[1][0] - elem_[0][0]*elem_[1][2]) * elem_[2][1]
          + (elem_[0][0]*elem_[1][1] - elem_[0][1]*elem_[1][0]) * elem_[2][2]);
  };


  /// Transposed Matrix
  Mat3x3 transpose()
  {
    return( Mat3x3( elem_[0][0], elem_[1][0], elem_[2][0],
                    elem_[0][1], elem_[1][1], elem_[2][1],
                    elem_[0][2], elem_[1][2], elem_[2][2] ) );
  }


  /// Eigenvalues and Eigenvectors
  int symm_eigenv(Vec3d& eigenvals, 
                  Vec3d& eigenvec1, 
                  Vec3d& eigenvec2, 
                  Vec3d& eigenvec3);
  

private:

  /// matrix ments
  double elem_[3][3];
};


//=============================================================================
//}//namespace
//=============================================================================
#endif // MAT3X3_HH defined

