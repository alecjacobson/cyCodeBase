#ifndef IGL_BLUE_NOISE_H
#define IGL_BLUE_NOISE_H
#include <Eigen/Core>
namespace igl
{
  // Compute n blue noise samples on a given mesh
  //
  // Inputs:
  //   n  number of samples
  //   V  #V by dim list of mesh vertex points
  //   F  #F by 3 list of traingle indices into V
  // Outputs:
  //   B  n by 3 list of barycentric coordinates
  //   I  n by 1 list of indices into F
  //   X  n by dim list of 
  //
  inline void blue_noise(
    const int n,
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & B,
    Eigen::VectorXi & I,
    Eigen::MatrixXd & X);
}

// Implementation
#include "cyPoint.h"
#include "cySampleElim.h"
#include <igl/doublearea.h>
#include <igl/random_points_on_mesh.h>

namespace cy
{
  struct Point3f_id : public cy::Point3f
  {
    int i;
    Point3f_id( float v): cy::Point3f(v),i(-1){}
    Point3f_id(){}
    Point3f_id(const Point3f_id & that): cy::Point3f(that), i(that.i) {}
    Point3f_id(const Point3f & that): cy::Point3f(that), i(-1) {}
    Point3f_id operator-(Point3f_id const &that) const 
    {
      return cy::Point3f::operator-(that);
    };
  };
}

inline void igl::blue_noise(
  const int n,
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & B,
  Eigen::VectorXi & I,
  Eigen::MatrixXd & X)
{
  Eigen::MatrixXd X5;
  Eigen::MatrixXd B5;
  Eigen::VectorXi I5;
  igl::random_points_on_mesh(5*n,V,F,B5,I5,X5);
  std::vector< cy::Point3f_id > inputPoints(X5.rows());
  for ( size_t i=0; i<inputPoints.size(); i++ ) {
      inputPoints[i].x = X5(i,0);
      inputPoints[i].y = X5(i,1);
      inputPoints[i].z = X5(i,2);
      inputPoints[i].i = i;
  }
  cy::WeightedSampleElimination< cy::Point3f_id, float, 3, int > wse;
  std::vector< cy::Point3f_id > outputPoints(n);
  Eigen::VectorXd A;
  igl::doublearea(V,F,A);
  float area = A.array().sum()/2;
  float d_max = 2 * wse.GetMaxPoissonDiskRadius( 2, outputPoints.size(), area );
  wse.Eliminate( 
    inputPoints.data(), inputPoints.size(), 
    outputPoints.data(), outputPoints.size(),
    false,
    d_max, 2
    );
  //typedef float FType;
  //typedef cy::Point3f PointType;
  //FType beta  = 0.65;
  //FType gamma = 1.5;
  //typedef int SIZE_TYPE;
  //const auto GetWeightLimitFraction = [beta,gamma]( SIZE_TYPE inputSize, SIZE_TYPE outputSize )->FType
  //{
  //  FType ratio = FType(outputSize) / FType(inputSize);
  //  return ( 1 - std::pow( ratio, gamma ) ) * beta;
  //};
  //FType d_min = d_max * GetWeightLimitFraction( inputPoints.size(), outputPoints.size());
  //FType alpha = 8;
  //Eigen::RowVector3d minV = V.colwise().minCoeff();
  //Eigen::RowVector3d maxV = V.colwise().maxCoeff();
  //Eigen::RowVector3d meanV = V.colwise().mean();
  //const auto weight = 
  //  [d_min, alpha](
  //  PointType const & p0, 
  //  PointType const & p1, 
  //  FType d2, 
  //  FType d_max)
  //{
  //  FType d = cy::Sqrt(d2)*(3.0-2.0*(1.0-p0.Length()));
  //  return std::pow( FType(1) - d/d_max, alpha );
  //};
  //wse.Eliminate( 
  //  inputPoints.data(), inputPoints.size(), 
  //  outputPoints.data(), outputPoints.size(),
  //  false,
  //  3.0*d_max, 2,
  //  weight
  //  );
  X.resize(n,3);
  B.resize(n,3);
  I.resize(n,1);
  for ( size_t i=0; i<outputPoints.size(); i++ ) 
  {
    X.row(i) = X5.row(outputPoints[i].i);
    B.row(i) = B5.row(outputPoints[i].i);
    I(i) =         I5(outputPoints[i].i);
  }
}
#endif 
