#include "blue_noise.h"

#ifdef MEX
#include <mex.h>
#endif
#include <igl/cumsum.h>
#include <igl/histc.h>
#include <igl/slice.h>
#include <igl/C_STR.h>
#ifdef MEX
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#endif
#include <cmath>

#ifdef MEX
void mexFunction(
         int          nlhs,
         mxArray      *plhs[],
         int          nrhs,
         const mxArray *prhs[]
         )
{
  using namespace std;
  using namespace igl;
  using namespace igl::matlab;
  using namespace Eigen;
  igl::matlab::MexStream mout;        
  std::streambuf *outbuf = std::cout.rdbuf(&mout);
  Eigen::MatrixXd I;

  mexErrMsgTxt(nrhs>=2,"nrhs should be >= 2");
  mexErrMsgTxt(
    mxIsDouble(prhs[0]) && mxGetM(prhs[0])==1 && mxGetN(prhs[0])==1,
    "should be scalar");
  int n = *mxGetPr(prhs[0]);
  parse_rhs_double(prhs+1,I);
  bool seed_provided = false;
  mxUint32 seed = 0;
  {
    int i = 2;
    while(i<nrhs)
    {
      mexErrMsgTxt(mxIsChar(prhs[i]),
        "Parameter names should be char strings");
      // Cast to char
      const char * name = mxArrayToString(prhs[i]);
      if(strcmp("Seed",name) == 0)
      {
        validate_arg_scalar(i,nrhs,prhs,name);
        // could be uint32
        seed = (mxUint32)*mxGetPr(prhs[++i]);
        seed_provided = true;
      }else
      {
        mexErrMsgTxt(false,
          C_STR("Unsupported parameter: "<<name));
      }
      i++;
    }
  }
  if(seed_provided)
  {
    srand(seed);
  }
#else
int main()
{
  Eigen::MatrixXd I = Eigen::MatrixXd::Zero(256,256);
  const int n = 10000;
#endif

  const int w = I.cols();
  const int h = I.rows();
  Eigen::MatrixXd U(n*5,2);
  Eigen::VectorXd UI(U.rows());
  // Should pre-bias samples toward I... by rejection?
  for(int i = 0;i<U.rows();i++)
  { 
    while(true)
    {
      const double x = ((double)rand()/(double)RAND_MAX)*(w-1);
      const double y = ((double)rand()/(double)RAND_MAX)*(h-1);
      // Nearest neighbor
      //const double ui = I(int(y),int(x));
      // Bilinear interpolation
      double xid;
      const double xf = std::modf(x,&xid);
      int xi = int(xid);
      double yid;
      const double yf = std::modf(y,&yid);
      int yi = int(yid);
      const double ui = 
        (I(yi+0,xi+0)*(1.-yf) + 
         I(yi+1,xi+0)*(   yf))*(1.-xf)+
        (I(yi+0,xi+1)*(1.-yf) + 
         I(yi+1,xi+1)*(   yf))*(   xf);

      const double r = ((double)rand()/(double)RAND_MAX);
      //if(r < 0.01+0.99*(1.0-ui))
      if(r < (1.0-ui))
      {
        U(i,0) = x;
        U(i,1) = y;
        UI(i) = ui;
        break;
      }
    }
  }
  //// This leads to repitition and crashing
  //Eigen::VectorXd UI(U.rows());
  //for(int i = 0;i<U.rows();i++)
  //{
  //  // Nearest neighbor
  //  UI(i) = I(int(U(i,1)),int(U(i,0)));
  //}
  //Eigen::VectorXd C;
  //Eigen::VectorXd A0(UI.size()+1);
  //A0(0) = 0;
  //A0.bottomRightCorner(UI.size(),1) = 1.0-UI.array();
  //// Even faster would be to use the "Alias Table Method"
  //igl::cumsum(A0,1,C);
  //const double Cmax = C(C.size()-1);
  //for(int i = 0;i<C.size();i++) { C(i) = C(i)/Cmax; }
  //const Eigen::VectorXd R = (Eigen::VectorXd::Random(n*5,1).array() + 1.)/2.;
  //Eigen::VectorXi J;
  //igl::histc(R,C,J);
  //igl::slice(Eigen::MatrixXd(U),J,1,U);
  //igl::slice(Eigen::VectorXd(UI),J,1,UI);

  std::vector< cy::Point3f_id > inputPoints(U.rows());
  for ( size_t i=0; i<inputPoints.size(); i++ ) {
      inputPoints[i].x = U(i,0);
      inputPoints[i].y = U(i,1);
      inputPoints[i].z = (U.cols()>=3)?U(i,2):0;
      inputPoints[i].i = i;
  }
  cy::WeightedSampleElimination< cy::Point3f_id, float, 3, int > wse;
  std::vector< cy::Point3f_id > outputPoints(n);
  float area = w*h;
  float d_max = 2 * wse.GetMaxPoissonDiskRadius( 2, outputPoints.size(), area );
  //wse.Eliminate( 
  //  inputPoints.data(), inputPoints.size(), 
  //  outputPoints.data(), outputPoints.size(),
  //  false,
  //  d_max, 2
  //  );
  typedef float FType;
  typedef cy::Point3f_id PointType;
  FType beta  = 0.65;
  FType gamma = 1.5;
  typedef int SIZE_TYPE;
  const auto GetWeightLimitFraction = [beta,gamma]( SIZE_TYPE inputSize, SIZE_TYPE outputSize )->FType
  {
    FType ratio = FType(outputSize) / FType(inputSize);
    return ( 1 - std::pow( ratio, gamma ) ) * beta;
  };
  FType d_min = d_max * GetWeightLimitFraction( inputPoints.size(), outputPoints.size());
  FType alpha = 8;
  const auto weight = 
    [d_min, alpha, &UI](
    PointType const & p0, 
    PointType const & p1, 
    FType d2, 
    FType d_max)
  {
    FType d = p0.i>=0 ? cy::Sqrt(d2)*(3.0-2.0*(UI(p0.i))) : cy::Sqrt(d2);
    //FType d = cy::Sqrt(d2);
    //if ( d < d_min ) d = d_min;
    return std::pow( FType(1) - d/(3.*d_max), alpha );
  };
  wse.Eliminate( 
    inputPoints.data(), inputPoints.size(), 
    outputPoints.data(), outputPoints.size(),
    false,
    d_max, 2,
    weight
    );

  Eigen::MatrixXd X(n,2);
  for ( size_t i=0; i<outputPoints.size(); i++ ) 
  {
    X.row(i) = U.row(outputPoints[i].i);
  }
#ifdef MEX
  switch(nlhs)
  {
    case 3:
      prepare_lhs_double(UI,plhs+2);
    case 2:
      prepare_lhs_double(U,plhs+1);
    case 1:
      prepare_lhs_double(X,plhs+0);
    default:break;
  }
  
  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}
#else
}
#endif
