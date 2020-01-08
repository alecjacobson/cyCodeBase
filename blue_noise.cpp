#include "blue_noise.h"

#include <mex.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#undef assert
#define assert( isOK ) ( (isOK) ? (void)0 : (void) ::mexErrMsgTxt(C_STR(__FILE__<<":"<<__LINE__<<": failed assertion `"<<#isOK<<"'"<<std::endl) ) )

#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>

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
  Eigen::MatrixXd V,X,B;
  Eigen::MatrixXi F;
  Eigen::VectorXi I;

  mexErrMsgTxt(nrhs>=3,"nrhs should be >= 3");
  mexErrMsgTxt(
    mxIsDouble(prhs[0]) && mxGetM(prhs[0])==1 && mxGetN(prhs[0])==1,
    "should be scalar");
  int n = *mxGetPr(prhs[0]);
  mexErrMsgTxt(
    mxIsDouble(prhs[1]) && mxGetN(prhs[1])==3,
    "V should be #V by 3 ");
  parse_rhs_double(prhs+1,V);
  parse_rhs_index(prhs+2,F);
  bool seed_provided = false;
  mxUint32 seed = 0;
  {
    int i = 3;
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


  igl::blue_noise(n,V,F,B,I,X);

  switch(nlhs)
  {
    case 3:
      prepare_lhs_double(B,plhs+2);
    case 2:
      prepare_lhs_index(I,plhs+1);
    case 1:
      prepare_lhs_double(X,plhs+0);
    default:break;
  }
  
  // Restore the std stream buffer Important!
  std::cout.rdbuf(outbuf);
  return;
}
