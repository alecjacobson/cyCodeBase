#!/bin/bash
/*/../bin/ls > /dev/null
# BEGIN BASH SCRIPT
export PS4=""
set -o xtrace
source ~/.profile
echo $DYLD_FALLBACK_LIBRARY_PATH
TEMP="$0.cpp"
printf "//" | cat - $0 >$TEMP
#/usr/local/opt/llvm/bin/clang++ -fopenmp=libiomp5 -I /usr/local/include/libiomp/ -liomp5 -g -O0 -std=c++11 -ferror-limit=4 -o .main $TEMP -msse4.2 \
export LIBIGL=/usr/local/libigl/
export STAN=/Users/ajx/Dropbox/math/
#export LIBIGL=/Users/ajx/Dropbox/interactive-segmentation/libigl/
#-DIGL_STATIC_LIBRARY -L"$LIBIGL"/lib -ligl -ligl_matlab -ligl_opengl -ligl_opengl_glfw -ligl_cgal \
#g++-4.8   -g -O0 -std=c++11 -fopenmp -o .main $TEMP -DNDEBUG -msse4.2 -m32 \
#clang++ -g -O0 -std=c++14 -o .main $TEMP -msse4.2 \
clang++ -g -O3 -std=c++11 -o .main $TEMP -DNDEBUG -msse4.2 \
  -I. \
  -I"$STAN" \
  -I/Applications/MATLAB_R2019a.app/extern/include -L/Applications/MATLAB_R2019a.app//bin/maci64/ -leng -lmat -lmex -lmx \
  -I"$LIBIGL"/external/glad/include/ \
  -I"$LIBIGL"/external/eigen \
  -I"$LIBIGL"/external/glfw/include \
  -I"$LIBIGL"/include \
  -I"$LIBIGL"/external/AntTweakBar/include \
  -I"$LIBIGL"/external/AntTweakBar/src \
  -I"$LIBIGL"/external/tetgen \
  -I"$LIBIGL"/external/tinyxml2/ \
  -I"$LIBIGL"/external/Singular_Value_Decomposition/ \
  -framework Carbon -framework QuartzCore -framework IOKit \
  -I"$LIBIGL"/external/nanogui/include \
  -I"$LIBIGL"/external/nanogui/ext \
  -I"$LIBIGL"/external/nanogui/ext/nanovg/src \
  -I"$LIBIGL"/external/embree/ \
  -I"$LIBIGL"/external/embree/include \
  -framework OpenGL \
  -framework AppKit \
  -I/Users/ajx/Dropbox/ \
  -I/usr/local/libigl/external/triangle \
  -L"$LIBIGL"/build-release/ -lglfw3 \
  -L"$LIBIGL"/build-release/ -lglad \
  -L"$LIBIGL"/build-release/embree -lembree3 -lsys -lmath -llexers -lsimd -lembree_avx -lembree_avx2 -lembree_sse42 -ltasking \
&& /Applications/Xcode.app/Contents/Developer/usr/bin/lldb -b -o r ./.main -- "$@"
#&& ./.main "$@"
#-L/usr/local/lib -lCGAL -lCGAL_Core -lgmp -lmpfr -lboost_thread-mt -lboost_system-mt \
#-L"$LIBIGL"/lib/ -lglad \
#-L/usr/local/lib -lboost_thread-mt -lboost_system-mt \
#-L/usr/local/lib -lboost_program_options-mt -fno-math-errno \
# /usr/local/libigl/lib/libtriangle.a \
#-I~/Documents/eigen/ \
#rm -f .main
#rm -f $TEMP
# END BASH SCRIPT
exit
*/

#include "blue_noise.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <set>
#include <queue>
#include <list>
#include <vector>
#include <time.h>
#include <random>
#include <chrono>

int main(int argc, char * argv[])
{
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::read_triangle_mesh(argv[1],V,F);
  igl::opengl::glfw::Viewer v;
  v.data().set_mesh(V,F);
  int n = strtol(argv[2], NULL, 10);
  Eigen::MatrixXd B,P;
  Eigen::VectorXi FI;
  igl::blue_noise(n,V,F,B,FI,P);

  v.data().set_points(P,Eigen::RowVector3d(0.1,0.1,0.1));
  v.data().show_lines = false;
  v.data().point_size = 5;
  v.launch();
}
