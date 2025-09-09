#!/bin/bash

rm *.a 
rm *.o
rm *.m

# Define system-specific variables for easy modification
MWRAP_PATH="../../../mwrap"
MEX_PATH="/usr/local/MATLAB/R2023a/bin/mex"
GCC_LIB_PATH="/usr/lib/gcc/x86_64-linux-gnu/"
SRC_PATH="../src"

# Compiler settings
FC="gfortran"
FFLAGS="-c -O3 -march=native -std=legacy -fopenmp"
AR="ar rcs"

# List of Fortran source files
SRC_FILES=(
    "$SRC_PATH/common/corrand.f"
    "$SRC_PATH/common/cumsum.f"
    "$SRC_PATH/common/dfft.f"
    "$SRC_PATH/common/fmmcommon.f"
    "$SRC_PATH/common/prini.f"
    "$SRC_PATH/common/legeexps.f"
    "$SRC_PATH/common/pts_tree3d.f"
    "$SRC_PATH/common/tree_routs3d.f"
    "$SRC_PATH/common/yrecursion.f"
    "$SRC_PATH/common/zrotate.f"
    "$SRC_PATH/common/ztoolbox.f"
    "$SRC_PATH/helmholtz/besseljs3d.f"
    "$SRC_PATH/helmholtz/h3dcommon.f"
    "$SRC_PATH/helmholtz/h3dterms.f"
    "$SRC_PATH/helmholtz/hank101.f"
    "$SRC_PATH/helmholtz/hndiv_fast.f"
    "$SRC_PATH/helmholtz/zh3dtrans.f"
    "$SRC_PATH/helmholtz/zhfmm3d.f"
    "$SRC_PATH/laplace/l3dterms.f"
    "$SRC_PATH/laplace/lndiv_fast.f"
    "$SRC_PATH/laplace/zl3dtrans.f"
    "$SRC_PATH/laplace/zlfmm3d.f"
)

# Compile Fortran files
for file in "${SRC_FILES[@]}"; do
    obj_file=$(basename "$file" .f).o
    $FC $FFLAGS "$file" -o "$obj_file"
done

# Archive Fortran object files into a static library
$AR zfmm3dfortran.a *.o

# Run MWrap to generate wrapper source files
$MWRAP_PATH/mwrap -c99complex -list -mex zhfmm3dmex -mb zhfmm3d.mw
$MWRAP_PATH/mwrap -c99complex -list -mex zh3devaldirectmex -mb zh3devaldirect.mw
$MWRAP_PATH/mwrap -c99complex -list -mex zlfmm3dmex -mb zlfmm3d.mw
$MWRAP_PATH/mwrap -c99complex -list -mex zl3devaldirectmex -mb zl3devaldirect.mw

$MWRAP_PATH/mwrap -c99complex -mex zhfmm3dmex -c zhfmm3dmex.c zhfmm3d.mw
$MWRAP_PATH/mwrap -c99complex -mex zh3devaldirectmex -c zh3devaldirectmex.c zh3devaldirect.mw
$MWRAP_PATH/mwrap -c99complex -mex zlfmm3dmex -c zlfmm3dmex.c zlfmm3d.mw
$MWRAP_PATH/mwrap -c99complex -mex zl3devaldirectmex -c zl3devaldirectmex.c zl3devaldirect.mw

# Compile MEX files
for mex_file in zhfmm3dmex zh3devaldirectmex zlfmm3dmex zl3devaldirectmex; do
    $MEX_PATH -v "$mex_file.c" zfmm3dfortran.a \
        -largeArrayDims -DMWF77_UNDERSCORE1 -output "$mex_file" \
        -lgfortran -L$GCC_LIB_PATH -lgomp
done


  # # Compile Fortran object files
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/common/corrand.f -o corrand.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/common/cumsum.f -o cumsum.o
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/common/dfft.f -o dfft.o
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/common/fmmcommon.f -o fmmcommon.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/common/prini.f -o prini.o  
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/common/legeexps.f -o legeexps.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/common/pts_tree3d.f -o pts_tree3d.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/common/tree_routs3d.f -o tree_routs3d.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/common/yrecursion.f -o yrecursion.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/common/zrotate.f -o zrotate.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/common/ztoolbox.f -o ztoolbox.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/helmholtz/besseljs3d.f -o besseljs3d.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/helmholtz/h3dcommon.f -o h3dcommon.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/helmholtz/h3dterms.f -o h3dterms.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/helmholtz/hank101.f -o hank101.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/helmholtz/hndiv_fast.f -o hndiv_fast.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/helmholtz/zh3dtrans.f -o zh3dtrans.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/helmholtz/zhfmm3d.f -o zhfmm3d.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/laplace/l3dterms.f -o l3dterms.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/laplace/lndiv_fast.f -o lndiv_fast.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/laplace/zl3dtrans.f -o zl3dtrans.o 
  # gfortran -c -O3 -march=native -std=legacy -fopenmp ../src/laplace/zlfmm3d.f -o zlfmm3d.o 

  # # Create a static Fortran library
  # ar rcs zfmm3dfortran.a corrand.o cumsum.o dfft.o fmmcommon.o prini.o legeexps.o pts_tree3d.o  \
  # tree_routs3d.o yrecursion.o zrotate.o ztoolbox.o besseljs3d.o h3dcommon.o h3dterms.o  \
  # hank101.o hndiv_fast.o zh3dtrans.o zhfmm3d.o l3dterms.o lndiv_fast.o zl3dtrans.o zlfmm3d.o 

  # # Generate MEX wrapper files using MWrap
  # ../../mwrap/mwrap -c99complex -list -mex zhfmm3dmex -mb zhfmm3d.mw
  # ../../mwrap/mwrap -c99complex -list -mex zh3devaldirectmex -mb zh3devaldirect.mw
  # ../../mwrap/mwrap -c99complex -list -mex zlfmm3dmex -mb zlfmm3d.mw
  # ../../mwrap/mwrap -c99complex -list -mex zl3devaldirectmex -mb zl3devaldirect.mw

  # ../../mwrap/mwrap -c99complex -mex zhfmm3dmex -c zhfmm3dmex.c zhfmm3d.mw
  # ../../mwrap/mwrap -c99complex -mex zh3devaldirectmex -c zh3devaldirectmex.c zh3devaldirect.mw
  # ../../mwrap/mwrap -c99complex -mex zlfmm3dmex -c zlfmm3dmex.c zlfmm3d.mw
  # ../../mwrap/mwrap -c99complex -mex zl3devaldirectmex -c zl3devaldirectmex.c zl3devaldirect.mw

  # # Compile MEX files using Linux MATLAB Path
  # MATLAB_BIN="/usr/local/MATLAB/R2023a/bin/mex"
  # GFORTRAN_LIB_PATH="/usr/lib/gcc/x86_64-linux-gnu/"

  # $MATLAB_BIN -v zhfmm3dmex.c zfmm3dfortran.a -largeArrayDims -DMWF77_UNDERSCORE1 -output zhfmm3dmex -lgfortran -L$GFORTRAN_LIB_PATH -lgomp
  # $MATLAB_BIN -v zh3devaldirectmex.c zfmm3dfortran.a -largeArrayDims -DMWF77_UNDERSCORE1 -output zh3devaldirectmex -lgfortran -L$GFORTRAN_LIB_PATH -lgomp
  # $MATLAB_BIN -v zlfmm3dmex.c zfmm3dfortran.a -largeArrayDims -DMWF77_UNDERSCORE1 -output zlfmm3dmex -lgfortran -L$GFORTRAN_LIB_PATH -lgomp
  # $MATLAB_BIN -v zl3devaldirectmex.c zfmm3dfortran.a -largeArrayDims -DMWF77_UNDERSCORE1 -output zl3devaldirectmex -lgfortran -L$GFORTRAN_LIB_PATH -lgomp