  rm *.o 
  rm *.a 
  rm *.m
  
  gfortran -c -O3 -march=native -std=legacy ../src/common/corrand.f -o corrand.o 
  gfortran -c -O3 -march=native -std=legacy ../src/common/cumsum.f -o cumsum.o
  gfortran -c -O3 -march=native -std=legacy ../src/common/fmmcommon2d.f -o fmmcommon2d.o 
  gfortran -c -O3 -march=native -std=legacy ../src/common/prini.f -o prini.o  
  gfortran -c -O3 -march=native -std=legacy ../src/common/pts_tree2d_isep.f -o pts_tree2d_isep.o 
  gfortran -c -O3 -march=native -std=legacy ../src/common/toolbox.f -o toolbox.o 
  gfortran -c -O3 -march=native -std=legacy ../src/common/tree_routs2d_isep.f -o tree_routs2d_isep.o 
  gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/cdjseval2d.f -o cdjseval2d.o 
  gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/h2dcommon.f -o h2dcommon.o 
  gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/h2dterms.f -o h2dterms.o 
  gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/hank101.f -o hank101.o 
  gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/hank103.f -o hank103.o 
  gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/hndiv2d.f -o hndiv2d.o 
  gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/zh2dtrans.f -o zh2dtrans.o 
  gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/zhfmm2d.f -o zhfmm2d.o 
  gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/next235.f -o next235.o 
  gfortran -c -O3 -march=native -std=legacy ../src/laplace/l2dterms.f -o l2dterms.o 
  gfortran -c -O3 -march=native -std=legacy ../src/laplace/lndiv2d.f -o lndiv2d.o 
  gfortran -c -O3 -march=native -std=legacy ../src/laplace/zl2dtrans.f -o zl2dtrans.o 
  gfortran -c -O3 -march=native -std=legacy ../src/laplace/zlfmm2d.f -o zlfmm2d.o 

  # Create a static Fortran library
  ar rcs zfmm2dfortran.a corrand.o cumsum.o fmmcommon2d.o prini.o pts_tree2d_isep.o  \
  toolbox.o tree_routs2d_isep.o cdjseval2d.o h2dcommon.o h2dterms.o hank101.o hank103.o \
  hndiv2d.o zh2dtrans.o zhfmm2d.o next235.o l2dterms.o l2dterms.o lndiv2d.o zl2dtrans.o zlfmm2d.o 

  # Generate MEX wrapper files using MWrap
  ../../mwrap/mwrap -c99complex -list -mex zhfmm2dmex -mb zhfmm2d.mw
  ../../mwrap/mwrap -c99complex -list -mex zh2devaldirectmex -mb zh2devaldirect.mw
  ../../mwrap/mwrap -c99complex -list -mex zlfmm2dmex -mb zlfmm2d.mw
  ../../mwrap/mwrap -c99complex -list -mex zl2devaldirectmex -mb zl2devaldirect.mw

  ../../mwrap/mwrap -c99complex -mex zhfmm2dmex -c zhfmm2dmex.c zhfmm2d.mw
  ../../mwrap/mwrap -c99complex -mex zh2devaldirectmex -c zh2devaldirectmex.c zh2devaldirect.mw
  ../../mwrap/mwrap -c99complex -mex zlfmm2dmex -c zlfmm2dmex.c zlfmm2d.mw
  ../../mwrap/mwrap -c99complex -mex zl2devaldirectmex -c zl2devaldirectmex.c zl2devaldirect.mw

  # Compile MEX files using Linux MATLAB Path
  MATLAB_BIN="/usr/local/MATLAB/R2023a/bin/mex"
  GFORTRAN_LIB_PATH="/usr/lib/gcc/x86_64-linux-gnu/"

  $MATLAB_BIN -v zhfmm2dmex.c zfmm2dfortran.a -largeArrayDims -DMWF77_UNDERSCORE1 -output zhfmm2dmex -lgfortran -L$GFORTRAN_LIB_PATH 
  $MATLAB_BIN -v zh2devaldirectmex.c zfmm2dfortran.a -largeArrayDims -DMWF77_UNDERSCORE1 -output zh2devaldirectmex -lgfortran -L$GFORTRAN_LIB_PATH 
  $MATLAB_BIN -v zlfmm2dmex.c zfmm2dfortran.a -largeArrayDims -DMWF77_UNDERSCORE1 -output zlfmm2dmex -lgfortran -L$GFORTRAN_LIB_PATH 
  $MATLAB_BIN -v zl2devaldirectmex.c zfmm2dfortran.a -largeArrayDims -DMWF77_UNDERSCORE1 -output zl2devaldirectmex -lgfortran -L$GFORTRAN_LIB_PATH 




# #!/bin/bash

# # Define system-specific variables for easy modification
# MWRAP_PATH="../../mwrap"
# MEX_PATH="/usr/local/MATLAB/R2023a/bin/mex"
# GCC_LIB_PATH="/usr/lib/gcc/x86_64-linux-gnu/"
# SRC_PATH="../src"

# # Compiler settings
# FC="gfortran"
# FFLAGS="-c -O3 -march=native -std=legacy"
# AR="ar rcs"

# # List of Fortran source files
# SRC_FILES=(
#     "$SRC_PATH/common/corrand.f"
#     "$SRC_PATH/common/cumsum.f"
#     "$SRC_PATH/common/dfft.f"
#     "$SRC_PATH/common/fmmcommon.f"
#     "$SRC_PATH/common/prini.f"
#     "$SRC_PATH/common/legeexps.f"
#     "$SRC_PATH/common/pts_tree3d.f"
#     "$SRC_PATH/common/tree_routs3d.f"
#     "$SRC_PATH/common/yrecursion.f"
#     "$SRC_PATH/common/zrotate.f"
#     "$SRC_PATH/common/ztoolbox.f"
#     "$SRC_PATH/helmholtz/besseljs3d.f"
#     "$SRC_PATH/helmholtz/h3dcommon.f"
#     "$SRC_PATH/helmholtz/h3dterms.f"
#     "$SRC_PATH/helmholtz/hank101.f"
#     "$SRC_PATH/helmholtz/hndiv_fast.f"
#     "$SRC_PATH/helmholtz/zh3dtrans.f"
#     "$SRC_PATH/helmholtz/zhfmm3d.f"
#     "$SRC_PATH/laplace/l3dterms.f"
#     "$SRC_PATH/laplace/lndiv_fast.f"
#     "$SRC_PATH/laplace/zl3dtrans.f"
#     "$SRC_PATH/laplace/zlfmm3d.f"
# )

# # Compile Fortran files
# for file in "${SRC_FILES[@]}"; do
#     obj_file=$(basename "$file" .f).o
#     $FC $FFLAGS "$file" -o "$obj_file"
# done

# # Archive Fortran object files into a static library
# $AR zfmm3dfortran.a *.o

# # Run MWrap to generate wrapper source files
# $MWRAP_PATH/mwrap -c99complex -list -mex zhfmm3dmex -mb zhfmm3d.mw
# $MWRAP_PATH/mwrap -c99complex -list -mex zh3devaldirectmex -mb zh3devaldirect.mw
# $MWRAP_PATH/mwrap -c99complex -list -mex zlfmm3dmex -mb zlfmm3d.mw
# $MWRAP_PATH/mwrap -c99complex -list -mex zl3devaldirectmex -mb zl3devaldirect.mw

# $MWRAP_PATH/mwrap -c99complex -mex zhfmm3dmex -c zhfmm3dmex.c zhfmm3d.mw
# $MWRAP_PATH/mwrap -c99complex -mex zh3devaldirectmex -c zh3devaldirectmex.c zh3devaldirect.mw
# $MWRAP_PATH/mwrap -c99complex -mex zlfmm3dmex -c zlfmm3dmex.c zlfmm3d.mw
# $MWRAP_PATH/mwrap -c99complex -mex zl3devaldirectmex -c zl3devaldirectmex.c zl3devaldirect.mw

# # Compile MEX files
# for mex_file in zhfmm3dmex zh3devaldirectmex zlfmm3dmex zl3devaldirectmex; do
#     $MEX_PATH -v "$mex_file.c" zfmm3dfortran.a \
#         -largeArrayDims -DMWF77_UNDERSCORE1 -output "$mex_file" \
#         -lgfortran -L$GCC_LIB_PATH 
# done