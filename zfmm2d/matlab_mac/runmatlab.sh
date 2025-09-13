  rm *.o 
  rm *.a 
  rm *.m



  #!/bin/bash
  rm *.a 
  rm *.m
  rm *.o

  # Define system-specific variables for easy modification
  MWRAP_PATH="../../../mwrap"
  MEX_PATH="/Applications/MATLAB_R2024a.app/bin/mex"
  GCC_LIB_PATH="/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/current"
  SRC_PATH="../src"

  # Compiler settings
  FC="gfortran"
  FFLAGS="-c -O3 -march=native -std=legacy -fopenmp"
  AR="ar rcs"

  # List of Fortran source files
  SRC_FILES=(
      "$SRC_PATH/common/corrand.f"
      "$SRC_PATH/common/cumsum.f"
      "$SRC_PATH/common/fmmcommon2d.f"
      "$SRC_PATH/common/prini.f"
      "$SRC_PATH/common/pts_tree2d.f"
      "$SRC_PATH/common/tree_routs2d.f"
      "$SRC_PATH/common/toolbox.f"
      "$SRC_PATH/helmholtz/cdjseval2d.f"
      "$SRC_PATH/helmholtz/h2dcommon.f"
      "$SRC_PATH/helmholtz/h2dterms.f"
      "$SRC_PATH/helmholtz/hank101.f"
      "$SRC_PATH/helmholtz/hank103.f"
      "$SRC_PATH/helmholtz/hndiv2d.f"
      "$SRC_PATH/helmholtz/zh2dtrans.f"
      "$SRC_PATH/helmholtz/zhfmm2d.f"
      "$SRC_PATH/helmholtz/next235.f"
      "$SRC_PATH/laplace/l2dterms.f"
      "$SRC_PATH/laplace/lndiv2d.f"
      "$SRC_PATH/laplace/zl2dtrans.f"
      "$SRC_PATH/laplace/zlfmm2d.f"
  )

  # Compile Fortran files
  for file in "${SRC_FILES[@]}"; do
      obj_file=$(basename "$file" .f).o
      $FC $FFLAGS "$file" -o "$obj_file"
  done

  # Archive Fortran object files into a static library
  $AR zfmm2dfortran.a *.o

  # Run MWrap
  $MWRAP_PATH/mwrap -c99complex -list -mex zhfmm2dmex -mb zhfmm2d.mw
  $MWRAP_PATH/mwrap -c99complex -list -mex zh2devaldirectmex -mb zh2devaldirect.mw
  $MWRAP_PATH/mwrap -c99complex -list -mex zlfmm2dmex -mb zlfmm2d.mw
  $MWRAP_PATH/mwrap -c99complex -list -mex zl2devaldirectmex -mb zl2devaldirect.mw

  $MWRAP_PATH/mwrap -c99complex -mex zhfmm2dmex -c zhfmm2dmex.c zhfmm2d.mw
  $MWRAP_PATH/mwrap -c99complex -mex zh2devaldirectmex -c zh2devaldirectmex.c zh2devaldirect.mw
  $MWRAP_PATH/mwrap -c99complex -mex zlfmm2dmex -c zlfmm2dmex.c zlfmm2d.mw
  $MWRAP_PATH/mwrap -c99complex -mex zl2devaldirectmex -c zl2devaldirectmex.c zl2devaldirect.mw

  # Compile MEX files
  for mex_file in zhfmm2dmex zh2devaldirectmex zlfmm2dmex zl2devaldirectmex; do
      $MEX_PATH -v "$mex_file.c" zfmm2dfortran.a \
          -largeArrayDims -DMWF77_UNDERSCORE1 -output "$mex_file" \
          -lgfortran -L$GCC_LIB_PATH -lgomp
  done





  
  # gfortran -c -O3 -march=native -std=legacy ../src/common/corrand.f -o corrand.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/common/cumsum.f -o cumsum.o
  # gfortran -c -O3 -march=native -std=legacy ../src/common/fmmcommon2d.f -o fmmcommon2d.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/common/prini.f -o prini.o  
  # gfortran -c -O3 -march=native -std=legacy ../src/common/pts_tree2d.f -o pts_tree2d.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/common/toolbox.f -o toolbox.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/common/tree_routs2d.f -o tree_routs2d.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/cdjseval2d.f -o cdjseval2d.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/h2dcommon.f -o h2dcommon.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/h2dterms.f -o h2dterms.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/hank101.f -o hank101.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/hank103.f -o hank103.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/hndiv2d.f -o hndiv2d.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/zh2dtrans.f -o zh2dtrans.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/zhfmm2d.f -o zhfmm2d.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/helmholtz/next235.f -o next235.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/laplace/l2dterms.f -o l2dterms.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/laplace/lndiv2d.f -o lndiv2d.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/laplace/zl2dtrans.f -o zl2dtrans.o 
  # gfortran -c -O3 -march=native -std=legacy ../src/laplace/zlfmm2d.f -o zlfmm2d.o 

  # # Create a static Fortran library
  # ar rcs zfmm2dfortran.a corrand.o cumsum.o fmmcommon2d.o prini.o pts_tree2d.o  \
  # toolbox.o tree_routs2d.o cdjseval2d.o h2dcommon.o h2dterms.o hank101.o hank103.o \
  # hndiv2d.o zh2dtrans.o zhfmm2d.o next235.o l2dterms.o l2dterms.o lndiv2d.o zl2dtrans.o zlfmm2d.o 

  # # Generate MEX wrapper files using MWrap
  # ../../../mwrap/mwrap -c99complex -list -mex zhfmm2dmex -mb zhfmm2d.mw
  # ../../../mwrap/mwrap -c99complex -list -mex zh2devaldirectmex -mb zh2devaldirect.mw
  # ../../../mwrap/mwrap -c99complex -list -mex zlfmm2dmex -mb zlfmm2d.mw
  # ../../../mwrap/mwrap -c99complex -list -mex zl2devaldirectmex -mb zl2devaldirect.mw

  # ../../../mwrap/mwrap -c99complex -mex zhfmm2dmex -c zhfmm2dmex.c zhfmm2d.mw
  # ../../../mwrap/mwrap -c99complex -mex zh2devaldirectmex -c zh2devaldirectmex.c zh2devaldirect.mw
  # ../../../mwrap/mwrap -c99complex -mex zlfmm2dmex -c zlfmm2dmex.c zlfmm2d.mw
  # ../../../mwrap/mwrap -c99complex -mex zl2devaldirectmex -c zl2devaldirectmex.c zl2devaldirect.mw

  # # Compile MEX files using Linux MATLAB Path
  # MATLAB_BIN="/Applications/MATLAB_R2024a.app/bin/mex"
  # GFORTRAN_LIB_PATH="/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/current"

  # $MATLAB_BIN -v zhfmm2dmex.c zfmm2dfortran.a -largeArrayDims -DMWF77_UNDERSCORE1 -output zhfmm2dmex -lgfortran -L$GFORTRAN_LIB_PATH 
  # $MATLAB_BIN -v zh2devaldirectmex.c zfmm2dfortran.a -largeArrayDims -DMWF77_UNDERSCORE1 -output zh2devaldirectmex -lgfortran -L$GFORTRAN_LIB_PATH 
  # $MATLAB_BIN -v zlfmm2dmex.c zfmm2dfortran.a -largeArrayDims -DMWF77_UNDERSCORE1 -output zlfmm2dmex -lgfortran -L$GFORTRAN_LIB_PATH 
  # $MATLAB_BIN -v zl2devaldirectmex.c zfmm2dfortran.a -largeArrayDims -DMWF77_UNDERSCORE1 -output zl2devaldirectmex -lgfortran -L$GFORTRAN_LIB_PATH 



