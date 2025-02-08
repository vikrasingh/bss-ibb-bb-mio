%===============================================================================================
% Authors: Vikram Singh and Min Sun
%
% Copyright (2025): Vikram Singh and Min Sun
% This code is distributed under the terms of the GNU General Public License
% 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire 
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose.
================================================================================================
This package provides the testing and algorithm files to solve the following Best Subset Selection (BSS) problem

 min_{b in R^p} ||y-Xb||_{2}^{2}  subject to ||b||_{0}<=k , l<=b<=u  

where l,u,b are in R^p, y in R^n, X in R^{nxp}, ||.||_{2} is l-2 norm, ||.||_{0} is pseudo-norm which counts the no. of non-zero entries of b
Three algorithms has been implemented to solve (BSS) problem: IBB, BB, and MIO

"main" directory
contains the m-files ibb.m, bb.m, mio.m for the three algorithms.

"quad_min" directory
contains the m-files to find the OLS solutions in the original or the reduced space, using the method of parallel tangent for boxed quadratic
minimization, and the method of conjugate gradient for unboxed quadratic minimization.

"topmost" directory
contains the wrapper m-files to generate the example data, call each algorithm, and save the output in a spreadsheet, as well as in matlab data files
with the following dependencies
serialRunCall.m -> runLLS.m -> runAnEgLLS.m -> setupForIntvalAlgo.m -> ibb.m, bb.m, mio.m

"test_files" directory
contains m-files to test the three algorithms  
p28pairBox7a2c3dOne24Ex20t2kOd2Tm_V9.m to test example 1 in OD case
p28pairBox7a2c3dOne24Ex20t2kUd2Tm_V9.m to test example 1 in UD case
p28pairBox7a2c3dTwo24Ex20t2kOd2Tm_V9.m to test example 2 in OD case
p28pairBox7a2c3dTwo24Ex20t2kUd2Tm_V9.m to test example 2 in UD case
p28pairBox7a2c3dRan24Ex20t2kOd2Tm_V9.m to test example 3 in OD case
p28pairBox7a2c3dRan24Ex20t2kUd2Tm_V9.m to test example 3 in UD case

                                                                        
