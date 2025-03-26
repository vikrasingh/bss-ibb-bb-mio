# bss-ibb-bb-mio
This package provides the testing and algorithm files to solve the following **Best Subset Selection (BSS)** problem
```math
\min_{ \beta \in R^p} \|y-X \beta \|_{2}^{2}  \text{ subject to } \|\beta\|_{0}\leq k \text{ , } l\leq \beta \leq u
``` 
where $l$, $u$, $\beta \in R^p$, $y \in R^n$, $X \in R^{n \times p}$, $\| \cdot \|_{2}$ is $l_{2}$ norm, $\|\cdot\|_{0}$ is pseudo-norm which counts the no. of non-zero entries of $\beta$.\
Three algorithms have been implemented to solve the (BSS) problem: IBB, BB, and MIO

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

                                                                        
