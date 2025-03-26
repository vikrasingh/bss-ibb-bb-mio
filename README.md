# bss-ibb-bb-mio
This package provides the testing and algorithm files to solve the following **Best Subset Selection (BSS)** problem
```math
\min_{ \beta \in R^p} \|y-X \beta \|_{2}^{2}  \text{ subject to } \|\beta\|_{0}\leq k \text{ , } l\leq \beta \leq u
``` 
where $l$, $u$, $\beta \in R^p$, $y \in R^n$, $X \in R^{n \times p}$, $`\| \cdot \|_{2}`$ is $`l_{2}`$ norm, $`\|\cdot\|_{0}`$ is a pseudo-norm which counts the no. of non-zero entries of $`\beta`$.\
Three algorithms have been implemented to solve the (BSS) problem: IBB, BB, and MIO.

"main" directory
contains the m-files ibb.m, bb.m, mio.m for the three algorithms.

"quad_min" directory
contains the m-files to find the OLS solutions in the original or the reduced space, using the method of parallel tangent for boxed quadratic
minimization, and the method of conjugate gradient for unboxed quadratic minimization.

"topmost" directory
contains the wrapper m-files to generate the example data, call each algorithm, and save the output in a spreadsheet, as well as in MATLAB data files
with the following dependencies
serialRunCall.m -> runLLS.m -> runAnEgLLS.m -> setupForIntvalAlgo.m -> ibb.m, bb.m, mio.m

"test_files" directory
contains m-files to test the three algorithms  
Eg1Od.m to test example 1 in OD case\
Eg1Ud.m to test example 1 in UD case\
Eg2Od.m to test example 2 in OD case\
Eg2Ud.m to test example 2 in UD case\
Eg3Od.m to test example 3 in OD case\
Eg3Ud.m to test example 3 in UD case

                                                                        
