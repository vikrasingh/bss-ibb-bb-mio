# bss-ibb-bb-mio
This package provides the testing and algorithm files to solve the following **Best Subset Selection (BSS)** problem
```math
\min_{ \beta \in R^p} f(\beta):=\|y-X \beta \|_{2}^{2} \text{ subject to } \|\beta\|_{0}\leq k \text{ , } l\leq \beta \leq u
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

```matlab
// Pseudo-Code for ibb.m
FUNCTION ibb
  INPUT p n k y X
  OUTPUT fopt xopt 

  INITIALIZE L = Root Node
  INITIALIZE xopt = Initial Feasible Point
  INITIALIZE fopt = f(xopt)
  INITIALZIZE nb = 1  // number of nodes in L

  WHILE nb > 0
    SET Y = Select a Node from L  // choose a node to branch on
    SET nb = nb - 1 
    FIND V1 V2  // child nodes of Y
    FOR i FROM 1 TO 2
        IF Vi satisfies any of the deletion conditions
           DISCARD Vi
           CONTINUE with next iteration of FOR loop

        SET xtilde = Feasible point of Vi

        IF f(xtilde) < fxopt
            UPDATE xopt = xtilde
            UPDATE fxopt = f(xtilde)

        IF fxopt < lb f(Vi)  // If the current best solution is better than the lower bound f(Vi)
           DISCARD Vi
           CONTINUE with next iteration of FOR loop

        Add Vi to the list L
        SET nb = nb + 1

     END FOR

  END WHILE

END FUNCTION
                                                                        
