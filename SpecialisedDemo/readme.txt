

Dependencies:

Sparse null space and orthogonal
Available:
https://se.mathworks.com/matlabcentral/fileexchange/27550-sparse-null-space-and-orthogonal
Notes:
The Matlab function 'spnull.m' is required.

minFunc
Available:
https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html
Notes:
The Matlab function 'WolfeLineSearch.m' is required.

SPOTless
Available:
https://github.com/spot-toolbox/spotless
Notes:
please replace the following functions from the standard SPOTless repository with the modified functions of the same name found in SpecializedDemo/SpotlessFunctions (i.e. this repository),

1) spotopt/@spotprogsol/spotprogsol.m - corrects a (possible) bug that expects Gram matrices when none are defined.

2) spotopt/@spotsosprog/spotsosprog.m - makes the monomial generating function public.

3) @msspoly/msspoly.m - adds 'prod' functionality for matrices.

SDP Solvers:
Note that the function SpecializedDemo/SpotlessFunctions/lr_init_lmiPLusLin_spot.m requires a SDP solver compatible with SPOTless to solve the feasibility problem that initializes the decision variables in the specialized algorithm.
In all experiments, we used the solver MOSEK, which is available at:
https://www.mosek.com/
To use other solvers, please modify line 35 of lr_init_lmiPLusLin_spot.m, e.g., to use the 'SeDuMi' solver, replace @spot_mosek with @spot_sedumi, c.f., the SPOTless documentation for further details.
