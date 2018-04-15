

Dependencies:

Sparse null space and orthogonal
Available:
https://se.mathworks.com/matlabcentral/fileexchange/27550-sparse-null-space-and-orthogonal
Notes:
The Matlab function 'spnull' is required.

minFunc
Available:
https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html
Notes:
The Matlab function 'WolfeLineSearch' is required.

Mosek
Available:
https://www.mosek.com/


SPOTless
Available:
https://github.com/spot-toolbox/spotless
Notes:
please replace the following functions from the standard SPOTless repository with the equivalent functions from SpecializedDemo/SpotlessFunctions (i.e. this repository),

1) spotopt/@spotprogsol/spotprogsol - corrects a (possible) bug that expects Gram matrices when none are defined.

2) spotopt/@spotsosprog/spotsosprog - makes the monomial generating function public.

3) @msspoly/msspoly - adds 'prod' functionality for matrices.

