# re3q3 - Robust 3Q3 solver
This is a re-implementation of the 3Q3 solver (three quadratics in three unknowns) from *Kukelova et al., Efficient Intersection of Three Quadrics and Applications in Computer Vision,  CVPR 2016*. The code is adapted from Jan Heller's original implementation. It also includes the additional tricks for improving stability based on choosing the elimination variable at runtime (see *Zhou et al., A Stable Algebraic Camera Pose Estimation for Minimal Configurations of 2D/3D Point and Line Correspondences, ACCV 2018*). Additionally we do a random affine change of variables to handle further degeneracies.

The function takes a 3x10 matrix corresponding to the coefficients of the three quadratics. The order of the monomials is
> x^2, xy, xz, y^2, yz, z^2, x, y, z, 1.0;

There is also a mex-interface for matlab in re3q3_mex.cpp which can be compiled by running
> mex('-I/usr/include/eigen3','re3q3_mex.cpp')

TODO:
* Replace companion matrix solver with sturm-sequences or similar root-bracketing method
* Add more examples.