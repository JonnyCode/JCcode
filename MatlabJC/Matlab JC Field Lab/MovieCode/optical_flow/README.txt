Only "example.m" and "optical_flow.m" have been properly documented.
The code is an implementation of the phase-based optical flow algorithm
described in

Gautama, T. and Van Hulle, M.M.  (2002).  A Phase-based Approach to the Estimation
of the Optical Flow Field Using Spatial Filtering, IEEE Trans. Neural Networks,
13(5), 1127--1136.

I have modified the algorithm by leaving out the Goal Programming
Network and have replaced it by an analytical solution to the
intersection-of-constraints (the original algorithm is still available
as optical_flow_orig.m)

The code has been tested on Solaris (Matlab 5.3 R11) and Linux (Matlab
6.0 R12), and it works fine.

To check the performance, run the example.m script:

Coverage: 32.02 %
Angular error (full): 4.1771 deg (3.2145 deg)

If the errors deviate from these values, there is probably something
wrong.  Feel free to hack in or, alternatively, to email me.

For more information, contact:
temu@neuro.kuleuven.ac.be
simone.neuro.kuleuven.ac.be/temu


example.m		example script

optical_flow.m		optical flow algorithm

eval_flow.m		evaluate obtained OF (only if desired OF is available)

vis_flow.m		visualise OF

read_sequence.m		read image sequence
read_flow.m		read OF (Barron's format)
rastread.m		read Sun raster file

yosemite/		contains yosemite sequence + desired flow

% OBSOLETE:

optical_flow_orig.m	original optical flow algorithm
gpn_ada.m		Goal Programming Network (calls the mex function)
gpn_ada_mex.c			source
gpn_ada_mex.mexglx		precompiled for linux
gpn_ada_mex.mexsol		precompiled for Solaris
gpn_ada_mex.dll			precompiled for windows
gpn_ada.m.windows               non-mex version (slow, slow, slow)



