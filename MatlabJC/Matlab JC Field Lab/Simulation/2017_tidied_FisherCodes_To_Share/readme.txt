These scripts implement the theoretical Fisher Info calculations from Fig. 7 of Zylberberg, Cafaro, Turner, et al. Neuron 2016 “Direction selective circuits shape noise to ensure a precise population code”

Please contact Joel Zylberberg (joel.zylberberg@ucdenver.edu) to report any bugs.

The script TC.m contains the tuning curve shapes. It is called by do_FI_calc_stimdep.m, which generates responses from the tuning curves, and computes the covariance matrices using Eq. 2 of the paper (assuming Poisson noise). It then computes the Fisher info for the case of stim dependent correlations, and for the case of “matched” constant correlations.

There are two looper scripts that repeat the calculation for different population sizes and correlation strengths. looper_FI_HOMOG.m specifies homogeneous tuning curves, whereas looper_FI_HETEROG.m specifies heterogeneous ones.

To Run for heterogeneous tuning curves (note, this may take several minutes to run):
>> looper_FI_HETEROG.m

To Run for homogeneous tuning curves (note, this may take several minutes to run):
>> looper_FI_HOMOG.m