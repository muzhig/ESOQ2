ESOQ2.1
=======
Optimal quaternion estimation

[ESOQ-2 Single-Point Algorithm for Fast Optimal Spacecraft Attitude Determination](http://mortari.tamu.edu/Attitude-Estimation/J05.pdf "ESOQ2 Publication of Daniele Mortari")

    ESOQ2.1 - Attitude EStimate OPtimization (F. L. Markley, 7/15/99).
    Variation of Daniele Mortari's ESOQ2 (Paper AAS 97-167, AAS/AIAA Space Flight Mechanics Meeting, Huntsville, AL, February 10-12, 1997), with new singularity-avoidance and lambda_max computation logic.

* input: obs(3,n) - array of n observation unit vectors
* ref(3,n) - array of n reference unit vectors
* wt(n) - row array of n measurement weights (ommited in C++ implementation)

* output: q(4) - optimal quaternion [w,x,y,z]
* loss - optimized value of Wahba's loss function

Implementations:
* *esoq2p1.m* - original matlab source of ESOQ 2.1
* *esoq2p1.py* - python port. Depends on _numpy_
* *ESOQ2.h* & *ESOQ2.cpp* - C++ port, depency: [LineArduino](https://github.com/muzhig/linearduino "linear algebra library for arduino") 


