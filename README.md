ESOQ2.1
=======
Optimal quaternion estimation given measured vectors and vectors in reference system.

prof. Daniele Mortari's publication (pdf):
["ESOQ-2 Single-Point Algorithm for Fast Optimal Spacecraft Attitude Determination"](http://mortari.tamu.edu/Attitude-Estimation/J05.pdf)



* input: obs(3,n) - array of n observation unit vectors
* ref(3,n) - array of n reference unit vectors
* wt(n) - row array of n measurement weights (ommited in C++ implementation)
* output: q(4) - optimal quaternion [w,x,y,z]
* loss - optimized value of Wahba's loss function


Implementations:
* *esoq2p1.m* - original matlab source provided by prof. Daniele Mortari
* *esoq2p1.py* - python port. Depends on _numpy_
* *ESOQ2.h* & *ESOQ2.cpp* - C++ port, Arduino compatible, depends on [LineArduino](https://github.com/muzhig/linearduino "linear algebra library for arduino") 


