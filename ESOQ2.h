/*
 * ESOQ2.h
 *
 *  Created on: Sep 3, 2013
 *      Author: muzhig
 */

#ifndef ESOQ2_H_
#define ESOQ2_H_

#include <Matrix.h>

void esoq2(Matrix& obs, const Matrix& ref, double& loss, Matrix& quaternion);


#endif /* ESOQ2_H_ */
