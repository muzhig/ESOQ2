#include <Matrix.h>
#include "ESOQ2.h"
#include <iostream>
double sina = 0.237557359899;
double cosa = -0.971373512485;

double ref_[] = {
		1.0, 0.0,
		0.0, 0.0,
		0.0, 1.0
};
Matrix ref(3,2,ref_);

void test_zero_rotation() {

	double obs_[] = {
			1.0, 0.0,
			0.0, 0.0,
			0.0, 1.0
	};
	Matrix obs(3,2,obs_);
	Matrix q;
	double loss;
	esoq2(obs, ref, loss, q);

	std::cout << loss << " (" << q.get(0,0) << ", " << q.get(1,0) << ", " << q.get(2,0) << ", " << q.get(3,0) << ")\n";

}

void test_rotation_Z_90() {

	double obs_[] = {
			0.0, 0.0,
			1.0, 0.0,
			0.0, 1.0
	};
	Matrix obs(3,2,obs_);
	Matrix q;
	double loss;
	esoq2(obs, ref, loss, q);

	std::cout << loss << " (" << q.get(0,0) << ", " << q.get(1,0) << ", " << q.get(2,0) << ", " << q.get(3,0) << ")\n";

}

void test_rotation_X_90() {

	double obs_[] = {
			1.0, 0.0,
			0.0, -1.0,
			0.0, 0.0
	};
	Matrix obs(3,2,obs_);
	Matrix q;
	double loss;
	esoq2(obs, ref, loss, q);
	//0.70710678, -0.70710678, -0.        ,  0.
	std::cout << loss << " (" << q.get(0,0) << ", " << q.get(1,0) << ", " << q.get(2,0) << ", " << q.get(3,0) << ")\n";

}

void test_rotation_Y_90() {

	double obs_[] = {
			0.0, 1.0,
			0.0, 0.0,
			-1.0, 0.0
	};
	Matrix obs(3,2,obs_);
	Matrix q;
	double loss;
	esoq2(obs, ref, loss, q);

	std::cout << loss << " (" << q.get(0,0) << ", " << q.get(1,0) << ", " << q.get(2,0) << ", " << q.get(3,0) << ")\n";

}


int main() {
	test_zero_rotation();
	test_rotation_X_90();
	test_rotation_Y_90();
	test_rotation_Z_90();
	return 0;
}
