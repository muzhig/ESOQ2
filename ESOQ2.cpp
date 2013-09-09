#include <math.h>
#include "ESOQ2.h"
#include <string.h>

void minimax(double* arr, int n, int& imin, int& imax) {
	imax = imin = 0;
	for (int i=0; i<n; i++) {
		if (arr[i] < arr[imin])
			imin = i;
		if (arr[i] > arr[imax])
			imax = i;
	}
}

void esoq2(Matrix& obs, const Matrix& ref, double& loss, Matrix& quaternion) {
	/* *
	 * vectors - matrix of columns, normalized measurement vectors.
	 * ref - matrix of normalized column vectors in reference orientation
	 * loss - output variable for loss metric
	 * quaternion - output matrix for one column, 4D vector.
	 * */
	double lam = obs.n; // number of observations
	bool two_observations = obs.n == 2;

	Matrix& B = obs.dotSelf(ref.transposed()); // 3xN * Nx3 = 3x3

	double trB = B.trace();
	double diag[] = {B(0,0), B(1,1), B(2,2), trB};

	int imin;
	int imax;
	minimax(diag, 4, imin, imax);
	int irot = imin;
	double Bmin = diag[irot];
	int col1;
	int col2;

	if (irot == 0) {
		col1 = 1;
		col2 = 2;
	} else if (irot == 1) {
		col1 = 0;
		col2 = 2;
	} else if (irot == 2) {
		col1 = 0;
		col2 = 1;
	}
	trB = 2 * Bmin - trB;
	for(int i=0; i<3; i++) {
		B(i,col1) *= -1;
		B(i,col2) *= -1;
	}

	double S11 = 2 * B(0, 0);
	double S23 = B(1, 2) + B(2, 1);
	double S22 = 2 * B(1, 1);
	double S31 = B(2, 0) + B(0, 2);
	double S33 = 2 * B(2, 2);
	double S12 = B(0, 1) + B(1, 0);

	double z_[] = {
			B(1, 2) - B(2, 1),
			B(2, 0) - B(0, 2),
			B(0, 1) - B(1, 0)
	};
	B.release(); // free memory
	Matrix z(1, 3, z_);//row

	double z12 = z_[0] * z_[0];
	double z22 = z_[1] * z_[1];
	double z32 = z_[2] * z_[2];



	if (two_observations) {
		double lam0 = lam;
		double trB2 = trB * trB;
		double SzTSz = 0.0;
		{
			double Sz_[] = {S11, S12, S31, S12, S22, S23, S31, S23, S33};
			Matrix Sz = Matrix(3, 3, Sz_).dot(z);
			for (int i = 0; i< Sz.m; i++) {
				SzTSz += Sz(i, 0) * Sz(i, 0);
			}
		}

		double aa = trB2 - S22 * S33 + S23 * S23 - S11 * S33 + S31 * S31 - S22 * S11 + S12 * S12;
		double bb = trB2 + z12 + z22 + z32;
		double c2 = -aa - bb;

		double u = 2 * sqrt(aa * bb - SzTSz);

		lam = (sqrt(u - c2) + sqrt(-u - c2)) * 0.5;
		loss = lam0 - lam;
	}
	double tml = trB - lam;
	double tpl = trB + lam;

	double M11 = tml * (S11 - tpl) - z12;
	double M23 = tml * S23 - z_[1] * z_[2];
	double M22 = tml * (S22 - tpl) - z22;
	double M31 = tml * S31 - z_[2] * z_[0];
	double M33 = tml * (S33 - tpl) - z32;
	double M12 = tml * S12 - z_[0] * z_[1];
	double e_[] = {
		M22 * M33 - M23 * M23,
		M11 * M33 - M31 * M31,
		M11 * M22 - M12 * M12
	};
	{
		double e_abs[] = {
			fabs(e_[0]),
			fabs(e_[1]),
			fabs(e_[2])
		};
		minimax(e_abs, 3, imin, imax);
	}

	if (imax == 0) {
		e_[1] = M31 * M23 - M12 * M33;
		e_[2] = M12 * M23 - M31 * M22;
	} else if (imax == 1) {
		e_[0] = M31 * M23 - M12 * M33;
		e_[2] = M12 * M31 - M11 * M23;
	} else if (imax == 2) {
		e_[0] = M12 * M23 - M31 * M22;
		e_[1] = M12 * M31 - M11 * M23;
	}
	Matrix e(3,1,e_); // column vector
	if (!two_observations) {
		double m1_[] = {M11, M12, M31};
		double m2_[] = {M12, M22, M23};
		double m3_[] = {M31, M23, M33};
		double n1_[] = {(S11 - 2 * lam), S12, S31};
		double n2_[] = {S12, (S22 - 2 * lam), S23};
		double n3_[] = {S31, S23, (S33 - 2 * lam)};
		Matrix m1(1,3,m1_);
		Matrix m2(1,3,m2_);
		Matrix m3(1,3,m3_);
		Matrix n1(1,3,n1_);
		Matrix n2(1,3,n2_);
		Matrix n3(1,3,n3_);

		Matrix& a = (imax==0)?(m2):(imax==1)?(m3):(m1);
		Matrix& b = (imax==0)?(n3):(imax==1)?(n1):(n2);
		Matrix& c = (imax==0)?(m3):(imax==1)?(m1):(m2);
		Matrix& d = (imax==0)?(n2):(imax==1)?(n3):(n1);
		Matrix v = a.cross(b) - c.cross(d); // row vector

		Matrix& m = (imax==0)?(m1):(imax==1)?(m2):(m3); // row
		Matrix& n = (imax==0)?(n1):(imax==1)?(n2):(n3); // row
		v.transpose(); // column
		loss = -(m.dot(e).get(0,0)) / (n.dot(e) + m.dot(v)).get(0, 0);
		tml = tml + loss;
		v *= loss;
		e += v;
	}
	double q[4];
	q[3] = -(z.dot(e).get(0,0));
	e *= tml;
	q[0] = e(0,0);
	q[1] = e(1,0);
	q[2] = e(2,0);



	if (irot == 0){
		double tmp[] = {-q[0], q[3], -q[2], q[1]};
		memcpy(&q, tmp, sizeof(tmp));
	} else if (irot == 1) {
		double tmp[] = {-q[1], q[2], q[3], -q[0]};
		memcpy(&q, tmp, sizeof(tmp));
	} else if (irot == 2) {
		double tmp[] = {-q[2], -q[1], q[0], q[3]};
		memcpy(&q, tmp, sizeof(tmp));
	}

	quaternion = Matrix(4, 1, q);
	quaternion.normalize();
}
