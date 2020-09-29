
#include "mex.h"
#include <Eigen/Dense>
#include "re3q3/re3q3.h"

void printUsage() {
	mexPrintf("[x] = re3q3(C);\n");
	mexPrintf("    coefficient order is [ x^2, xy, xz, y^2, yz, z^2, x, y, z, 1.0 ]\n");
}

void mexFunction(int nlhs,mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	if (nrhs != 1) {
		printUsage();
		mexErrMsgTxt("Please, specify 3 x 10N coefficient matrix");
	}
	if (nlhs != 1) {
		printUsage();
		mexErrMsgTxt("Require one output argument.");
	}
	
	if ((mxGetM(prhs[0]) != 3) || (mxGetN(prhs[0]) % 10 != 0)) {
		printUsage();
		mexErrMsgTxt("One input 3 x 10N matrix is required.");
	}

	int n_instances = mxGetN(prhs[0]) / 10;
	
	double *data = mxGetPr(prhs[0]);

	Eigen::Matrix<double, 3, 10> coeffs;
	Eigen::Matrix<double, 3, 8> solutions;

	if (n_instances == 1) {
		// we are solving a single instance, output is 3xNsols
		coeffs = Eigen::Map<Eigen::Matrix<double, 3, 10>>(data);

		int n_sols = re3q3(coeffs, &solutions);
		plhs[0] = mxCreateDoubleMatrix(3, n_sols, mxREAL);

		Eigen::Map<Eigen::MatrixXd> output_matrix = Eigen::Map<Eigen::MatrixXd>(mxGetPr(plhs[0]), 3, n_sols);
		output_matrix = solutions.block(0, 0, 3, n_sols);
	} else {
		// we are solving multiple instances, output is 3x8N
		plhs[0] = mxCreateDoubleMatrix(3, 8 * n_instances, mxREAL);
		double *output = mxGetPr(plhs[0]);

		Eigen::Map<Eigen::MatrixXd> output_matrix = Eigen::Map<Eigen::MatrixXd>(mxGetPr(plhs[0]), 3, 8*n_instances);
		output_matrix.setZero();

		for (int i = 0; i < n_instances; i++) {
			coeffs = Eigen::Map<Eigen::Matrix<double, 3, 10>>(data + 30 * i);
			solutions.setZero();
			int n_sols = re3q3(coeffs, &solutions);
			output_matrix.block(0, i * 8, 3, n_sols) = solutions.block(0,0,3,n_sols);
		}
	}
}
