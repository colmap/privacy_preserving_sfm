#include <Eigen/Dense>
#include <re3q3.h>
#include <iostream>
#include <time.h>

#define TEST(FUNC) if(!FUNC()) { std::cout << #FUNC"\033[1m\033[31m FAILED!\033[0m\n"; } else { std::cout << #FUNC"\033[1m\033[32m PASSED!\033[0m\n"; passed++;} num_tests++; 
#define REQUIRE(COND) if(!(COND)) { std::cout << "Failure: "#COND" was not satisfied.\n"; return false; }


bool verify_solutions(const Eigen::Matrix<double, 3, 10> &coeffs,
					const Eigen::Matrix<double, 3, 8> &solution,
					int n_sol, double tol) {

	Eigen::Matrix<double, 10, 1> mons;

	bool ok = true;

	for (int i = 0; i < n_sol; i++) {
		double x = solution(0, i);
		double y = solution(1, i);
		double z = solution(2, i);

		mons << x * x, x * y, x * z, y * y, y * z, z * z, x, y, z, 1.0;

		Eigen::Matrix<double, 3, 1> residuals = coeffs * mons;
		//std::cout << "solution: " << i << " res: " << residuals << "\n";
		ok &= (residuals.cwiseAbs().maxCoeff() < tol);
	}
	return  ok;
}


bool test_random_coefficients() {
	Eigen::Matrix<double, 3, 10> coeffs;
	Eigen::Matrix<double, 3, 8> solutions;

	coeffs.setRandom();

	//std::cout << "coeffs: " << coeffs << "\n";

	int n_sols = re3q3(coeffs, &solutions);
	
	return verify_solutions(coeffs, solutions, n_sols, 1e-8);
}


bool test_degenerate_for_x() {
	Eigen::Matrix<double, 3, 10> coeffs;
	Eigen::Matrix<double, 3, 8> solutions;

	coeffs.setRandom();
	coeffs.col(3) = 0.5 * (coeffs.col(5) + coeffs.col(4));

	int n_sols = re3q3(coeffs, &solutions);

	return verify_solutions(coeffs, solutions, n_sols, 1e-8);
}


bool test_degenerate_for_y() {
	Eigen::Matrix<double, 3, 10> coeffs;
	Eigen::Matrix<double, 3, 8> solutions;

	coeffs.setRandom();
	coeffs.col(0) = 0.5 * (coeffs.col(5) + coeffs.col(2));

	int n_sols = re3q3(coeffs, &solutions);

	return verify_solutions(coeffs, solutions, n_sols, 1e-8);
}

bool test_degenerate_for_z() {
	Eigen::Matrix<double, 3, 10> coeffs;
	Eigen::Matrix<double, 3, 8> solutions;

	coeffs.setRandom();
	coeffs.col(0) = 0.5 * (coeffs.col(1) + coeffs.col(3));

	int n_sols = re3q3(coeffs, &solutions);

	return verify_solutions(coeffs, solutions, n_sols, 1e-8);
}

bool test_degenerate_for_xy() {
	Eigen::Matrix<double, 3, 10> coeffs;
	Eigen::Matrix<double, 3, 8> solutions;

	coeffs.setRandom();
	coeffs.col(0) = 0.5 * (coeffs.col(5) + coeffs.col(2));
	coeffs.col(3) = 0.5 * (coeffs.col(5) + coeffs.col(4));

	int n_sols = re3q3(coeffs, &solutions);	

	return verify_solutions(coeffs, solutions, n_sols, 1e-8);
}


bool test_pure_squares() {
	Eigen::Matrix<double, 3, 10> coeffs;
	Eigen::Matrix<double, 3, 8> solutions;

	coeffs.setZero();
	coeffs(0, 0) = 1.0;
	coeffs(0, 9) = -1.0;
	coeffs(1, 3) = 1.0;
	coeffs(1, 9) = -1.0;
	coeffs(2, 5) = 1.0;
	coeffs(2, 9) = -1.0;

	int ok = 0;
	for(int test = 0; test < 1000; ++test) {
		int n_sols = re3q3(coeffs, &solutions);

		REQUIRE(n_sols == 8);

		if(verify_solutions(coeffs, solutions, n_sols, 1e-8))
			++ok;
	}

	return ok > 950; // require 95% successrate at least
}


int main() {
	
	unsigned int seed = (unsigned int)time(0);	
	srand(seed);

	std::cout << "Running tests... (seed = " << seed << ")\n\n";

	int passed = 0;
	int num_tests = 0;
	
	TEST(test_random_coefficients);
	TEST(test_degenerate_for_x);
	TEST(test_degenerate_for_y);
	TEST(test_degenerate_for_z);
	TEST(test_degenerate_for_xy);
	TEST(test_pure_squares);

	std::cout << "\nDone! Passed " << passed << "/" << num_tests << " tests.\n";
}
