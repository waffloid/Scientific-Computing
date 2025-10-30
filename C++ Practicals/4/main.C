#include<iostream>
#include<vector>
#include<array>

using Matrix3x3 = std::array<std::array<double, 3>, 3>;
using Matrix2x2 = std::array<std::array<double, 2>, 2>;

double determinant_2x2(Matrix2x2 matrix) {
	return matrix.at(1).at(1) * matrix.at(0).at(0) - matrix.at(1).at(0) * matrix.at(0).at(1);
}

void fprint_matrix_2x2(Matrix2x2 matrix) {
	for (int i = 0; i<2; ++i) {
		std::cout << matrix.at(i).at(0) << " " << matrix.at(i).at(1) << std::endl;
	}
}

double determinant_3x3(const Matrix3x3& matrix) {
	double sum = 0;
	int sign = 1;

	for (int i = 0; i < 3; ++i) {
		Matrix2x2 sub_matrix;

		bool seen_q = false;

		for (int q = 0; q < 3; q++) {
			if (q == i) {
				continue;
			}

			for (int p = 1; p < 3; p++) {
				int output_p = p - 1;
				int output_q = seen_q ? 1 : 0;

				sub_matrix.at(output_p).at(output_q) = matrix.at(p).at(q);
			}

			seen_q = true;
		}
		double sub_det = determinant_2x2(sub_matrix);

		sum += sign * matrix.at(0).at(i) * sub_det;
		sign *= -1;
	}

	return sum;
}

int main() {
	const Matrix3x3 my_mat = {{
		{1.0, 2.0, 1.0},
		{1.0, 2.0, 0.0},
		{1.0, 2.0, 3.0}
	}};

	std::cout << determinant_3x3(my_mat) << std::endl;

	return 0;
}