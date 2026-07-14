/*
* Описание проекта
*/

#include "../heads/ortodromy.h"
#include "../heads/init_data.h"


int dispatcher() {
	init_data init{};
	init.make_header();
	init.disp_config();
	
	return 0;
}

int test() {
	Eigen::Matrix<double, 3, 3> test_mtr = Eigen::Matrix<double, 3, 3>::Random();
	test_mtr = test_mtr + test_mtr;

	for (int row = 0; row < test_mtr.rows(); row++) {
		for (int cols = 0; cols < test_mtr.cols(); cols++) {
			std::cout << test_mtr(row, cols) << " ";
		}
		std::cout << "\n";
	}
	return 0;
}

int main() {
	switch (test()) {
	case 0:
		break;

	case 1:
		break;

	default:
		break;
	}

	system ("pause");
	return 0;
}

