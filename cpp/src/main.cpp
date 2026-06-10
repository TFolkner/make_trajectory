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



int main() {
	switch (dispatcher()) {
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

