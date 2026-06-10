//heads
#include <iostream>
#include <fstream>
#include <math.h>

#include "../heads/init_data.h"




class main_dispatcher {
private:

public:
	int operator()() {
		init_data init{};
		init.make_header();
		init.disp_config();
		return 0;
	}
};








