/*
* Класс чтения конфигурационного файла
*
* Примечание:
* -
*/

#ifndef INIT_DATA
#define INIT_DATA

#include <iostream>
#include <fstream>
#include <string>

#include "../libs/tomlplusplus-master/toml.hpp"

#define init_file_name_ "config.cfg"


class init_data {
private:
	std::string 
		author,
		project_name,
		version
		;
	std::string
		name_start, 
		name_end
		;
	double
		phi_start, lambda_start, // [deg]
		phi_end, lambda_end, // [deg]
		W, // vehicle linear speed [m/s]
		t, // travel time [s]
		h, // travel height [m]
		Hz, // model frequency [^-1]
		W_Earth, // Earth angular spead [1/sec] 
		R_Earth // Earth radius [m]
		;

public:
	init_data();

	int make_header();
	int disp_config();


	// getters
	const std::string& get_author();
	const std::string& get_project_name();
	const std::string& get_version();

	const std::string& get_name_start();
	const double& get_phi_start();
	const double& get_lambda_start();

	const std::string& get_name_end();
	const double& get_phi_end();
	const double& get_lambda_end();

	const double& get_W();
	const double& get_t();
	const double& get_h();
	const double& get_Hz();

	const double& get_W_Earth();
	const double& get_R_Earth();

};


#endif