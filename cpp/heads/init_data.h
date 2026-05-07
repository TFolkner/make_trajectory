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

// ошибки и статусы
enum cfg_read_states {
	// ошибки
	dont_opend_file,
	// статусы
	cfg_read_success
};

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
	std::string& get_author();
	std::string& get_project_name();
	std::string& get_version();

	std::string& get_name_start();
	double& get_phi_start();
	double& get_lambda_start();

	std::string& get_name_end();
	double& get_phi_end();
	double& get_lambda_end();

	double& get_W();
	double& get_t();
	double& get_h();
	double& get_Hz();

	double& get_W_Earth();
	double& get_R_Earth();

};


#endif