#include "../heads/init_data.h"


// getters ==========================================================
std::string& init_data::get_author() { return author; }
std::string& init_data::get_project_name() { return project_name; }
std::string& init_data::get_version() { return version; }
std::string& init_data::get_name_start() { return name_start; }
double& init_data::get_phi_start() { return phi_start; }
double& init_data::get_lambda_start() { return lambda_start; }
std::string& init_data::get_name_end() { return name_end; }
double& init_data::get_phi_end() { return phi_end; }
double& init_data::get_lambda_end() { return lambda_end; }
double& init_data::get_W() { return W; }
double& init_data::get_t() { return t; }
double& init_data::get_h() { return h; }
double& init_data::get_Hz() { return Hz; }
double& init_data::get_W_Earth() { return W_Earth; }
double& init_data::get_R_Earth() { return R_Earth; }

init_data::init_data() {
	try {
		toml::table init_file = toml::parse_file("config.toml");

		author = init_file["author"].value_or("Default author");
		project_name = init_file["project_name"].value_or("Default project");
		version = init_file["version"].value_or("Default version");

		name_start = init_file["name_start"].value_or("Default start nme");
		phi_start = init_file["phi_start"].value_or(999.);
		lambda_start = init_file["lambda_start"].value_or(999.);

		name_end = init_file["name_end"].value_or("Default end name");
		phi_end = init_file["phi_end"].value_or(999.);
		lambda_end = init_file["lambda_end"].value_or(999.);

		W = init_file["W"].value_or(999.);
		t = init_file["t"].value_or(999.);
		h = init_file["h"].value_or(999.);
		Hz = init_file["Hz"].value_or(999.);

		W_Earth = init_file["W_Earth"].value_or(999.);
		R_Earth = init_file["R_Earth"].value_or(999.);
	}
	catch (const toml::parse_error& err) {
		std::cerr << "TOML parse error:\n"<< err << '\n';
	}
}


int init_data::make_header() {
	std::cout << "==========================" << std::endl;
	std::cout << "Author - " << get_author() << std::endl;
	std::cout << "Project - " << get_project_name() << std::endl;
	std::cout << "Version - " << get_version() << std::endl;
	std::cout << std::endl;
	return 0;
}


int init_data::disp_config() {
	std::cout << "==========================" << std::endl;
	std::cout << "Start -> " << "\"" << get_name_start() << "\"" << " " << get_phi_start() << ", " << get_lambda_start() << std::endl;
	std::cout << "End -> " << "\"" << get_name_end() << "\"" << " " << get_phi_end() << ", " << get_lambda_end() << std::endl;
	std::cout << W << std::endl;
	std::cout << t << std::endl;
	std::cout << h << std::endl;
	std::cout << Hz << std::endl;
	std::cout << W_Earth << std::endl;
	std::cout << R_Earth << std::endl;
	std::cout << std::endl;
	return 0;
}

