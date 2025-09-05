#include "utils/filesystem.hpp"
#include <iostream>

#ifndef TEST_DIRECTORY
#define TEST_DIRECTORY "./"
#endif

std::string cfd::test_directory_file(std::string path){
	return TEST_DIRECTORY + path;
}

std::string cfd::tmp_directory_file(std::string path){
	return TEST_DIRECTORY + std::string("../tmp_data/") + path;
}
