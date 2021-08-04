#ifndef _STATIC_STRUCTURE_
#define _STATIC_STRUCTURE_

#include "variable.hpp"
#include <vector>
#include <iostream>
#include <string>

using namespace std;

class Fkt{
private:
	Variables *vars;
	vector<string> split(string& input_line, char delimiter);
	void calc_fskt_dat(int N,double density,string filename);
	void calc_fskt_simple_dat(int N,double density,string filename);
	void calc_fskt_simple_fraction_dat(int N,double density,string filename);
	void calc_fkt_simple_dat(int N,double density,string filename);
	void calc_fkt_simple_fraction_dat(int N,double density,string filename);
	void calc_fkt_dat(int N,double density,string filename);
	void get_para(string filename, int &N, double &density);
	void adjust_periodic(double &dx, double &dy, double L);

public:
	Fkt();
	~Fkt();
	void calc(string filename);
};

#endif
