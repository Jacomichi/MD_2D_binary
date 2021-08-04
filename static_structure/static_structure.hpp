#ifndef _STATIC_STRUCTURE_
#define _STATIC_STRUCTURE_

#include "variable.hpp"
#include <vector>
#include <iostream>
#include <string>

using namespace std;

class Sk{
private:
	Variables *vars;
	vector<string> split(string& input_line, char delimiter);
	void get_csv(string filename);
	int p_num(void);
	void calc_sk(int N,double density,string filename);
	void calc_sk_dat(int N,double density,string filename);
	void calc_sk_fraction_dat(int N,double density,string filename);
	void get_para(string filename, int &N, double &density);

public:
	Sk();
	~Sk();
	void calc(string filename);
};

#endif
