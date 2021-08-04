#ifndef _STATIC_STRUCTURE_
#define _STATIC_STRUCTURE_

#include "variable.hpp"
#include <vector>
#include <string>

using namespace std;

class Gr{
private:
	Variables *vars;
	vector<string> split(string& input_line, char delimiter);
	void get_csv(string filename);
	int p_num(void);
	void calc_gr(int N,double density, string filename);
	void calc_gr_fraction(int N,double density, string filename);
	void adjust_periodic(double &dx, double &dy,double L);
	void get_para(string filename, int &N, double &density);

public:
	Gr();
	~Gr();
	void calc(string filename);
};

#endif
