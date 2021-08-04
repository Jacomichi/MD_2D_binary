#ifndef _VARIABLE_HPP_
#define _VARIABLE_HPP_

#include <vector>

using namespace std;

struct DataFrame {
	double qx, qy;
	double vx, vy;
	int id;
};

class Variables {
public:
	vector<DataFrame> dataframe;
	void add_DataFrame(int id,double x,double y,double vx,double vy);
	int number_of_data(void);
};

#endif