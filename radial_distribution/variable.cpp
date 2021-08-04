#include <iostream>
#include <fstream>
#include "variable.hpp"

using namespace std;

void Variables::add_DataFrame(int id,double x, double y, double vx, double vy) {
	DataFrame a;
	a.id = id;
	a.qx = x;
	a.qy = y;
	a.vx = vx;
	a.vy = vy;
	dataframe.push_back(a);
}


int Variables::number_of_data(void) {
	return static_cast<int> (dataframe.size());
}
