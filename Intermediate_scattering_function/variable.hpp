#ifndef _VARIABLE_HPP_
#define _VARIABLE_HPP_

#include <vector>

using namespace std;

struct Atom {
	double qx, qy;
};

class Variables {
public:
	vector<Atom> atoms;
	void add_atoms(double x,double y);
	int number_of_atoms(void);
};

#endif