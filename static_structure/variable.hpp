#ifndef _VARIABLE_HPP_
#define _VARIABLE_HPP_

#include <vector>

using namespace std;

struct Atom {
	double qx, qy;
	double vx, vy;
	int id;
};

class Variables {
public:
	vector<Atom> atoms;
	void add_atoms(int id,double x,double y,double vx,double vy);
	int number_of_atoms(void);
};

#endif