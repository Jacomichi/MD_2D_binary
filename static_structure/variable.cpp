#include <iostream>
#include <fstream>
#include "variable.hpp"

using namespace std;

void Variables::add_atoms(int id,double x, double y, double vx, double vy) {
	Atom a;
	a.id = id;
	a.qx = x;
	a.qy = y;
	a.vx = vx;
	a.vy = vy;
	atoms.push_back(a);
}


int Variables::number_of_atoms(void) {
	return static_cast<int> (atoms.size());
}
