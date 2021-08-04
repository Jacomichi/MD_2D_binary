#include <iostream>
#include <fstream>
#include "variable.hpp"

using namespace std;

void Variables::add_atoms(double x, double y) {
	Atom a;
	a.qx = x;
	a.qy = y;
	atoms.push_back(a);
}


int Variables::number_of_atoms(void) {
	return static_cast<int> (atoms.size());
}
