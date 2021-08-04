#ifndef _VARIABLE_HPP_
#define _VARIABLE_HPP_

#include <vector>

using namespace std;

struct Atom {
	double qx, qy;
	double vx, vy;
	double sigma;
};

struct Pair {
	int i, j;
};

class Variables {
public:
	vector<Atom> atoms;
	double time;
	Variables(void) {time = 0.0;}
	void add_atoms(double x, double y);
	int number_of_atoms(void);
	void set_initial_velocity(const double);
	void set_radius(double ratio, double s1, double s2);
	void export_cdview(void);
	void export_trajectory_csv(int id, int N, double density);
	void export_pos_dat(int id, int N, double density, double step,double s2);
	Atom ave_variable(void);
};

#endif