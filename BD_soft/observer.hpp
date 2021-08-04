#ifndef _OBSERVER_HPP_
#define _OBSERVER_HPP_

#include "variable.hpp"

class Observer {
public:
	double kinetic_energy(Variables *vars);
	double potential_energy(Variables *vars, vector<Pair> &pairs,double L);
	double potential_energy(Variables *vars,double L);
	double potential_energy_harmonic(Variables *vars, vector<Pair> &pairs, double L);
	double temperature(Variables *vars) {
		return kinetic_energy(vars);
	}
	double total_energy(Variables *vars, vector<Pair> &pairs,double L) {
		return kinetic_energy(vars) + potential_energy(vars, pairs,L);
	}
	void export_energy_dat(Variables *vars, vector<Pair> &pairs, int id,double L);
	void export_energy_harmonic_dat(Variables *vars, vector<Pair> &pairs, int id, double L);
	void export_energy_dat(Variables *vars, int id,double L);
private:
	void adjust_periodic(double &dx, double &dy, double L);
};

#endif