#include <iostream>
#include <fstream>
#include "systemparam.hpp"
#include "observer.hpp"

using namespace std;

void Observer::adjust_periodic(double &dx, double &dy, double L) {
	const double LH = L * 0.5;
	if (dx < -LH) dx += L;
	if (dx > LH ) dx -= L;
	if (dy < -LH) dy += L;
	if (dy > LH ) dy -= L;
}

double Observer::kinetic_energy(Variables *vars) {
	double kinetic = 0;
	for (auto &a : vars->atoms) {
		kinetic += a.vx * a.vx + a.vy * a.vy;
	}
	kinetic /= static_cast<double>(vars->number_of_atoms());
	return kinetic * 0.5;
}

double Observer::potential_energy(Variables *vars, vector<Pair> &pairs, double L) {
	double potential = 0.0;
	const int p_pair_num = pairs.size();
	const int p_num = vars->number_of_atoms();
	Atom *atoms = vars->atoms.data(); //配列の先頭のポインタを返す
	for (int k = 0 ; k < p_pair_num; k++) {
		const int i = pairs[k].i;
		const int j = pairs[k].j;
		double dx = atoms[j].qx - atoms[i].qx;
		double dy = atoms[j].qy - atoms[i].qy;
		adjust_periodic(dx, dy, L);
		double s = (atoms[i].sigma + atoms[j].sigma) / 2.0;
		double s_sq = s * s;
		double r2 = (dx * dx + dy * dy) / s_sq;
		if (r2 > CL2) {
			continue;
		}
		double r6 = r2 * r2 * r2;
		double r12 = r6 * r6;
		potential += (1.0 / r12) + C0;
	}

	potential /= static_cast<double>(p_num);//大きさがほぼ0の無視できるのを消しただけだから総数は0
	return potential;
}

double Observer::potential_energy_harmonic(Variables *vars, vector<Pair> &pairs, double L) {
	double potential = 0.0;
	const int p_pair_num = pairs.size();
	const int p_num = vars->number_of_atoms();
	Atom *atoms = vars->atoms.data();
	for (int k = 0 ; k < p_pair_num; k++) {
		const int i = pairs[k].i;
		const int j = pairs[k].j;
		double dx = atoms[j].qx - atoms[i].qx;
		double dy = atoms[j].qy - atoms[i].qy;
		adjust_periodic(dx, dy, L);
		double s = (atoms[i].sigma + atoms[j].sigma) / 2.0;
		double s_sq = s * s;
		double r2 = (dx * dx + dy * dy) / s_sq;
		if (r2 > 1.0) {
			continue;
		}
		double r = sqrt(r2);
		potential += (1.0-r)*(1.0-r);
	}

	potential /= static_cast<double>(p_num);
	return potential;
}

double Observer::potential_energy(Variables *vars, double L) {
	double potential = 0.0;
	const int p_num = vars->number_of_atoms();
	Atom *atoms = vars->atoms.data(); //配列の先頭のポインタを返す
	for (int i = 0 ; i < p_num - 1; i++) {
		for (int j = i + 1; j < p_num; j++) {
			double dx = atoms[j].qx - atoms[i].qx;
			double dy = atoms[j].qy - atoms[i].qy;
			adjust_periodic(dx, dy, L);
			double s = (atoms[i].sigma + atoms[j].sigma) / 2.0;
			double s_sq = s * s;
			double r2 = (dx * dx + dy * dy)/s_sq;
			if (r2 > CL2) {
				continue;
			}
			double r6 = r2 * r2 * r2;
			double r12 = r6 * r6;
			potential += (1.0 / r12 ) + C0;
		}
	}

	potential /= static_cast<double>(p_num);//大きさがほぼ0の無視できるのを消しただけだから総数は0
	return potential;
}

void Observer::export_energy_dat(Variables *vars, vector<Pair> &pairs, int id, double L) {
	char filename[256];
	sprintf(filename, "energy%02d.dat", id);
	ofstream ofs(filename, ios_base::app);
	double kinetic = kinetic_energy(vars);
	double potential = potential_energy(vars, pairs, L);
	ofs << vars->time << " " << kinetic << " " << potential << " " << kinetic + potential << " " << temperature(vars) << "\n";
}

void Observer::export_energy_harmonic_dat(Variables *vars, vector<Pair> &pairs, int id, double L) {
	char filename[256];
	sprintf(filename, "energy%02d.dat", id);
	ofstream ofs(filename, ios_base::app);
	double kinetic = kinetic_energy(vars);
	double potential = potential_energy_harmonic(vars, pairs, L);
	ofs << vars->time << " " << kinetic << " " << potential << " " << kinetic + potential << " " << temperature(vars) << "\n";
}

void Observer::export_energy_dat(Variables *vars, int id, double L) {
	char filename[256];
	sprintf(filename, "energy%02d.dat", id);
	ofstream ofs(filename, ios_base::app);
	double kinetic = kinetic_energy(vars);
	double potential = potential_energy(vars, L);
	ofs << vars->time << " " << kinetic << " " << potential << " " << kinetic + potential << " " << temperature(vars) << "\n";
}
