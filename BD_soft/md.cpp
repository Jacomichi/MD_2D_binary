#include <iostream>
#include <assert.h>
#include <cmath>
#include <random>
#include "md.hpp"
#include "systemparam.hpp"
#include "observer.hpp"
#include "variable.hpp"

using namespace std;

random_device rnd;
mt19937 mt(rnd());

MD::MD(int N, double L, double density, double rate, double s1, double s2, int ID) {
	rho = density;
	system_size = L;
	N_particles = N;
	vars = new Variables();
	obs = new Observer();
	mesh = new MeshList(L);
	margin_length = 0.0;
	md_id = ID;
	ratio = rate;
	sigma1 = s1;
	sigma2 = s2;
}

MD::~MD() {
	delete vars;
	delete obs;
	delete mesh;
}

void MD::adjust_periodic(double &dx, double &dy) {
	const double LH = system_size * 0.5;
	if (dx < -LH) dx += system_size;
	if (dx > LH ) dx -= system_size;
	if (dy < -LH) dy += system_size;
	if (dy > LH ) dy -= system_size;
}

void MD::makeconfSquare(void) {
	const double lattice_const = 1.0;
	const int cell_num = static_cast<int>(sqrt(N_particles)) + 1;//端でかぶるのを避けるために1足しておく
	//cout << cell_num << endl;
	for (int iy = 0; iy < cell_num; ++iy) {
		for (int ix = 0; ix < cell_num; ++ix) {
			if (iy * cell_num + ix > N_particles - 1) break;
			vars->add_atoms(ix * lattice_const, iy * lattice_const);
		}
	}
	vars->set_initial_velocity(Temperature);
}

void MD::makeconf(void) {
	const int cell_num = static_cast<int>(sqrt(N_particles)) + 1;
	for (int i = 0; i < N_particles; ++i) {
		vars->add_atoms((system_size - 6.0) * (1.0 / cell_num * (i / cell_num)) + 3, (system_size - 6.0) * (i % cell_num) / cell_num + 3);
	}
	vars->set_initial_velocity(Temperature);
}

void MD::makeRandomconf(void) {
	uniform_real_distribution<double> dist(0, system_size);
	for (int i = 0; i < N_particles; ++i) {
		vars->add_atoms(dist(mt), dist(mt));
	}
	vars->set_initial_velocity(0);
}

void MD::update_position(void) {
	for (auto &a : vars->atoms) {
		a.qx += a.vx * dt;
		a.qy += a.vy * dt;
	}
}

void MD::calculate_force(void) {
	const int p_num = vars->number_of_atoms();
	Atom *atoms = vars->atoms.data();
	for (int i = 0; i < p_num - 1; i++) {
		for (int j = i + 1; j < p_num; j++) {
			double dx = atoms[j].qx - atoms[i].qx;
			double dy = atoms[j].qy - atoms[i].qy;
			adjust_periodic(dx, dy);
			double s = (atoms[i].sigma + atoms[j].sigma) / 2.0;
			double s_sq = s * s;
			double r2 = (dx * dx + dy * dy) / s_sq;
			if (r2 > CL2) continue;
			double r6 = r2 * r2 * r2;
			double df =  - 12.0 / (r6 * r6 * r2) * dt;
			atoms[i].vx += df * dx;
			atoms[i].vy += df * dy;
			atoms[j].vx -= df * dx;
			atoms[j].vy -= df * dy;
		}
	}
}

void MD::periodic(void) {
	for (auto &a : vars->atoms) {
		if (a.qx < 0.0) a.qx += system_size;
		if (a.qy < 0.0) a.qy += system_size;
		if (a.qx > system_size) a.qx -= system_size;
		if (a.qy > system_size) a.qy -= system_size;
		//cout << a.qx << " " << a.qy << endl;
		assert(a.qx < system_size && a.qx > 0.0);
		assert(a.qy < system_size && a.qy > 0.0);
	}
}


void MD::make_pair(void) {
	pairs.clear();
	const int p_num = vars->number_of_atoms();
	//cout << p_num  << endl;
	Atom *atoms = vars->atoms.data();
	for (int i = 0; i < p_num - 1; i++) {
		for (int j = i + 1; j < p_num; j++) {
			double dx = atoms[j].qx - atoms[i].qx;
			double dy = atoms[j].qy - atoms[i].qy;
			adjust_periodic(dx, dy);
			double r2 = (dx * dx + dy * dy);
			if (r2 > ML2) continue;
			Pair p;
			p.i = i;
			p.j = j;
			pairs.push_back(p);
		}
	}
}

void MD::check_pairlist(void) {
	double vmax2 = 0.0;
	for (auto &a : vars->atoms) {
		double v2 = (a.vx * a.vx + a.vy * a.vy);
		if (vmax2 < v2) vmax2 = v2;
	}
	double vmax = sqrt(vmax2);
	margin_length -= vmax * 2.0 * dt;
	if (margin_length < 0.0) {
		margin_length = MARGIN;
		mesh->make_pair(vars, pairs);
	}
}


void MD::calculate_force_soft(void) {
	const int p_pair_num = pairs.size();
	Atom *atoms = vars->atoms.data();
	for (int k = 0; k < p_pair_num; k++) {
		const int i = pairs[k].i;
		const int j = pairs[k].j;
		double dx = atoms[j].qx - atoms[i].qx;
		double dy = atoms[j].qy - atoms[i].qy;
		adjust_periodic(dx, dy);
		double s = (atoms[i].sigma + atoms[j].sigma) / 2.0;
		double s_sq = s * s;
		double r2 = (dx * dx + dy * dy) / s_sq;
		//cout << dx <<  " " << dy << endl;
		if (r2 > CL2) continue;
		double r6 = r2 * r2 * r2;
		double df = - 12.0 / (r6 * r6 * r2) * dt;
		atoms[i].vx += df * dx;
		atoms[i].vy += df * dy;
		atoms[j].vx -= df * dx;
		atoms[j].vy -= df * dy;
		if (abs(df) > 100) {
			cout << df  << endl;
		}
	}
}

void MD::calculate_force_harmonic(void) {
	const int p_pair_num = pairs.size();
	Atom *atoms = vars->atoms.data();
	for (int k = 0; k < p_pair_num; k++) {
		const int i = pairs[k].i;
		const int j = pairs[k].j;
		double dx = atoms[j].qx - atoms[i].qx;
		double dy = atoms[j].qy - atoms[i].qy;
		adjust_periodic(dx, dy);
		double s = (atoms[i].sigma + atoms[j].sigma) / 2.0;
		double s_sq = s * s;
		double r2 = (dx * dx + dy * dy) / s_sq;
		if (r2 > 1.0) {
			continue;
		}
		double r = sqrt(r2);
		double df = - (1.0 - r) / r * dt;
		atoms[i].vx += df * dx;
		atoms[i].vy += df * dy;
		atoms[j].vx -= df * dx;
		atoms[j].vy -= df * dy;
	}
}


void MD::langevin(const double aimed_temperature) {
	const double gamma = 1.0;
	const double T = aimed_temperature;
	const double D = sqrt(2.0 * gamma * T / dt);
	normal_distribution<double> nd(0.0, D);
	for (auto &a : vars->atoms) {
		a.vx += (-gamma * a.vx + nd(mt)) * dt;
		a.vy += (-gamma * a.vy + nd(mt)) * dt;
		//cout << a.vx << "\t" << a.vy << "\t"  << endl;
	}
}

void MD::time_ave_var(Atom ave_var, int time) {
	Atom ave = vars->ave_variable();
	ave_var.vx += ave.vx;
	ave_var.vy += ave.vy;
	ave_var.qx += ave.qx;
	ave_var.qy += ave.qy;
	cout << ave_var.vx / static_cast<double>(time) << " " << ave_var.vy / static_cast<double>(time) << " " << ave_var.qx / static_cast<double>(time) << " " << ave_var.qy / static_cast<double>(time) << "\n";
}

void MD::calculate() {
	//update_position();
	check_pairlist();
	calculate_force_soft();
	langevin(Temperature);
	update_position();
	periodic();
	vars->time += dt;
}

void MD::equillibration(double Temp) {
	//update_position();
	check_pairlist();
	calculate_force_harmonic();
	langevin(Temp);
	update_position();
	periodic();
	vars->time += dt;
}

void MD::run(void) {
	makeRandomconf();
	cout << sigma1 << " " << sigma2 << endl;
	vars->set_radius(ratio, sigma1, sigma2);
	mesh->set_number_of_atoms(vars->number_of_atoms());
	mesh->make_pair(vars, pairs);


	int EQ_STEPS = 500 * tau;
	for (int t = 0; t < EQ_STEPS; t++) {
		if ( (t % 100) == 0) {
			//obs->export_energy_harmonic_dat(vars, pairs, md_id, system_size);
			//vars->export_pos_dat(md_id, N_particles, rho);
		}
		equillibration(1.0 - 1.0 * (1.0 - 1.0 / (t + 1) ));
	}

	cout << "fin create random configuration" << endl;

	EQ_STEPS = 2000 * tau;
	for (int t = 0; t < EQ_STEPS; t++) {
		if ( (t % 100) == 0) {
			//obs->export_energy_harmonic_dat(vars, pairs, md_id, system_size);
			//vars->export_pos_dat(md_id, N_particles, rho);
		}
		calculate();
	}

	cout << "fin equilibration" << endl;

	int STEPS = 100 * tau;
	int OBSERVE_TIMES = 1;
	for (int t = 0; t < STEPS; t++) {
		if ( (t % OBSERVE_TIMES) == 0) {
			//obs->export_energy_dat(vars,pairs, md_id,system_size);
			vars->export_pos_dat(md_id, N_particles, rho, sigma2, 0.01);
		}
		calculate();
	}

	STEPS = 1000 * tau;
	OBSERVE_TIMES = 10;
	for (int t = 0; t < STEPS; t++) {
		if ( (t % OBSERVE_TIMES) == 0) {
			//obs->export_energy_dat(vars,pairs, md_id,system_size);
			vars->export_pos_dat(md_id, N_particles, rho, sigma2, 0.1);
		}
		calculate();
	}

	STEPS = 10000 * tau;
	OBSERVE_TIMES = 100;

	for (int t = 0; t < STEPS; t++) {
		if ( (t % OBSERVE_TIMES) == 0) {
			//obs->export_energy_dat(vars,pairs, md_id,system_size);
			vars->export_pos_dat(md_id, N_particles, rho, sigma2, 1.0);
		}
		calculate();
		//time_ave_var(ave_var,t);
		//cout << ave_var_vx / static_cast<double>(t) << " " << ave_var_vy / static_cast<double>(t) << " " << ave_var_qx / static_cast<double>(t) << " " << ave_var_qy / static_cast<double>(t) << "\n";
	}


}
