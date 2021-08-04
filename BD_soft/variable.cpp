#include <iostream>
#include <fstream>
#include <random>
#include "systemparam.hpp"
#include "variable.hpp"

using namespace std;

void Variables::add_atoms(double x, double y) {
	Atom a;
	a.qx = x;
	a.qy = y;
	a.vx = 0;
	a.vy = 0;
	a.sigma = 0;
	atoms.push_back(a);
}

void Variables::set_initial_velocity(const double V0) {
	mt19937 mt(2);
	uniform_real_distribution<double> ud(0.0, 1.0);
	double ave_vx = 0.0;
	double ave_vy = 0.0;
	for (auto &a : atoms) {
		//向きはランダムで、速さは与えた大きさV0にする
		double phi = 2.0 * ud(mt) * M_PI;
		double vx = V0 * cos(phi);
		double vy = V0 * sin(phi);
		a.vx = vx;
		a.vy = vy;
		ave_vx += vx;
		ave_vy += vy;
	}
	const int p_num = atoms.size();
	//粒子数で割る
	ave_vx /= static_cast<double>(p_num);
	ave_vy /= static_cast<double>(p_num);
	for (auto &a : atoms) {
		a.vx -= ave_vx;
		a.vy -= ave_vy;
	}
}

void Variables::set_radius(double ratio,double s1,double s2) {
	int N = number_of_atoms();
	int N_1 = N * ratio;
	int i = 0;
	for (auto &a : atoms) {
		if (i < N_1) {
			a.sigma = s1;
		} else {
			a.sigma = s2;
		}
		++i;
	}
}


int Variables::number_of_atoms(void) {
	return static_cast<int> (atoms.size());
}

void Variables::export_cdview(void) {
	static int count = 0;
	char filename[256];
	sprintf(filename, "conf%03d.cdv", count);
	++count;
	ofstream ofs(filename);
	int i = 0;
	for (auto &a : atoms) {
		ofs << i << " " << "0" << " ";
		ofs << a.qx << " ";
		ofs << a.qy << "\n";
		++i;
	}
}

void Variables::export_trajectory_csv(int id, int N, double density) {
	char filename[256];
	sprintf(filename, "N%05drho%4.2fT%4.2f_%02d_soft.csv", N, density, Temperature, id);
	ofstream ofs(filename, ios_base::app);
	int i = 0;
	for (auto &a : atoms) {
		ofs << i << "," << a.qx << "," << a.qy << "," << a.vx << "," << a.vy << "\n";
		++i;
	}
}

void Variables::export_pos_dat(int id, int N, double density,double step,double s2) {
	char filename[256];
	sprintf(filename, "binary_pos_N%05drho%4.2fT%4.2f_%4.2ftau_s%4.2f_%02d_soft.dat", N, density, Temperature,step,s2, id);
	ofstream ofs(filename, ios_base::app);
	int i = 0;
	for (auto &a : atoms) {
		ofs << a.qx << " " << a.qy << "\n";
		++i;
	}
}

Atom Variables::ave_variable(void) {
	Atom ave;
	ave.vx = 0.0;
	ave.vy = 0.0;
	ave.qx = 0.0;
	ave.qy = 0.0;
	int p_num = atoms.size();
	p_num  = static_cast<double>(p_num);
	for (auto &a : atoms) {
		ave.vx += a.vx / p_num;
		ave.vy += a.vy / p_num;
		ave.qx += a.qx / p_num;
		ave.qy += a.qy / p_num;
	}
	//cout << ave.vx << " " << ave.vy << " " << ave.qx << " " << ave.qy << "\n";
	return ave;
}