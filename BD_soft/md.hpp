#ifndef _MD_HPP_
#define _MD_HPP_

#include "variable.hpp"
#include "observer.hpp"
#include "meshlist.hpp"

class MD {
private:
	Variables *vars;
	Observer *obs;
	MeshList *mesh;
	vector<Pair> pairs;
	double system_size;
	int N_particles;
	double rho;
	double margin_length;
	double ratio;
	double sigma1;
	double sigma2;
	int md_id;
	void adjust_periodic(double &dx, double &dy);
	void makeconfSquare(void);
	void makeconf(void);
	void makeRandomconf(void);
	void update_position(void);
	void calculate_force(void);
	void calculate_force_soft(void);
	void calculate_force_harmonic(void);
	void periodic(void);
	void calculate(void);
	void equillibration(double Temp);
	void make_pair(void);
	void check_pairlist(void);
	void langevin(const double aimed_temperature);
	void time_ave_var(Atom ave_var, int time);

public:
	MD(int N, double L, double density,double rate,double s1,double s2, int ID);
	~MD();
	void run(void);
};

#endif