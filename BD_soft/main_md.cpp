#include <iostream>
#include <stdio.h>
#include <cmath>
#include "md.hpp"

void simulation(int N, double density, int id) {
	double ratio = 0.5;
	double s1 = 1.0;
	double s2 = 1.4;
	int N1 = ratio * N;
	int N2 = (1.0 - ratio) * N;
	double v1 = M_PI / 4.0 * s1 * s1;
	double v2 = M_PI / 4.0 * s2 * s2;
	double vtot = v1 * N1 + v2 * N2;
	double V = vtot / density;
	double L = sqrt(V);
	MD *md = new MD(N, L, density, ratio, s1, s2, id);
	md->run();
	delete md;
}

int main(int argc, char *argv[]) {
	int N, SIMULATION_ID;
	double density;
	if (argc != 4) {
		cout << "usage:" << argv[0] << " N(int) density(double) ID(int) " << endl;
		exit(EXIT_FAILURE);
	} else {
		char *result;
		N = strtol(argv[1], &result, 0);
		if (*result != '\0' || N <= 0) {
			cout << "The number of particles is not positive int" << endl;
			exit(EXIT_FAILURE);
		}
		density = strtod(argv[2], &result);
		if (*result != '\0' || density < 0) {
			cout << "The density is not positive double" << endl;
			exit(EXIT_FAILURE);
		}

		SIMULATION_ID = strtol(argv[3], &result, 0);
		if (*result != '\0' || SIMULATION_ID < 0) {
			cout << "The SIMULATION_ID should be positive." << endl;
			exit(EXIT_FAILURE);
		}
	}
	simulation(N, density, SIMULATION_ID);
}
