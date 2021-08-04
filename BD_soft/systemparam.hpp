#ifndef _SYSTEMPARAM_HPP_
#define _SYSTEMPARAM_HPP_

#include <cmath>

//int N;
//double density;
const int SIMULATION_ID = 0;
//const double L = sqrt(static_cast<double>(N)/density);
const double dt = 0.01;
const double tau = 1.0/dt;
const double Temperature = 0.5;
const double CUTOFF = 3.0;
const double MARGIN = 0.5;
const double SEARCH_LENGTH = (MARGIN + CUTOFF);
const double ML2 = SEARCH_LENGTH * SEARCH_LENGTH;
const double CL2 = (CUTOFF*CUTOFF);
const double RC2 = 1.0 / CL2;
const double RC6 = RC2 * RC2 * RC2;
const double RC12 = RC6 * RC6;
const double C0 = - RC12;
/*
inline void adjust_periodic(double &dx, double &dy){
	const double LH = L * 0.5;
	if (dx < -LH) dx += L;
	if (dx > LH ) dx -= L;
	if (dy < -LH) dy += L;
	if (dy > LH ) dy -= L;
}
*/

#endif