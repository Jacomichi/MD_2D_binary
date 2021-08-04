#pragma once

#include <vector>
#include "variable.hpp"

using namespace std;

class MeshList {
private:
	double mesh_size;
	int mesh_per_axis;
	int num_of_mesh;
	double system_size;
	vector<int> count;//どの番地に何個の粒子がいるか
	vector<int> indexes;//番地番号にソートした時の、その番地の一番最初になる粒子のインデックス、n番地の住人がはいる一番最初の場所
	vector<int> sorted_buffer;//番地番号でソートした原子のインデックス、
	void search(int id, Variables *vars, vector<Pair> &pairs);
	void search_other(int id, int ix, int iy, Variables *vars, vector<Pair> &pairs);
	void adjust_periodic(double &dx, double &dy);
public:
	MeshList(double L);
	void make_pair(Variables *vars,vector<Pair> &pairs);
	void set_number_of_atoms(int p_num){
		sorted_buffer.resize(p_num);
	}
};