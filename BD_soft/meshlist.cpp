#include <iostream>
#include <assert.h>
#include <algorithm>
#include "systemparam.hpp"
#include "meshlist.hpp"

using namespace std;

MeshList::MeshList(double L) {
	system_size =L;
	mesh_per_axis = static_cast<int>(system_size / SEARCH_LENGTH) - 1;
	mesh_size = system_size / mesh_per_axis;
	assert(mesh_per_axis > 2);
	assert(mesh_size > SEARCH_LENGTH);
	num_of_mesh = mesh_per_axis * mesh_per_axis;
	count.resize(num_of_mesh);
	indexes.resize(num_of_mesh);
}

void MeshList::adjust_periodic(double &dx, double &dy) {
	const double LH = system_size * 0.5;
	if (dx < -LH) dx += system_size;
	if (dx > LH ) dx -= system_size;
	if (dy < -LH) dy += system_size;
	if (dy > LH ) dy -= system_size;
}

void MeshList::make_pair(Variables *vars, vector<Pair> &pairs) {
	pairs.clear();
	const int p_num = vars->number_of_atoms();
	//cout << p_num  << endl;
	Atom *atoms = vars->atoms.data();
	vector<int> particle_position(p_num);//粒子iのいるメッシュを管理するベクトル
	vector<int> pointer(p_num);//その番地に何こ粒子がいるかを管理するベクトル
	fill(particle_position.begin(), particle_position.end(), 0);
	fill(count.begin(), count.end(), 0);
	fill(pointer.begin(), pointer.end(), 0);

	double im = 1.0 / mesh_size;
	for (int i = 0; i < p_num; i++) {
		int ix = static_cast<int>(atoms[i].qx * im);
		int iy = static_cast<int>(atoms[i].qy * im);
		if (ix < 0) ix += mesh_per_axis;
		if (ix >= mesh_per_axis) ix -= mesh_per_axis;
		if (iy < 0) iy += mesh_per_axis;
		if (iy >= mesh_per_axis) iy -= mesh_per_axis;

		int index = ix + iy * mesh_per_axis;
		//cout << index << " " << ix << " " << iy << " " << p_num <<endl;
		assert(index >= 0);
		assert(index < num_of_mesh);
		count[index]++;//number of particles in the mesh
		particle_position[i] = index; //which meshes the paticles are in.
	}

	indexes[0] = 0;
	int sum = 0;
	for (int i = 0; i < num_of_mesh - 1 ; i++) {
		sum += count[i];
		indexes[i + 1] = sum;//the index of the first particle in the mesh
	}
	for (int i = 0; i < p_num; i++) {
		int pos = particle_position[i];
		int j = indexes[pos] + pointer[pos];
		sorted_buffer[j] = i;
		++pointer[pos];
	}
	for (int i = 0; i < num_of_mesh; i++) {
		search(i, vars, pairs);
	}
}

void MeshList::search_other(int id, int ix, int iy, Variables *vars, vector<Pair> &pairs) {
	if (ix < 0) ix += mesh_per_axis;
	if (ix >= mesh_per_axis) ix -= mesh_per_axis;
	if (iy < 0) iy += mesh_per_axis;
	if (iy >= mesh_per_axis) iy -= mesh_per_axis;
	int id2 = ix + iy * mesh_per_axis;
	Atom *atoms = vars->atoms.data();
	for (int k = indexes[id]; k < indexes[id] + count[id]; k++) {
		for (int l = indexes[id2]; l < indexes[id2] + count[id2]; l++) {
			int i = sorted_buffer[k];
			int j = sorted_buffer[l];
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


void MeshList::search(int id, Variables *vars, vector<Pair> &pairs) {
	int ix = id % mesh_per_axis;
	int iy = (id / mesh_per_axis) % mesh_per_axis;
	search_other(id, ix + 1, iy    ,vars, pairs);
	search_other(id, ix - 1, iy + 1,vars, pairs);
	search_other(id, ix    , iy + 1,vars, pairs);
	search_other(id, ix + 1, iy + 1,vars, pairs);

	int si = indexes[id];
	int n = count[id];
	Atom *atoms = vars->atoms.data();
	for (int k = si; k < si + n - 1; k++) {
		for (int l = k + 1; l < si + n; l++) {
			int i = sorted_buffer[k];
			int j = sorted_buffer[l];
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
