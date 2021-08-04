#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <complex>
#include <regex>
#include <assert.h>
#include "static_structure.hpp"

using namespace std;

Sk::Sk() {
    vars = new Variables();
}

Sk::~Sk() {
    delete vars;
}

void Sk::get_para(string filename, int &N, double &density) {
    //This function gets this type filename "N01000rho0.50T1.00_01_soft.csv";
    vector<double> parameter;
    smatch match;
    regex rx2(R"(([0-9.]+))");
    regex_search(filename, match, rx2);
    auto it = filename.cbegin();
    while (regex_search(it, filename.cend(), match, rx2)) {
        // i==0の全体マッチはスキップする
        for (int i = 1; i < static_cast<int>(match.size()); ++i) {
            //printf(match.str(i).c_str());
            //printf("\n%d\n", i);
            char *result;
            parameter.push_back(strtod(match.str(i).c_str(), &result));
        }
        // マッチした位置+長さ分イテレータをずらす
        it += match.position(0) + match.length(0);
    }
    N = static_cast<int>(parameter[0]);
    density = parameter[1];
    cout << "N:" << N << " rho:" << density << endl;
}

vector<string> Sk::split(string& input_line, char delimiter)
{
    istringstream file_input(input_line);
    string line_in_file;
    vector<string> result;
    while (getline(file_input, line_in_file, delimiter)) {
        result.push_back(line_in_file);
    }
    return result;
}


void Sk::get_csv(string filename) {
    ifstream ifs(filename);
    string line;
    while (getline(ifs, line)) {
        vector<string> strvec = split(line, ',');
        for (int i = 0; i < static_cast<int>(strvec.size() / 3); i++) {
            int p_id = stoi(strvec.at(5 * i));
            double x = stof(strvec.at(5 * i + 1));
            double y = stof(strvec.at(5 * i + 2));
            //double vx = stof(strvec.at(5 * i + 3));
            //double vy = stof(strvec.at(5 * i + 4));
            vars->add_atoms(p_id, x, y, 0, 0/*, vx, vy*/);
            //printf("%5d\t%5f\t%5f\t%5f\t%5f\n", i, atoms[i].qx,atoms[i].qy,atoms[i].vx,atoms[i].vy);
        }
    }
}


int Sk::p_num(void) {
    Atom *atoms = vars->atoms.data();
    int max_num = 0;
    for (int i = 0; i < vars->number_of_atoms(); ++i) {
        max_num = max(atoms[i].id, max_num);
    }
    return max_num + 1;
}

void Sk::calc_sk(int N , double density, string filename) {
    double L = sqrt(static_cast<double>(N) / density);
    double k_min = 2.0 * M_PI / static_cast<double>(L);
    double k_max = M_PI;
    int kcount = static_cast<int>(k_max / k_min);
    double S[kcount][kcount];
    for (int nx = 0; nx < kcount; ++nx) {
        for (int ny = 0; ny < kcount; ++ny) {
            S[nx][ny] = 0.0;
        }
    }

    Atom *atoms = vars->atoms.data();
    int T_count = vars->number_of_atoms() / N;
    for (int i = 0; i < T_count; i++) {
        for (int nx = 0; nx < kcount; nx++) {
            double kx = k_min * static_cast<double>(nx);
            for (int ny = 0; ny < kcount; ny++) {
                double a = 0.0;
                double b = 0.0;
                double ky = k_min * static_cast<double>(ny);
                for (int j = 0; j < N; ++j) {
                    a += cos(atoms[j + i * N].qx * kx + atoms[j + i * N].qy * ky);
                    b += sin(atoms[j + i * N].qx * kx + atoms[j + i * N].qy * ky);
                }
                S[nx][ny] += (a * a + b * b) / static_cast<double>(N * T_count);
            }
        }
    }

    double dk = 0.01;
    int k_MAXCount = (int)(sqrt(2.0) * k_max / dk);
    //cout << k_MAXCount << endl;
    double s[k_MAXCount];
    long int counter[k_MAXCount];

    for (int i = 0; i < k_MAXCount; ++i) {
        s[i] = 0;
        counter[i] = 0;
    }


    for (int nx = 0; nx < kcount; ++nx) {
        for (int ny = 0; ny < kcount; ++ny) {
            int distance = (int)(sqrt(nx * nx + ny * ny) * k_min / dk);
            if (distance < k_MAXCount) {
                s[distance] += S[nx][ny];
                counter[distance] += 1;
            }
        }
    }

    /*
    for (int i = 1; i < k_MAXCount; ++i) {
        if (counter[i] != 0) {
            cout << static_cast<double>(i)*dk << " " << s[i] / static_cast<double>(counter[i]) << "\n";
        }
    }
    */

    filename = "Sk_" + filename;
    ofstream ofs(filename, ios_base::app);
    for (int i = 1; i < k_MAXCount; ++i) {
        if (counter[i] != 0) {
            ofs << static_cast<double>(i)*dk << " " << s[i] / static_cast<double>(counter[i]) << "\n";
        }
    }
}

void Sk::calc_sk_dat(int N , double density, string filename) {
    double L = sqrt(static_cast<double>(N) / density);
    const double k_min = 2.0 * M_PI / static_cast<double>(L);
    const double k_max = 2.0 * M_PI + 2.0;
    int kcount = static_cast<int>(k_max / k_min);
    double S[kcount][kcount];
    Atom input_data[N];
    for (int i = 0; i < N ; ++i) {
        input_data[N].qx = 0.0;
        input_data[N].qy = 0.0;
    }
    for (int nx = 0; nx < kcount; ++nx) {
        for (int ny = 0; ny < kcount; ++ny) {
            S[nx][ny] = 0.0;
        }
    }

    ifstream input_file(filename);
    if (input_file.fail()) {
        cout << "入力ファイルをオープンできません" << endl;
        assert(!input_file.fail());
    }

    int count = 0;

    while (!input_file.eof()) {
        for (int n = 0; n < N; ++n) {
            input_file >> input_data[n].qx >> input_data[n].qy;
        }
        //cout << count << endl;
        for (int nx = 0; nx < kcount; nx++) {
            double kx = k_min * static_cast<double>(nx);
            for (int ny = 0; ny < kcount; ny++) {
                double a = 0.0;
                double b = 0.0;
                double ky = k_min * static_cast<double>(ny);
                //cout << kx << " " << ky << endl;
                for (int j = 0; j < N; ++j) {
                    a += cos(input_data[j].qx * kx + input_data[j].qy * ky);
                    b += sin(input_data[j].qx * kx + input_data[j].qy * ky);
                }
                S[nx][ny] += (a * a + b * b) ;
            }
        }

        count++;
    }
    cout << count << endl;

    double dk = k_min;
    int k_MAXCount = (int)(sqrt(2.0) * k_max / dk);
    //cout << k_MAXCount << endl;
    double s[k_MAXCount];
    long int counter[k_MAXCount];

    for (int i = 0; i < k_MAXCount; ++i) {
        s[i] = 0;
        counter[i] = 0;
    }

    for (int nx = 0; nx < kcount; ++nx) {
        for (int ny = 0; ny < kcount; ++ny) {
            int distance = (int)(sqrt(nx * nx + ny * ny) * k_min / dk);
            if (distance < k_MAXCount) {
                s[distance] += S[nx][ny] / static_cast<double>(N * count);
                counter[distance] += 1;
            }
        }
    }


    filename = "Sk_" + filename;
    ofstream ofs(filename);
    for (int i = 1; i < k_MAXCount; ++i) {
        if (counter[i] != 0) {
            ofs << static_cast<double>(i)*dk << " " << s[i] / static_cast<double>(counter[i]) << "\n";
        }
    }
}

void Sk::calc_sk_fraction_dat(int N , double density, string filename) {
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
    const double k_min = 2.0 * M_PI / static_cast<double>(L);
    const double k_max = 2.0 * M_PI + 2.0;
    int kcount = static_cast<int>(k_max / k_min);
    double S[kcount][kcount];
    Atom input_data[N];
    for (int i = 0; i < N ; ++i) {
        input_data[N].qx = 0.0;
        input_data[N].qy = 0.0;
    }
    for (int nx = 0; nx < kcount; ++nx) {
        for (int ny = 0; ny < kcount; ++ny) {
            S[nx][ny] = 0.0;
        }
    }

    ifstream input_file(filename);
    if (input_file.fail()) {
        cout << "入力ファイルをオープンできません" << endl;
        assert(!input_file.fail());
    }

    int count = 0;

    while (!input_file.eof()) {
        for (int n = 0; n < N; ++n) {
            input_file >> input_data[n].qx >> input_data[n].qy;
        }
        //cout << count << endl;
        for (int nx = 0; nx < kcount; nx++) {
            double kx = k_min * static_cast<double>(nx);
            for (int ny = 0; ny < kcount; ny++) {
                double a = 0.0;
                double b = 0.0;
                double ky = k_min * static_cast<double>(ny);
                //cout << kx << " " << ky << endl;
                for (int j = 0; j < N; ++j) {
                    a += cos(input_data[j].qx * kx + input_data[j].qy * ky);
                    b += sin(input_data[j].qx * kx + input_data[j].qy * ky);
                }
                S[nx][ny] += (a * a + b * b) ;
            }
        }

        count++;
    }
    cout << count << endl;

    double dk = k_min;
    int k_MAXCount = (int)(sqrt(2.0) * k_max / dk);
    //cout << k_MAXCount << endl;
    double s[k_MAXCount];
    long int counter[k_MAXCount];

    for (int i = 0; i < k_MAXCount; ++i) {
        s[i] = 0;
        counter[i] = 0;
    }

    for (int nx = 0; nx < kcount; ++nx) {
        for (int ny = 0; ny < kcount; ++ny) {
            int distance = (int)(sqrt(nx * nx + ny * ny) * k_min / dk);
            if (distance < k_MAXCount) {
                s[distance] += S[nx][ny] / static_cast<double>(N * count);
                counter[distance] += 1;
            }
        }
    }


    filename = "Sk_" + filename;
    ofstream ofs(filename);
    for (int i = 1; i < k_MAXCount; ++i) {
        if (counter[i] != 0) {
            ofs << static_cast<double>(i)*dk << " " << s[i] / static_cast<double>(counter[i]) << "\n";
        }
    }
}


void Sk::calc(string filename) {
    int N;
    double density;
    get_para(filename, N, density);
    calc_sk_fraction_dat(N, density, filename);
}

