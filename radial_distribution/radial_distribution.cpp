#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <regex>
#include <string>
#include <assert.h>
#include "radial_distribution.hpp"

using namespace std;

Gr::Gr() {
    vars = new Variables();
}

Gr::~Gr() {
    delete vars;
}

void Gr::get_para(string filename, int &N, double &density) {
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
}

vector<string> Gr::split(string& input_line, char delimiter)
{
    istringstream file_input(input_line);
    string line_in_file;
    vector<string> result;
    while (getline(file_input, line_in_file, delimiter)) {
        result.push_back(line_in_file);
    }
    return result;
}

void Gr::get_csv(string filename) {
    ifstream ifs(filename);
    string line;
    while (getline(ifs, line)) {
        vector<string> strvec = split(line, ',');
        for (int i = 0; i < static_cast<int>(strvec.size() / 5); i++) {
            int p_id = stoi(strvec.at(5 * i));
            double x = stof(strvec.at(5 * i + 1));
            double y = stof(strvec.at(5 * i + 2));
            double vx = stof(strvec.at(5 * i + 3));
            double vy = stof(strvec.at(5 * i + 4));
            vars->add_DataFrame(p_id, x, y, vx, vy);
        }
    }
}


int Gr::p_num(void) {
    DataFrame *dataframe = vars->dataframe.data();
    int max_num = 0;
    for (int i = 0; i < vars->number_of_data(); ++i) {
        max_num = max(dataframe[i].id, max_num);
    }
    return max_num + 1;
}

void Gr::adjust_periodic(double &dx, double &dy, double L) {
    const double LH = L * 0.5;
    if (dx < -LH) dx += L;
    if (dx > LH ) dx -= L;
    if (dy < -LH) dy += L;
    if (dy > LH ) dy -= L;
}

void Gr::calc_gr(int N, double density, string filename) {
    double L = sqrt(static_cast<double>(N) / density);
    double dr = 0.01;
    int count = (int)(L / 2.0 / dr);//MDセルの半分しか、
    unsigned long n[count];
    double r[count];

    for (int i = 0; i < count; ++i) {
        n[i] = 0;
        r[i] = i * dr;
    }

    DataFrame *dataframe = vars->dataframe.data();
    int T_count = vars->number_of_data() / N;
    for (int i = 0; i < T_count; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                double dx = dataframe[k + i * N].qx - dataframe[j + i * N].qx;
                double dy = dataframe[k + i * N].qy - dataframe[j + i * N].qy;
                adjust_periodic(dx, dy, L);
                //cout << dx << "\t" << dy << "\t" << dz << endl;
                double r2 = sqrt(dx * dx + dy * dy);
                if (r2 < r[count - 1] && j != k) {
                    n[(int)(r2 / dr)]++;
                }
            }
        }
    }

    filename = "gr_" + filename;
    ofstream ofs(filename);

    for (int i = 0; i < count; ++i) {
        double para = 2.0 * M_PI * density * r[i] * N * dr * T_count;
        ofs << r[i] << " " << n[i] / para << "\n";
        //fprintf(f, "%f\t%f\n", r[i], n[i] / para);
    }
}

void Gr::calc_gr_fraction(int N, double density, string filename) {
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
    double dr = 0.01;
    int count = (int)(L / 2.0 / dr);//MDセルの半分しか、
    unsigned long n[count];
    double r[count];
    DataFrame input_data[N];
    for (int i = 0; i < N ; ++i) {
        input_data[N].qx = 0.0;
        input_data[N].qy = 0.0;
    }

    for (int i = 0; i < count; ++i) {
        n[i] = 0;
        r[i] = i * dr;
    }

    ifstream input_file(filename);
    if (input_file.fail()) {
        cout << "入力ファイルをオープンできません" << endl;
        assert(!input_file.fail());
    }

    int Tcount = 0;
    while (!input_file.eof()) {
        for (int n = 0; n < N; ++n) {
            input_file >> input_data[n].qx >> input_data[n].qy;
        }
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                double dx = input_data[k].qx - input_data[j].qx;
                double dy = input_data[k].qy - input_data[j].qy;
                adjust_periodic(dx, dy, L);
                //cout << dx << "\t" << dy << endl;
                double r2 = sqrt(dx * dx + dy * dy);
                if (r2 < r[count - 1] && j != k) {
                    n[(int)(r2 / dr)]++;
                }
            }
        }
        //cout << count << endl;
        Tcount++;
    }

    /*
    DataFrame *dataframe = vars->dataframe.data();
    int T_count = vars->number_of_data() / N;
    for (int i = 0; i < T_count; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                double dx = dataframe[k + i * N].qx - dataframe[j + i * N].qx;
                double dy = dataframe[k + i * N].qy - dataframe[j + i * N].qy;
                adjust_periodic(dx, dy, L);
                cout << dx << "\t" << dy << endl;
                double r2 = sqrt(dx * dx + dy * dy);
                if (r2 < r[count - 1] && j != k) {
                    n[(int)(r2 / dr)]++;
                }
            }
        }
    }
    */

    filename = "gr_" + filename;
    ofstream ofs(filename);

    for (int i = 0; i < count; ++i) {
        double para = 2.0 * M_PI * density * r[i] * N * dr * Tcount;
        ofs << r[i] << " " << n[i] / para << "\n";
        //fprintf(f, "%f\t%f\n", r[i], n[i] / para);
    }
}


void Gr::calc(string filename) {
    int N;
    double density;
    get_para(filename, N, density);
    //cout << N << " " << density << endl;
    //get_csv(filename);
    calc_gr(N, density, filename);
    //calc_gr_fraction(N, density, filename);
}



