#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <complex>
#include <regex>
#include <assert.h>
#include "Intermediate_scattering_function.hpp"

using namespace std;

Fkt::Fkt() {
    vars = new Variables();
}

Fkt::~Fkt() {
    delete vars;
}

void Fkt::get_para(string filename, int &N, double &density) {
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

vector<string> Fkt::split(string& input_line, char delimiter)
{
    istringstream file_input(input_line);
    string line_in_file;
    vector<string> result;
    while (getline(file_input, line_in_file, delimiter)) {
        result.push_back(line_in_file);
    }
    return result;
}


void Fkt::calc_fskt_dat(int N , double density, string filename) {
    double L = sqrt(static_cast<double>(N) / density);
    double k = 2 * M_PI;
    double steps = 0.01;
    double t_MAX = 10.0 / steps;

    string filename_out = "Fskt_k=" + to_string(k) + "_" + filename;
    ofstream ofs(filename_out);
    if (ofs.fail()) {
        cout << "出力ファイルをオープンできません" << endl;
        assert(!ofs.fail());
    }

    for (int t = 1; t <= t_MAX; t++) {
        ifstream input_file(filename);
        if (input_file.fail()) {
            cout << "入力ファイルをオープンできません" << endl;
            assert(!input_file.fail());
        }

        Atom buf[N][t];

        for (int j = 0; j < t; ++j) {
            for (int i = 0; i < N; ++i) {
                input_file >> buf[i][j].qx >> buf[i][j].qy;
                //cout << buf[i][j].qx << " " << buf[i][j].qy << endl;
            }
        }
        cout << t << endl;

        double a = 0.;
        int sampleNum = 50;
        for (int s = 0; s < 50; ++s) {
            int t1 = s % t;
            int t2 = (s + t - 1) % t;
            //cout << t1 << " " << t2 << endl;
            for (int i = 0; i < N; ++i) {
                double dx = buf[i][t1].qx - buf[i][t2].qx;
                double dy = buf[i][t1].qy - buf[i][t2].qy;
                //cout << buf[i][t1].qx << " " << buf[i][t1].qy << endl;
                adjust_periodic(dx, dy, L);
                a += cos(k * dx) + cos(k * dy);
                input_file >> buf[i][t1].qx >> buf[i][t1].qy;
            }
        }
        a = a / 2.0 / static_cast<double>(sampleNum) / N;
        ofs << t - 1 << " " << a << "\n";
    }
}

void Fkt::calc_fskt_simple_dat(int N , double density, string filename) {
    double L = sqrt(static_cast<double>(N) / density);
    double k = 2 * M_PI;
    double steps = 0.01;
    double Max_time = 10.0;
    int t_MAX = static_cast<int>(Max_time / steps);

    double a[t_MAX];
    for (int n = 0; n < t_MAX; ++n) {
        a[n] = 0.0;
    }

    string filename_out = "Fskt_simp_k=" + to_string(k) + "_" + filename;
    ofstream ofs(filename_out);
    if (ofs.fail()) {
        cout << "出力ファイルをオープンできません" << endl;
        assert(!ofs.fail());
    }

    int Nsample = 50;
    for (int t_ave = 0; t_ave < Nsample; ++t_ave) {
        ifstream input_file(filename);
        if (input_file.fail()) {
            cout << "入力ファイルをオープンできません" << endl;
            assert(!input_file.fail());
        }

        Atom buf[N][t_MAX];

        //throw trash away the data which we don't use
        int counter = 0;
        double trash;
        for (int n = 0; n < counter; ++n) {
            for (int j = 0; j < t_MAX; ++j) {
                for (int i = 0; i < N; ++i) {
                    input_file >> trash >> trash;
                }
            }
        }

        //input sample Block
        for (int j = 0; j < t_MAX; ++j) {
            for (int i = 0; i < N; ++i) {
                input_file >> buf[i][j].qx >> buf[i][j].qy;
            }
        }
        counter++;

        for (int t = 0; t < t_MAX; ++t) {
            for (int i = 0; i < N; ++i) {
                double dx = buf[i][t].qx - buf[i][0].qx;
                double dy = buf[i][t].qy - buf[i][0].qy;
                //cout << buf[i][t1].qx << " " << buf[i][t1].qy << endl;
                adjust_periodic(dx, dy, L);
                a[t] += (cos(k * dx) + cos(k * dy)) / static_cast<double>(N);
            }
        }
    }

    for (int t = 0; t < t_MAX; ++t) {
        ofs << t * steps << " " << a[t] / static_cast<double>(Nsample) / 2.0 << "\n";
    }
}

void Fkt::calc_fskt_simple_fraction_dat(int N , double density, string filename) {
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
    double k = 2 * M_PI;
    double steps = 0.1;
    double Max_time = 10.0;
    int t_MAX = static_cast<int>(Max_time / steps);

    double a[t_MAX];
    for (int n = 0; n < t_MAX; ++n) {
        a[n] = 0.0;
    }

    string filename_out = "Fskt_simp_k=" + to_string(k) + "_" + filename;
    ofstream ofs(filename_out);
    if (ofs.fail()) {
        cout << "出力ファイルをオープンできません" << endl;
        assert(!ofs.fail());
    }

    int Nsample = 100;
    for (int t_ave = 0; t_ave < Nsample; ++t_ave) {
        ifstream input_file(filename);
        if (input_file.fail()) {
            cout << "入力ファイルをオープンできません" << endl;
            assert(!input_file.fail());
        }

        Atom buf[N][t_MAX];

        //throw trash away the data which we don't use
        int counter = 0;
        double trash;
        for (int n = 0; n < counter; ++n) {
            for (int j = 0; j < t_MAX; ++j) {
                for (int i = 0; i < N; ++i) {
                    input_file >> trash >> trash;
                }
            }
        }

        //input sample Block
        for (int j = 0; j < t_MAX; ++j) {
            for (int i = 0; i < N; ++i) {
                input_file >> buf[i][j].qx >> buf[i][j].qy;
            }
        }
        counter++;

        for (int t = 0; t < t_MAX; ++t) {
            for (int i = 0; i < N; ++i) {
                double dx = buf[i][t].qx - buf[i][0].qx;
                double dy = buf[i][t].qy - buf[i][0].qy;
                adjust_periodic(dx, dy, L);
                a[t] += (cos(k * dx) + cos(k * dy)) / static_cast<double>(N);
            }
        }
    }

    for (int t = 0; t < t_MAX; ++t) {
        ofs << t * steps << " " << a[t] / static_cast<double>(Nsample) / 2.0 << "\n";
    }
}

void Fkt::calc_fkt_simple_dat(int N , double density, string filename) {
    double L = sqrt(static_cast<double>(N) / density);
    double k = 2 * M_PI;
    double steps = 0.01;
    double Max_time = 10.0;
    int t_MAX = static_cast<int>(Max_time / steps);

    double a[t_MAX];
    double b[t_MAX];
    double c[t_MAX];
    double ax[t_MAX];
    double bx[t_MAX];
    double ay[t_MAX];
    double by[t_MAX];
    for (int n = 0; n < t_MAX; ++n) {
        a[n] = 0.0;
        b[n] = 0.0;
        c[n] = 0.0;
    }

    string filename_out = "Fkt_k=" + to_string(k) + "_" + filename;
    ofstream ofs(filename_out);
    if (ofs.fail()) {
        cout << "出力ファイルをオープンできません" << endl;
        assert(!ofs.fail());
    }

    int Nsample = 50;
    for (int t_ave = 0; t_ave < Nsample; ++t_ave) {
        ifstream input_file(filename);
        if (input_file.fail()) {
            cout << "入力ファイルをオープンできません" << endl;
            assert(!input_file.fail());
        }

        Atom buf[N][t_MAX];
        for (int n = 0; n < t_MAX; ++n) {
            ax[n] = 0.0;
            ay[n] = 0.0;
            bx[n] = 0.0;
            by[n] = 0.0;
        }

        //throw trash away the data which we don't use
        int counter = 0;
        double trash;
        for (int n = 0; n < counter; ++n) {
            for (int j = 0; j < t_MAX; ++j) {
                for (int i = 0; i < N; ++i) {
                    input_file >> trash >> trash;
                }
            }
        }

        //input sample Block
        for (int j = 0; j < t_MAX; ++j) {
            for (int i = 0; i < N; ++i) {
                input_file >> buf[i][j].qx >> buf[i][j].qy;
                ax[j] += cos(k * buf[i][j].qx);
                ay[j] += cos(k * buf[i][j].qy);
                bx[j] += sin(k * buf[i][j].qx);
                by[j] += sin(k * buf[i][j].qy);
            }
        }
        counter++;

        for (int t = 0; t < t_MAX; ++t) {
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    double dx = buf[i][t].qx - buf[j][0].qx;
                    double dy = buf[i][t].qy - buf[j][0].qy;
                    double dx_nobound = dx;
                    double dy_nobound = dy;
                    //cout << buf[i][t1].qx << " " << buf[i][t1].qy << endl;
                    adjust_periodic(dx, dy, L);
                    if (i == j) {
                        a[t] += (cos(k * dx) + cos(k * dy)) / static_cast<double>(N) / 2.0; //i=jだけ足す
                    }
                    //b[t] += (cos(k * dx) + cos(k * dy)) / static_cast<double>(N) / 2.0; //i,j両方たす
                    b[t] += (cos(k * dx_nobound) + cos(k * dy_nobound)) / static_cast<double>(N) / 2.0;
                    c[t] += (ax[0] * ax[t] + ay[0] * ay[t] + bx[0] * bx[t] + by[0] * by[t]) / 2.0 / static_cast<double>(N);
                    //c[t] += (exp(-imagnum * k * dx) + exp(-imagnum * k * dy)) / static_cast<double>(N) / 2.0;
                }
            }
        }
    }
    for (int t = 0; t < t_MAX; ++t) {
        ofs << t * steps << " " << a[t] / static_cast<double>(Nsample)  << " " << b[t] / static_cast<double>(Nsample) << " " <<  c[t] / static_cast<double>(Nsample) << "\n";
    }

}

void Fkt::calc_fkt_simple_fraction_dat(int N , double density, string filename) {
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
    double k = 2 * M_PI / 1.2;
    double steps = 0.01;
    double Max_time = 1.0;
    int t_MAX = static_cast<int>(Max_time / steps);

    double a[t_MAX];
    double b[t_MAX];
    double c[t_MAX];
    double ax[t_MAX];
    double bx[t_MAX];
    double ay[t_MAX];
    double by[t_MAX];
    for (int n = 0; n < t_MAX; ++n) {
        a[n] = 0.0;
        b[n] = 0.0;
        c[n] = 0.0;
    }

    string filename_out = "Fkt_simp_k=" + to_string(k) + "_" + filename;
    ofstream ofs(filename_out);
    if (ofs.fail()) {
        cout << "出力ファイルをオープンできません" << endl;
        assert(!ofs.fail());
    }

    int Nsample = 100;
    for (int t_ave = 0; t_ave < Nsample; ++t_ave) {
        ifstream input_file(filename);
        if (input_file.fail()) {
            cout << "入力ファイルをオープンできません" << endl;
            assert(!input_file.fail());
        }

        Atom buf[N][t_MAX];
        for (int n = 0; n < t_MAX; ++n) {
            ax[n] = 0.0;
            ay[n] = 0.0;
            bx[n] = 0.0;
            by[n] = 0.0;
        }

        //throw trash away the data which we don't use
        int counter = 0;
        double trash;
        for (int n = 0; n < counter; ++n) {
            for (int j = 0; j < t_MAX; ++j) {
                for (int i = 0; i < N; ++i) {
                    input_file >> trash >> trash;
                }
            }
        }

        //input sample Block
        for (int j = 0; j < t_MAX; ++j) {
            for (int i = 0; i < N; ++i) {
                input_file >> buf[i][j].qx >> buf[i][j].qy;
                ax[j] += cos(k * buf[i][j].qx);
                ay[j] += cos(k * buf[i][j].qy);
                bx[j] += sin(k * buf[i][j].qx);
                by[j] += sin(k * buf[i][j].qy);
            }
        }
        counter++;

        for (int t = 0; t < t_MAX; ++t) {
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    double dx = buf[i][t].qx - buf[j][0].qx;
                    double dy = buf[i][t].qy - buf[j][0].qy;
                    double dx_nobound = dx;
                    double dy_nobound = dy;
                    //cout << buf[i][t1].qx << " " << buf[i][t1].qy << endl;
                    adjust_periodic(dx, dy, L);
                    if (i == j) {
                        a[t] += (cos(k * dx) + cos(k * dy)) / static_cast<double>(N) / 2.0; //i=jだけ足す
                    }
                    //b[t] += (cos(k * dx) + cos(k * dy)) / static_cast<double>(N) / 2.0; //i,j両方たす
                    b[t] += (cos(k * dx_nobound) + cos(k * dy_nobound)) / static_cast<double>(N) / 2.0;
                    c[t] += (ax[0] * ax[t] + ay[0] * ay[t] + bx[0] * bx[t] + by[0] * by[t]) / 2.0 / static_cast<double>(N);
                    //c[t] += (exp(-imagnum * k * dx) + exp(-imagnum * k * dy)) / static_cast<double>(N) / 2.0;
                }
            }
        }
    }
    for (int t = 0; t < t_MAX; ++t) {
        ofs << t * steps << " " << a[t] / static_cast<double>(Nsample)  << " " << b[t] / static_cast<double>(Nsample) << " " <<  c[t] / static_cast<double>(Nsample) << "\n";
    }

}

void Fkt::adjust_periodic(double &dx, double &dy, double L) {
    const double LH = L * 0.5;
    if (dx < -LH) dx += L;
    if (dx > LH ) dx -= L;
    if (dy < -LH) dy += L;
    if (dy > LH ) dy -= L;
}

void Fkt::calc_fkt_dat(int N , double density, string filename) {
    //This is still Fs(k,t)
    double L = sqrt(static_cast<double>(N) / density);
    double k = 2 * M_PI;
    double t_MAX = 50;

    string filename_out = "Fkt_k=" + to_string(k) + "_" + filename;
    ofstream ofs(filename_out);
    if (ofs.fail()) {
        cout << "出力ファイルをオープンできません" << endl;
        assert(!ofs.fail());
    }

    for (int t = 1; t <= t_MAX; t++) {
        ifstream input_file(filename);
        if (input_file.fail()) {
            cout << "入力ファイルをオープンできません" << endl;
            assert(!input_file.fail());
        }

        Atom buf[N][t];

        for (int j = 0; j < t; ++j) {
            for (int i = 0; i < N; ++i) {
                input_file >> buf[i][j].qx >> buf[i][j].qy;
                //cout << buf[i][j].qx << " " << buf[i][j].qy << endl;
            }
        }
        cout << t << endl;

        double a = 0.;
        double b = 0.;
        for (int s = 0; s < 10; ++s) {
            int t1 = s % t;
            int t2 = (s + t - 1) % t;
            //cout << t1 << " " << t2 << endl;
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    double dx = buf[i][t1].qx - buf[j][t2].qx;
                    double dy = buf[i][t1].qy - buf[j][t2].qy;
                    adjust_periodic(dx, dy, L);
                    if (i == j) {
                        a += cos(k * dx) + cos(k * dy);
                    } else {
                        b += cos(k * dx) + cos(k * dy);
                    }
                    input_file >> buf[j][t2].qx >> buf[j][t2].qy;
                }
            }
        }
        a = a / 2.0 / 10.0 / N;
        b = b / 2.0 / 10.0 / N;
        ofs << t - 1 << " " << a << " " << b << "\n";
    }
}


void Fkt::calc(string filename) {
    int N;
    double density;
    get_para(filename, N, density);
    //get_csv(filename);
    //calc_sk(N, density, filename);
    //calc_fskt_dat(N, density, filename);
    calc_fkt_simple_dat(N, density, filename);
    //calc_fkt_simple_fraction_dat(N, density, filename);
}

