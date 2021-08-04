#include "static_structure.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

int main(int argc, char *argv[]) {
    string filename;
    FILE *fp;
    if (argc != 2) {
        cout << "usage:" << argv[0] << " filename(string) " << endl;
        exit(EXIT_FAILURE);
    } else {
        filename = string(argv[1]);
        fp = fopen( filename.c_str(), "r" );
        if (fp == NULL) {
            cout << "The file can not be opened." << endl;
            return -1;
        }
    }
    fclose(fp);
    Sk *sk = new Sk();
    sk->calc(filename);
    delete sk;
}