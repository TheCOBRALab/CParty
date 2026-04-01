#include "SHAPE.hh"
#include "part_func.hh"
#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>


SHAPEData::SHAPEData(const std::string &filename, cand_pos_t n, double slope, double intercept){
    this->slope = slope;
    this->intercept = intercept;
    this->n = n;
    this->calculated.resize(n+1,0);
    this->expcalculated.resize(n+1,1);
    if(exists(filename)){
        std::string line;
        std::ifstream in(filename);
        std::getline(in,line);
        in.close();
        std::istringstream ss(line);
        std::string token;
        std::string name;
        cand_pos_t length;
        std::getline(ss, name, '\t');
        std::getline(ss, token, '\t');
        length = std::stoi(token);
        std::getline(ss, token, '\t');
        if(length>n){
            std::cerr << "ERROR: Number of entries in file does not match length of sequence" << std::endl;
            exit(EXIT_FAILURE);
        }
        cand_pos_t i = 0;
        while (std::getline(ss, token, '\t')) {
            ++i;
            if(token == "NULL") continue;
            double reactivity = std::stod(token);
            calculated[i] = 100*calculate(reactivity);
            expcalculated[i] = calculated[i];
        }
    }

}

bool SHAPEData::exists(const std::string &filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

double SHAPEData::get_calculated(cand_pos_t i){
    return calculated[i];
}

double SHAPEData::get_expcalculated(cand_pos_t i){
    return expcalculated[i];
}

double SHAPEData::calculate(double reactivity){
    return slope*std::log(reactivity+1) + intercept;
}

void SHAPEData::rescale_calculate(double kT, double TT, int smooth){
    int pf_smooth = smooth;
    for(cand_pos_t i = 1; i<=n;++i){
        if(this->expcalculated[i] != 1){
            this->expcalculated[i] = RESCALE_BF(this->expcalculated[i], this->expcalculated[i]*3, TT, kT);
        }
    }
}