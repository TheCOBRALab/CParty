#ifndef RESULT_HEADER
#define RESULT_HEADER
#include "base_types.hh"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

class Result {
  public:
    // constructor
    Result(std::string sequence, std::string restricted, double restricted_energy, std::string final_structure, double final_energy,
           std::string final_structure_pf, pf_t pf_energy, std::string MEA_structure, pf_t MEA, std::string centroid_structure, pf_t distance, pf_t frequency);
    // destructor
    ~Result();

    // getter
    std::string get_sequence();
    std::string get_restricted();
    std::string get_final_structure();
    double get_restricted_energy();
    double get_final_energy();
    std::string get_final_structure_pf();
    pf_t get_pf_energy();
    std::string get_MEA_structure();
    pf_t get_MEA();
    std::string get_centroid_structure();
    pf_t get_distance();
    pf_t get_frequency();

    struct Result_comp {
        bool operator()(Result &x, Result &y) const {
            if (x.get_final_energy() != y.get_final_energy()) return x.get_final_energy() < y.get_final_energy();
            return x.get_restricted_energy() < y.get_restricted_energy();
        }
    } result_comp;

  private:
    std::string sequence;
    std::string restricted;
    double restricted_energy;
    std::string final_structure;
    double final_energy;
    std::string final_structure_pf;
    pf_t pf_energy;
    std::string MEA_structure;
    pf_t MEA;
    std::string centroid_structure;
    pf_t distance;
    pf_t frequency;
};

#endif
